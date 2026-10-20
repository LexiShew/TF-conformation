#!/bin/bash
# apply_2026-10-15_fixes.sh — regenerate everything invalidated by the
# 2026-10-15 fixes, in dependency order, on the cluster.
#
# Context: docs/codebase_review_2026-10-15.md (findings F1/F2/F6/F8/F9/F12/F13/F14),
#          docs/tbp_label_error_downstream.md, docs/retired_pilots.md
#
# The fixes are already in the code. What is still stale is the DATA on disk:
#   * 190 tbp npz carry NDT80's Y_pwm                      -> rebuild from stage 4
#   * 264 checkpoints were trained unseeded                -> retrain to get a paired design
#   * 22 id_benchmark_*.json hold normal-approximation p   -> re-evaluate
#   * mechanism CSVs hold non-monotonic BH q-values        -> re-run the analysis suite
#
# USAGE
#   ./scripts/pipeline/apply_2026-10-15_fixes.sh --help
#   ./scripts/pipeline/apply_2026-10-15_fixes.sh preflight          # free, read-only
#   ./scripts/pipeline/apply_2026-10-15_fixes.sh analysis           # ~minutes, CPU, no SLURM
#   ./scripts/pipeline/apply_2026-10-15_fixes.sh tbp      --go      # SLURM: stages 4-7
#   ./scripts/pipeline/apply_2026-10-15_fixes.sh ndt80    --go      # SLURM: stages 1-7
#   ./scripts/pipeline/apply_2026-10-15_fixes.sh reeval   --go      # SLURM: stage 7 only
#   ./scripts/pipeline/apply_2026-10-15_fixes.sh retrain  --go      # SLURM: stages 6-7, EXPENSIVE
#
# Phases that submit jobs are DRY-RUN by default and need --go to actually
# submit. Everything prints what it would do first.
#
# ORDER MATTERS. 'retrain' uses run_multiseed_pilot.sh, which requires stages
# 1-5 to already exist for that pilot, so:
#
#   preflight          always first; free
#   analysis           independent of everything else; safe any time
#   tbp     --go       MUST precede retrain for tbp   (rebuilds its stage 4-5)
#   ndt80   --go       MUST precede retrain for ndt80 (it has no stages 1-5 yet)
#   reeval  --go       optional shortcut: fixes the p-value arithmetic on the
#                      EXISTING unseeded checkpoints. Cheap, but it does not fix
#                      the pairing (F6). Skip it if you are going straight to
#                      retrain.
#   retrain --go       last; supersedes reeval for every pilot it covers
#
# Minimum path to trustworthy paired numbers:
#   preflight -> tbp --go -> ndt80 --go -> (wait) -> retrain --go -> analysis

set -eo pipefail

TFCONF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${TFCONF_DIR}"

PHASE="${1:-}"
GO=0
for a in "$@"; do [ "$a" = "--go" ] && GO=1; done
# SEEDS is a LIST; its length is what is passed to
# run_multiseed_pilot.sh, which always uses seeds 1..N.
SEEDS="${SEEDS:-1 2 3 4 5}"
PILOTS_ALL="csl egr1 engrailed err ets1 foxa hsf lef1 ndt80 nfat runx tbp"

say()  { printf '\n\033[1m== %s\033[0m\n' "$*"; }
note() { printf '   %s\n' "$*"; }
run()  {
    if [ "$GO" = "1" ]; then
        printf '   $ %s\n' "$*"
        eval "$@"
    else
        printf '   [dry-run] %s\n' "$*"
    fi
}

usage() {
    sed -n '2,32p' "$0" | sed 's/^# \{0,1\}//'
    exit 0
}
[ -z "${PHASE}" ] && usage
[ "${PHASE}" = "--help" ] || [ "${PHASE}" = "-h" ] && usage

# ---------------------------------------------------------------------------
case "${PHASE}" in

# ===========================================================================
preflight)
# ===========================================================================
    say "PHASE 0 — preflight (read-only, free)"

    note "0a. generated configs match config/pilots.tsv"
    python3 scripts/pipeline/gen_pilot_configs.py --check

    note "0b. all 12 pilots pass validation (C1-C10)"
    python3 scripts/pipeline/validate_pilot.py --all --quiet

    note "0c. tbp label is corrected"
    tbp_label=$(bash -c 'source config/pilots/tbp.sh >/dev/null 2>&1; echo $PWM_LABEL')
    if [ "${tbp_label}" = "TBP_HUMAN.H11MO.0.A" ]; then
        note "    OK: PWM_LABEL=${tbp_label}"
    else
        echo "    FAIL: tbp PWM_LABEL is '${tbp_label}', expected TBP_HUMAN.H11MO.0.A" >&2
        exit 1
    fi

    note "0d. contaminated tbp npz still on disk (these are what must be rebuilt)"
    for d in output/stage4_npz/tbp/output output/stage4_npz_dnarelax/tbp/output; do
        n=$(ls "${d}"/*MA0343* 2>/dev/null | wc -l)
        note "    ${d}: ${n} npz carrying NDT80's label"
    done

    note "0e. retired pilots fully gone from live config"
    if ls config/pilots/dux4*.sh config/pilots/irf*.sh >/dev/null 2>&1; then
        echo "    FAIL: retired pilot configs still present" >&2; exit 1
    fi
    note "    OK: no dux4/irf configs"

    note "0f. queue view of what is complete vs outstanding"
    python3 scripts/pipeline/pilot_queue.py --status || true

    say "preflight passed"
    ;;

# ===========================================================================
analysis)
# ===========================================================================
    say "PHASE 1 — regenerate the analysis layer (CPU, minutes, no SLURM)"
    note "Fixes the stale BH q-values (F12) and adds the missing FDR columns (F14)."
    note "Run order is derived from the scripts' read/write dependencies; later"
    note "scripts read earlier outputs, so order matters."

    # shellcheck disable=SC1091
    source "${TFCONF_DIR}/lib/common.sh" >/dev/null 2>&1 || true
    if command -v conda >/dev/null 2>&1; then
        conda activate "${DEEPPBS_ENV:-deeppbs}" || note "WARNING: could not activate ${DEEPPBS_ENV:-deeppbs}"
    fi

    M=analysis/analyses/mechanism
    note "backing up the current data dir"
    run "cp -r ${M}/data ${M}/data.pre_fix_$(date +%Y%m%d)"

    # Dependency order (verified from read_csv/to_csv pairs):
    #   mechanism_analysis -> mechanism_table, axis_correlations,
    #                         family_effects_seedlevel, regime_contrast
    #   mechanism_probe    -> needs mechanism_table + family_effects_seedlevel
    #                      -> selffam_effects, family_shape_effects, ...
    #   mechanism_confound -> needs mechanism_table + selffam_effects
    #                      -> same_vs_other_contrast  (the published q_bh column)
    #   pilot_audit        -> needs mechanism_table + selffam_effects
    #   make_samefamily    -> samefamily_both_arms (independent)
    #   make_dnaflex*      -> need mechanism_table + samefamily_both_arms
    #   make_mechanism_fig -> figures
    for s in mechanism_analysis mechanism_probe mechanism_confound pilot_audit \
             make_samefamily make_dnaflex make_dnaflex_labeled make_mechanism_fig; do
        note "--- ${s}.py"
        run "python ${M}/${s}.py"
    done

    note "verify the BH fix actually landed in the regenerated file"
    run "python3 - <<'EOF'
import pandas as pd, numpy as np, sys
sys.path.insert(0, 'analysis/common')
from stats import bh_fdr
d = pd.read_csv('${M}/data/same_vs_other_contrast.csv').sort_values('p')
ok = bool(np.all(np.diff(d.q_bh.to_numpy()) >= -1e-12))
print('  q_bh monotonic in p:', ok)
print(d[['tf','p','q_bh']].round(4).to_string(index=False))
assert ok, 'BH still non-monotonic - the fix did not take'
EOF"

    say "analysis layer regenerated"
    note "Re-read the correction boxes in ${M}/FINDINGS.md against the new numbers."
    ;;

# ===========================================================================
tbp)
# ===========================================================================
    say "PHASE 2 — rebuild tbp from stage 4 (label fix, F1/F2)"
    note "The 190 tbp npz carry NDT80's Y_pwm INSIDE them, so relabelling cannot"
    note "fix them - they must be rebuilt. Stages 1-3 are label-independent and"
    note "are reused; stage 4 onward is regenerated."
    note ""
    note "New npz will be named 1tgh_state_NNN_TBP_HUMAN.H11MO.0.A.npz, so old and"
    note "new cannot be confused."

    note "move the contaminated npz aside rather than deleting them"
    for arm_dir in output/stage4_npz/tbp output/stage4_npz_dnarelax/tbp; do
        [ -d "${arm_dir}/output" ] || continue
        run "mkdir -p ${arm_dir}/output.ndt80_mislabelled"
        run "mv ${arm_dir}/output/*MA0343* ${arm_dir}/output.ndt80_mislabelled/ 2>/dev/null || true"
    done

    note "stale stage 5-7 artifacts for tbp"
    run "rm -rf output/stage5_aug/combined_assembly_tbp output/stage5_aug_dnarelax/combined_assembly_tbp"
    run "rm -f  output/stage7_eval/id_benchmark_tbp.json output/stage7_eval/id_benchmark_tbp_dnarelax.json"
    note "KEEPING output/stage7_eval/tbp_label_audit/ - it is the evidence of the error"

    note "resubmit stages 4-7, both arms"
    run "./scripts/pipeline/run_pilot.sh tbp 4 7"
    run "ARM=dnarelax ./scripts/pipeline/run_pilot.sh tbp 4 7"

    say "tbp rebuild submitted"
    note "Afterwards, confirm the label took:"
    note "  ls output/stage4_npz/tbp/output | head -3     # expect _TBP_HUMAN...npz"
    note "  ls output/stage4_npz/tbp/output | grep -c MA0343   # expect 0"
    ;;

# ===========================================================================
ndt80)
# ===========================================================================
    say "PHASE 3 — run the ndt80 pilot end to end (never trained)"
    note "2euv IS NDT80, so structure, PWM_LABEL and the 5 test anchors agree."
    note "This is the pilot that replaces the mislabelled tbp experiment."
    note "Needs structures/source_chains/2euv_chains/ to exist."

    if [ ! -d structures/source_chains/2euv_chains ]; then
        echo "   ERROR: structures/source_chains/2euv_chains/ missing." >&2
        echo "          Prepare it first (2euv.cif, 2euv_chainA_protein.pdb, 2euv_dna.pdb)." >&2
        exit 1
    fi
    run "./scripts/pipeline/run_pilot.sh ndt80 1 7"
    run "ARM=dnarelax ./scripts/pipeline/run_pilot.sh ndt80 1 7"
    say "ndt80 submitted"
    ;;

# ===========================================================================
reeval)
# ===========================================================================
    say "PHASE 4 — re-evaluate existing checkpoints with the fixed statistics"
    note "CHEAP: no retraining. Re-scores the checkpoints already on disk so the"
    note "JSONs carry t-distribution p-values (F9) instead of normal-approximation"
    note "ones, plus q-values, Cohen's dz and n_structures."
    note ""
    note "IMPORTANT: this does NOT fix F6. Those checkpoints were trained unseeded,"
    note "so the per-seed PAIRING is still fictional - only the p-value arithmetic"
    note "is corrected. For a genuinely paired design you need 'retrain' below."

    STAMP_LOG="$(date +%Y%m%d_%H%M%S)"
    note "log prefix for this batch: ${STAMP_LOG}_eval_<tf>[_dnarelax]_<jobid>"
    note "archive the stale JSONs"
    run "mkdir -p output/stage7_eval/pre_fix_$(date +%Y%m%d)"
    run "mv output/stage7_eval/id_benchmark_*.json output/stage7_eval/pre_fix_$(date +%Y%m%d)/ 2>/dev/null || true"

    for tf in ${PILOTS_ALL}; do
        [ "${tf}" = "ndt80" ] && continue   # no checkpoints yet
        [ "${tf}" = "tbp" ]   && continue   # contaminated; phase 2 rebuilds it
        run "sbatch --export=ALL,TF_NAME=${tf},ARM=frozen \
             --output=slurm_output/${STAMP_LOG}_eval_${tf}_%j.out \
             --error=slurm_output/${STAMP_LOG}_eval_${tf}_%j.err \
             wrappers/eval_benchmark.sh"
        run "sbatch --export=ALL,TF_NAME=${tf},ARM=dnarelax \
             --output=slurm_output/${STAMP_LOG}_eval_${tf}_dnarelax_%j.out \
             --error=slurm_output/${STAMP_LOG}_eval_${tf}_dnarelax_%j.err \
             wrappers/eval_benchmark.sh"
    done
    say "re-evaluation submitted"
    ;;

# ===========================================================================
retrain)
# ===========================================================================
    say "PHASE 5 — retrain with real seeds (F6). EXPENSIVE."
    note "Every one of the 264 existing checkpoints was trained UNSEEDED:"
    note "driver.py ignored random_seed and its seeding block was gated on a flag"
    note "that was always False. So no baseline/augmented pair shares an"
    note "initialisation, and the paired design was nominal."
    note ""
    note "Stages 1-5 are unaffected by seeding, so this is stages 6-7 only."
    note "Scope: ${PILOTS_ALL// /, }"
    note "Seeds: ${SEEDS}  x 2 arms  -> $(echo ${SEEDS} | wc -w) x 2 chains per pilot"
    note ""
    if [ "$GO" != "1" ]; then
        note "This is a dry run. Re-run with --go to submit."
    fi

    note "archive the unseeded checkpoints (do NOT delete - they are the"
    note "provenance for every number published so far)"
    STAMP="unseeded_$(date +%Y%m%d)"
    run "mkdir -p output/stage6_train/${STAMP}"
    note "    Move per pilot as each retrain completes rather than all at once,"
    note "    so a failed submit does not leave you with nothing:"
    note "      mv output/stage6_train/{baseline,augmented}_<tf>_fold0*_s* \\"
    note "         output/stage6_train/${STAMP}/"
    note ""
    note "CORRECT ENTRY POINT: run_multiseed_pilot.sh, not run_pilot.sh."
    note "  run_pilot.sh's stage 5 builds only a DEFAULT-seed config, and SEED= on"
    note "  its own just renames the run (train_compare.sh uses it for the _sN"
    note "  suffix). Only run_multiseed_pilot.sh calls build_training_configs.py"
    note "  with --seed, which is what writes random_seed into the config - and"
    note "  random_seed is the key driver.py now actually reads."
    note ""
    note "It builds N paired-seed configs, submits 2N training jobs, then one eval"
    note "job gated afterok on all of them. ARM selects the arm."
    note ""

    NS=$(echo ${SEEDS} | wc -w)
    for tf in ${PILOTS_ALL}; do
        run "./scripts/pipeline/run_multiseed_pilot.sh ${tf} ${NS}"
        run "ARM=dnarelax ./scripts/pipeline/run_multiseed_pilot.sh ${tf} ${NS}"
    done

    say "retrain submitted"
    note "Verify the seeds actually took, per run:"
    note "  grep -h 'Seeded all RNGs' slurm_output/train_*.out | sort -u"
    note "  grep -l 'UNSEEDED' slurm_output/train_*.out    # must be empty"
    note "  python3 scripts/pipeline/validate_pilot.py --all   # C10 must pass"
    ;;

# ===========================================================================
queue)
# ===========================================================================
    say "Hand everything to the continuous queue instead"
    note "Lets the queue decide what to run next and keeps N chains in flight."
    run "MAX_INFLIGHT=${MAX_INFLIGHT:-6} python3 scripts/pipeline/pilot_queue.py --fill"
    note "For a standing cron entry:"
    note "  */30 * * * * cd ${TFCONF_DIR} && MAX_INFLIGHT=6 \\"
    note "    python3 scripts/pipeline/pilot_queue.py --fill >> slurm_output/queue.log 2>&1"
    ;;

*)
    echo "Unknown phase: ${PHASE}" >&2
    usage
    ;;
esac
