#!/bin/bash
# migrate_deeppbs_to_tfconformation.sh
#
# Copy the bare-minimum pipeline files from DeepPBS/ into a stage-aligned
# layout under TF-conformation/, then (optionally) delete the originals.
#
# Defaults to DRY-RUN — prints every action without touching disk.
#   --apply       actually perform the copies
#   --cleanup     also delete the originals after the copy phase
#                 (implies --apply)
#   --src DIR     override source root  (default: $LAB/DeepPBS)
#   --dst DIR     override dest root    (default: $LAB/TF-conformation)
#   -h, --help    show this header
#
# $LAB defaults to /home/lexishew/carc/lab. Override by exporting it.
#
# Re-running with --apply is idempotent for the copy phase: cp -a refreshes
# matching files; new files are added.

set -euo pipefail

LAB="${LAB:-/home/lexishew/carc/lab}"
SRC="${LAB}/DeepPBS"
DST="${LAB}/TF-conformation"
APPLY=0
CLEANUP=0

while [ $# -gt 0 ]; do
    case "$1" in
        --apply)   APPLY=1 ;;
        --cleanup) APPLY=1; CLEANUP=1 ;;
        --src)     SRC="$2"; shift ;;
        --dst)     DST="$2"; shift ;;
        -h|--help) sed -n '2,20p' "$0"; exit 0 ;;
        *) echo "Unknown arg: $1" >&2; exit 2 ;;
    esac
    shift
done

[ -d "$SRC" ] || { echo "Source not found: $SRC" >&2; exit 1; }
[ -d "$DST" ] || { echo "Dest not found:   $DST" >&2; exit 1; }

log() { printf '%s\n' "$*"; }

# src_rel -> dst_rel mapping. One pair per line. Trailing comments after `#` ok.
MAPPINGS=$(cat <<'EOF'
# --- shared (lib/, config/, top-level) ---
deeppbs                                              lib/deeppbs
setup.py                                             lib/setup.py
MANIFEST.in                                          lib/MANIFEST.in
cg_coefficients                                      lib/cg_coefficients
dependencies/bin                                     lib/dependencies/bin
x3dna-v2.3-linux-64bit                               lib/x3dna-v2.3-linux-64bit
run/jobs/lib/common.sh                               lib/common.sh
deeppbs_linux.yml                                    config/deeppbs_linux.yml
bioemu_pdb_classification.tsv                        config/bioemu_pdb_classification.tsv
run/jobs/config/pilots                               config/pilots
run/jobs/run_pilot.sh                                run_pilot.sh
run/jobs/run_legacy_ab.sh                            run_legacy_ab.sh
run/jobs/run_multiseed_pilot.sh                      run_multiseed_pilot.sh

# --- stage 1 — hpacker ---
run/jobs/lib/stage1_hpacker.sh                       stage1_hpacker/stage1_hpacker.sh

# --- stage 2 — redock ---
run/jobs/lib/stage2_redock.sh                        stage2_redock/stage2_redock.sh
run/jobs/scripts/stage2_redock.py                    stage2_redock/stage2_redock.py

# --- stage 3 — minimize ---
run/jobs/lib/stage3_array.sh                         stage3_minimize/stage3_array.sh
run/jobs/lib/stage3_recover.sh                       stage3_minimize/stage3_recover.sh
run/jobs/scripts/stage3_minimize.py                  stage3_minimize/stage3_minimize.py

# --- stage 4 — preprocess ---
run/jobs/lib/stage4_preprocess.sh                    stage4_preprocess/stage4_preprocess.sh
run/process_co_crystal.py                            stage4_preprocess/process_co_crystal.py
run/process/proc_source.sh                           stage4_preprocess/proc_source.sh
run/process/process_config.json                      stage4_preprocess/process_config.json
run/process/process_and_predict.sh                   stage4_preprocess/process_and_predict.sh
data/jaspar_h11mo_cluster_wise_dna_containing_dataset.npy   stage4_preprocess/jaspar_h11mo_cluster_wise_dna_containing_dataset.npy

# --- stage 5 — build_aug ---
run/jobs/lib/stage5_build_aug.sh                     stage5_build_aug/stage5_build_aug.sh
run/jobs/scripts/build_augmented_fold.py             stage5_build_aug/build_augmented_fold.py
run/jobs/scripts/build_combined_assembly.py          stage5_build_aug/build_combined_assembly.py
run/jobs/scripts/build_training_configs.py           stage5_build_aug/build_training_configs.py
run/jobs/scripts/build_legacy_training_configs.py    stage5_build_aug/build_legacy_training_configs.py
run/jobs/scripts/add_crystal_to_training.sh          stage5_build_aug/add_crystal_to_training.sh
run/folds                                            stage5_build_aug/folds
run/folds_aug                                        stage5_build_aug/folds_aug

# --- stage 6 — train ---
run/jobs/lib/train_compare.sh                        stage6_train/train_compare.sh
run/jobs/lib/train_legacy_aug.sh                     stage6_train/train_legacy_aug.sh
run/driver.py                                        stage6_train/driver.py
run/models/model_v2.py                               stage6_train/models/model_v2.py
run/config.json                                      stage6_train/config_template.json

# --- stage 7 — eval ---
run/jobs/lib/eval_benchmark.sh                       stage7_eval/eval_benchmark.sh
run/jobs/lib/eval_legacy_ab.sh                       stage7_eval/eval_legacy_ab.sh
run/jobs/scripts/evaluate_id_benchmark.py            stage7_eval/evaluate_id_benchmark.py
run/predict.py                                       stage7_eval/predict.py
run/make_benchmark_set.py                            stage7_eval/make_benchmark_set.py
run/make_dataset.py                                  stage7_eval/make_dataset.py
run/process/pred_configs                             stage7_eval/pred_configs
EOF
)

# Strip comments and blank lines into a clean list.
PAIRS=$(echo "$MAPPINGS" | sed -e 's/[[:space:]]*#.*//' -e '/^[[:space:]]*$/d')

do_cp() {
    local from="$1" to="$2"
    if [ ! -e "$SRC/$from" ]; then
        log "  [skip-missing] $from"
        return
    fi
    log "  cp -a $from -> $to"
    [ "$APPLY" = 1 ] || return 0
    mkdir -p "$(dirname "$DST/$to")"
    if [ -d "$SRC/$from" ]; then
        mkdir -p "$DST/$to"
        cp -a "$SRC/$from/." "$DST/$to/"
    else
        cp -a "$SRC/$from" "$DST/$to"
    fi
}

do_rm() {
    local rel="$1"
    if [ ! -e "$SRC/$rel" ]; then
        log "  [skip-missing-rm] $rel"
        return
    fi
    log "  rm -rf $SRC/$rel"
    [ "$CLEANUP" = 1 ] || return 0
    rm -rf "$SRC/$rel"
}

mode_label=$([ $APPLY = 1 ] && echo APPLY || echo DRY-RUN)
[ $CLEANUP = 1 ] && mode_label="${mode_label} + CLEANUP"
log "=== Mode: ${mode_label} ==="
log "Source:  $SRC"
log "Dest:    $DST"
log

#-------- copy phase --------
log "### COPY PHASE ###"
while IFS=' 	' read -r src dst _; do
    [ -z "$src" ] && continue
    do_cp "$src" "$dst"
done <<< "$PAIRS"

#-------- cleanup phase (only with --cleanup; only paths we copied) --------
if [ $CLEANUP = 1 ]; then
    log
    log "### CLEANUP PHASE (removing originals we just copied) ###"
    while IFS=' 	' read -r src _ _; do
        [ -z "$src" ] && continue
        do_rm "$src"
    done <<< "$PAIRS"

    log
    log "### Tidy empty parent dirs ###"
    # These will be empty after the copied items leave; remove if so.
    for d in \
        run/jobs/lib \
        run/jobs/scripts \
        run/jobs/config \
        run/jobs \
        run/process \
        run/models \
        run \
        dependencies \
        data ; do
        if [ -d "$SRC/$d" ]; then
            if [ -z "$(ls -A "$SRC/$d" 2>/dev/null || true)" ]; then
                log "  rmdir $SRC/$d  (empty after cleanup)"
                rmdir "$SRC/$d" || true
            else
                log "  keep   $SRC/$d  (non-empty; leftover files: $(ls -A "$SRC/$d" | tr '\n' ' '))"
            fi
        fi
    done
fi

log
log "Done."
if [ $APPLY = 1 ]; then
    cat <<'NEXT'

Next steps (do these manually):
  1. Reinstall the package from its new location:
       cd "$DST/lib" && pip install -e .
     (the old DeepPBS/.eggs / deeppbs.egg-info will be stale — pip handles it,
      but you can `pip uninstall deeppbs` first if you want a clean slate.)

  2. Update path env vars in $DST/lib/common.sh. The current values point at
     /project2/.../DeepPBS/run/jobs/...; rewrite to e.g.
       REPO_DIR="$TF_CONF_REPO"
       LIB_DIR="$REPO_DIR/lib"
       CONFIG_DIR="$REPO_DIR/config"
       STAGES_DIR="$REPO_DIR"
     and replace ${LIB_DIR}/stage*.sh / ${SCRIPTS_DIR}/*.py references inside
     run_pilot.sh and the stage wrappers with the new per-stage paths.

  3. Sanity-check run_pilot.sh: the run/jobs/ copy had a hard-coded
     JOBS_DIR=/project2/...; swap to derive from $(readlink -f ...).

  4. Update relative paths in stage4_preprocess/process_and_predict.sh
     ("../process_co_crystal.py", "../../dependencies/bin/").
NEXT
else
    log
    log "This was a dry run. Re-run with --apply to copy, or --cleanup to copy+delete."
fi
