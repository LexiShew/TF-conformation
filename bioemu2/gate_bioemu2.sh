#!/bin/bash
# gate_bioemu2.sh - run the fnat gate on RAW (unminimized) BioEmu2 frames.
#
# Purpose: measure the gate pass-rate BEFORE any stage-3 minimization, to decide
# whether BioEmu2 frames need minimizing at all. The BioEmu1 pass-rates quoted in
# the repo (ETS1 100, TBP 100, EGR1 91, engrailed 91, FOXA 75, LEF1 19) are all
# POST-minimization, so this is the like-for-like question: what does the gate say
# about frames straight out of the sampler?
#
# Two steps per pilot, both required:
#   1. extract_one_copy.py  -> states_1copy/
#      MANDATORY, not cosmetic. BioEmu2 imports the whole deposited assembly, so a
#      frame can hold two protein copies and two duplexes (1K79, 3HDD). score_stage3
#      measures fnat against the model's OWN DNA (use_model_dna=True, hard-wired),
#      and with several copies present it picks the wrong partner - the first ets1
#      run this way reported median fnat 0.000, which is arithmetic, not a result.
#   2. score_stage3.py --dir states_1copy
#      with the pilot's PROTEIN_CHAIN / DNA_CHAINS so the REFERENCE side is the
#      cognate copy rather than score_stage3's contact-count auto-pick.
#
# Writes per pilot, all under output/bioemu2/<pdb>_complex/ - NOTHING is written to
# output/stage3_min*/ , so the existing BioEmu1 arms are untouched:
#   <pdb>_fnat_raw.csv        per-state fnat
#   <pdb>_pass_raw.txt        states at or above FNAT_FLOOR
#   extract_one_copy_report.csv
#
# Usage:
#   sbatch --partition=qcbr --account=rohs_102 bioemu2/gate_bioemu2.sh
#   TF_LIST="ets1 tbp" bash bioemu2/gate_bioemu2.sh      # subset, interactively
#
#SBATCH --job-name=be2gate
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --output=/project2/rohs_102/shewchuk/TF-conformation/logs/bioemu2/be2gate_%j.out
set -eo pipefail

TFCONF_DIR="${TFCONF_DIR:-/project2/rohs_102/shewchuk/TF-conformation}"
cd "${TFCONF_DIR}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
FNAT_FLOOR="${FNAT_FLOOR:-0.5}"

source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh
conda activate deeppbs        # biopython + numpy, per fnat_gate/fnat_gate.sh

SUMMARY="${TFCONF_DIR}/output/bioemu2/gate_raw_summary.tsv"
printf "pilot\tpdb\tn_states\tn_pass\tpass_rate\tmedian_fnat\tfloor\n" > "${SUMMARY}"

if [ -z "${TF_LIST:-}" ]; then
    TF_LIST=$(for c in config/pilots/*.sh; do basename "$c" .sh; done | tr '\n' ' ')
fi

for tf in ${TF_LIST}; do
    cfg="config/pilots/${tf}.sh"
    pdb=$(bash -c "set -a; source '${cfg}' >/dev/null 2>&1; echo \${PDB_ID}" | tr 'A-Z' 'a-z')
    pchain=$(bash -c "set -a; source '${cfg}' >/dev/null 2>&1; echo \${PROTEIN_CHAIN}")
    dchains=$(bash -c "set -a; source '${cfg}' >/dev/null 2>&1; echo \${DNA_CHAINS}")
    ref="${TFCONF_DIR}/structures/source_chains/${pdb}_chains/${pdb}.cif"
    cx="${TFCONF_DIR}/output/bioemu2/${pdb}_complex"
    man="${cx}/sample_manifest.json"

    if ! grep -q '"n_frames_written": 100' "${man}" 2>/dev/null; then
        echo "[gate_raw/${tf}] SKIP - no complete BioEmu2 ensemble"
        continue
    fi

    echo "[gate_raw/${tf}] ${pdb}  protein_chain=${pchain} dna_chains=${dchains}"

    if [ ! -d "${cx}/states_1copy" ] || \
       [ "$(ls "${cx}/states_1copy" 2>/dev/null | wc -l)" -lt 100 ]; then
        echo "[gate_raw/${tf}] extracting one copy per frame"
        python bioemu2/extract_one_copy.py --tf-name "${tf}" \
            --states-dir "${cx}/states" \
            --out-dir "${cx}/states_1copy" \
            --report "${cx}/extract_one_copy_report.csv" || { echo "[gate_raw/${tf}] extract FAILED"; continue; }
    else
        echo "[gate_raw/${tf}] states_1copy already present"
    fi

    python fnat_gate/score_stage3.py \
        --ref "${ref}" \
        --dir "${cx}/states_1copy" \
        --pdb-id "${pdb}" \
        --floor "${FNAT_FLOOR}" \
        --protein-chain "${pchain}" \
        --dna-chains "${dchains}" \
        --out "${cx}/${pdb}_fnat_raw.csv" \
        --pass-out "${cx}/${pdb}_pass_raw.txt" || { echo "[gate_raw/${tf}] gate FAILED"; continue; }

    python - "$tf" "$pdb" "${cx}/${pdb}_fnat_raw.csv" "${FNAT_FLOOR}" "${SUMMARY}" <<'PYEOF'
import csv, statistics, sys
tf, pdb, path, floor, summary = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), sys.argv[5]
rows = list(csv.DictReader(open(path)))
key = next((k for k in rows[0] if k.strip().lower() == "fnat"), None)
vals = [float(r[key]) for r in rows if r.get(key) not in (None, "")]
npass = sum(1 for v in vals if v >= floor)
with open(summary, "a") as fh:
    fh.write(f"{tf}\t{pdb}\t{len(vals)}\t{npass}\t{npass/len(vals):.3f}\t"
             f"{statistics.median(vals):.4f}\t{floor}\n")
print(f"[gate_raw/{tf}] n={len(vals)} pass={npass} rate={npass/len(vals):.1%} "
      f"median_fnat={statistics.median(vals):.3f}")
PYEOF
done

echo; echo "=== raw (unminimized) gate summary ==="; column -t "${SUMMARY}"
