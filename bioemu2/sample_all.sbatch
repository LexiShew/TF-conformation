#!/bin/bash
#SBATCH --job-name=be2all
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32GB
#SBATCH --partition=qcbgpu
#SBATCH --account=rohs_102
#SBATCH --gres=gpu:l40s:1
#SBATCH --output=/project2/rohs_102/shewchuk/TF-conformation/logs/bioemu2/be2all_%A_%a.out
#
# SLURM ARRAY worker: sample ONE source_chains structure per array task.
#
# Submit with bioemu2/launch_all.sh (do not sbatch this directly - the array
# range must match the worklist length, and the launcher is what guarantees
# that). One task = one line of bioemu2/worklist/source_chains_ok.tsv, so array
# index N is structure N and the mapping is auditable after the fact.
#
# GPU choice. gpu:l40s:1 on qcbgpu, not the gpu:rtx5000:1 the numbered stage
# wrappers use and not the generic gpu:1 of sample_complex.sbatch. Rationale:
# BioEmu2's model card names an 80 GB A100 as its smallest tested device, and
# the one verified run on this cluster (job 6537427, 1k79/ets1, 100 samples)
# completed in 9:44 with 5.7 GB peak RSS on an L40S (48 GB). L40S is the
# largest card this account reaches in quantity (22 across qcbgpu), so it is
# the throughput choice; rtx5000 (16 GB) is 3x below that and untested here.
#
# Idempotent: a structure whose sample_manifest.json already reports the
# requested frame count is skipped, so a partially-finished array can be
# resubmitted wholesale without recomputing completed work.

set -eo pipefail

TFCONF_DIR="${TFCONF_DIR:-/project2/rohs_102/shewchuk/TF-conformation}"
cd "${TFCONF_DIR}"

WORKLIST="${WORKLIST:-${TFCONF_DIR}/bioemu2/worklist/source_chains_ok.tsv}"
NUM_SAMPLES="${NUM_SAMPLES:-100}"
SEED="${SEED:-1}"
OUT_ROOT="${OUT_ROOT:-${TFCONF_DIR}/output/bioemu2}"
STATUS_DIR="${OUT_ROOT}/_status"
mkdir -p "${STATUS_DIR}"

: "${SLURM_ARRAY_TASK_ID:?this script must be submitted as an array job; use bioemu2/launch_all.sh}"

# Worklist has a header, so task N is line N+1.
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${WORKLIST}")
[ -n "${LINE}" ] || { echo "no worklist line for task ${SLURM_ARRAY_TASK_ID}"; exit 0; }

PDB_ID=$(echo "${LINE}"  | cut -f1)
CIF=$(echo "${LINE}"     | cut -f2)
NTOK=$(echo "${LINE}"    | cut -f8)
OUT="${OUT_ROOT}/${PDB_ID}_complex"

echo "[be2all] task=${SLURM_ARRAY_TASK_ID} pdb=${PDB_ID} tokens_proxy=${NTOK}"
echo "[be2all] cif=${CIF}"
echo "[be2all] out=${OUT}"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

# ---- idempotence: skip a structure already sampled at this frame count ----
MAN="${OUT}/sample_manifest.json"
if [ -f "${MAN}" ]; then
    HAVE=$(python -c "import json,sys; print(json.load(open('${MAN}')).get('n_frames_written',0))" 2>/dev/null || echo 0)
    if [ "${HAVE}" = "${NUM_SAMPLES}" ]; then
        echo "[be2all] ${PDB_ID}: already has ${HAVE} frames; skipping"
        echo "skipped" > "${STATUS_DIR}/${PDB_ID}.status"
        exit 0
    fi
    echo "[be2all] ${PDB_ID}: manifest reports ${HAVE}/${NUM_SAMPLES} frames; re-running"
fi

# Shared ColabFold MSA cache. Without this every task fetches its own MSA from
# the PUBLIC api.colabfold.com, and >8 concurrent tasks earn HTTP 429s that kill
# the task ~31 s in, before the GPU does any work (20 tasks of array 6539978 died
# exactly this way on 2026-10-22). bioemu2/warm_msa_cache.py fills this dir
# serially; with it warm a sampling task never contacts ColabFold at all.
# Reduce allocator fragmentation. Does NOT rescue a system whose single pair-
# tensor allocation exceeds the card (see the OOM ceiling in README.md), but it
# buys headroom for systems sitting just under it.
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

export BIOEMU2_MSA_CACHE="${BIOEMU2_MSA_CACHE:-${TFCONF_DIR}/output/bioemu2/_msa_cache}"

# Stagger starts so a wave of array tasks does not burst the MSA endpoint at once
# on a cold cache. No-op once the cache is warm.
sleep $(( (SLURM_ARRAY_TASK_ID % 8) * 15 ))

source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh
conda activate "${BIOEMU2_ENV:-bioemu2}"

echo "running" > "${STATUS_DIR}/${PDB_ID}.status"
T0=$SECONDS
if python "${TFCONF_DIR}/bioemu2/sample_complex.py" \
        --ref "${CIF}" \
        --out "${OUT}" \
        --num-samples "${NUM_SAMPLES}" \
        --seed "${SEED}"; then
    echo "ok $((SECONDS - T0))s" > "${STATUS_DIR}/${PDB_ID}.status"
    echo "[be2all] ${PDB_ID} DONE in $((SECONDS - T0))s"
else
    rc=$?
    echo "failed rc=${rc} $((SECONDS - T0))s" > "${STATUS_DIR}/${PDB_ID}.status"
    echo "[be2all] ${PDB_ID} FAILED rc=${rc}"
    exit ${rc}
fi
