#!/bin/bash
#SBATCH --job-name=aligncmp
#SBATCH --partition=rohs
#SBATCH --account=rohs_102
#SBATCH --gres=gpu:rtx5000:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=/project2/rohs_102/shewchuk/TF-conformation/analysis/analyses/align_compare/aligncmp_%j.out
set -eo pipefail

cd /project2/rohs_102/shewchuk/TF-conformation
source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh

BASE=/project2/rohs_102/shewchuk/TF-conformation
WORK=$BASE/analysis/analyses/align_compare/_work
NSTATES=25          # subsample per mode to bound minimization time (both arms fresh, identical settings)
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
echo "GPU: ${CUDA_VISIBLE_DEVICES:-none}"; nvidia-smi -L 2>/dev/null | head -1

# pilot: name pdb binding_chain protein_chain dna_chains
PILOTS=(
  "ets1 1k79 A 4 2,3"
  "tbp  1tgh A 2 0,1"
  "lef1 2lef A 2 0,1"
)

dock_and_min () {
  local name=$1 pdb=$2 bc=$3 pc=$4 dc=$5
  local S1=$BASE/output/stage1_bioemu/${pdb}_chain${bc}_conformations
  local REF=$BASE/structures/source_chains/${pdb}_chains/${pdb}.cif
  local OUT=$WORK/$name
  rm -rf "$OUT"; mkdir -p "$OUT"

  conda activate bioemu
  for MODE in interface all; do
    mkdir -p "$OUT/${MODE}_docked" "$OUT/${MODE}_min"
    python stage2_redock/stage2_redock.py --pdb-id $pdb --ref "$REF" \
      --traj "$S1/samples_sidechain_rec.xtc" --top "$S1/samples_sidechain_rec.pdb" \
      --out-dir "$OUT/${MODE}_docked" --protein-chain $pc --dna-chains "$dc" \
      --align-mode $MODE >/dev/null 2>&1
    # minimize a subsample of NSTATES states
    i=0; t0=$(date +%s)
    for f in $(ls "$OUT/${MODE}_docked"/${pdb}_state_*.pdb | head -$NSTATES); do
      i=$((i+1))
      s=$(basename "$f")
      python stage3_minimize/stage3_minimize.py \
        --input-pdb "$f" --output-pdb "$OUT/${MODE}_min/$s" \
        --scratch-dir /scratch1/shewchuk/deeppbs_min_tmp > "$OUT/${MODE}_min/${s}.log" 2>&1 || echo "  min fail $name $MODE $s"
      if [ $i -eq 1 ]; then echo "  [$name/$MODE] first-state platform: $(grep -h Platform "$OUT/${MODE}_min/${s}.log" 2>/dev/null | head -1)"; fi
    done
    t1=$(date +%s)
    echo "[$name/$MODE] docked=$(ls $OUT/${MODE}_docked/*.pdb|wc -l) minimized=$(ls $OUT/${MODE}_min/*.pdb 2>/dev/null|wc -l) in $((t1-t0))s"
  done
}

for p in "${PILOTS[@]}"; do dock_and_min $p; done

# --- score fnat post-minimization, both modes, write one tidy CSV ---
conda activate deeppbs 2>/dev/null || conda activate bioemu
python - <<'PY'
import sys, glob, os, csv
sys.path.insert(0, "fnat_gate")
import interface_rmsd as ir
BASE="/project2/rohs_102/shewchuk/TF-conformation"
WORK=f"{BASE}/analysis/analyses/align_compare/_work"
PILOTS={"ets1":"1k79","tbp":"1tgh","lef1":"2lef"}
rows=[]
for name,pdb in PILOTS.items():
    ref=f"{BASE}/structures/source_chains/{pdb}_chains/{pdb}.cif"
    R=ir.ref_side(ir.load_models(ref)[0], 5.0,4.5,4)   # auto-pick bound monomer (robust)
    for mode in ("interface","all"):
        for f in sorted(glob.glob(f"{WORK}/{name}/{mode}_min/{pdb}_state_*.pdb")):
            try:
                m=ir.score(R, ir.load_models(f)[0], 4.5, use_model_dna=True)
                if m["fnat"]==m["fnat"]:
                    rows.append({"pilot":name,"mode":mode,"state":os.path.basename(f),
                                 "fnat":round(m["fnat"],4),"iRMSD_global":round(m["g"],4)})
            except Exception as e:
                pass
out=f"{BASE}/analysis/analyses/align_compare/fnat_by_alignmode.csv"
with open(out,"w",newline="") as fh:
    w=csv.DictWriter(fh, fieldnames=["pilot","mode","state","fnat","iRMSD_global"]); w.writeheader()
    for r in rows: w.writerow(r)
# summary
import statistics as st
print("=== fnat by alignment mode (post-minimization) ===")
for name in PILOTS:
    for mode in ("interface","all"):
        f=[r["fnat"] for r in rows if r["pilot"]==name and r["mode"]==mode]
        if f: print(f"{name:5s} {mode:9s} n={len(f)} median={st.median(f):.3f} mean={st.mean(f):.3f} pass>=0.5={sum(x>=0.5 for x in f)/len(f):.2f}")
print("wrote", out, "rows:", len(rows))
PY
echo "JOB DONE"
