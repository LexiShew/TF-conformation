#!/bin/bash
# render_ensemble.sh — overlay the crystal bound structure against a docked +
# minimized BioEmu ensemble for one pilot.
#
# Samples N states that made it through Stage 3, and renders, on a common DNA
# frame: the original crystal protein pose (opaque reference), the N Stage-2
# docked conformers, and the same N Stage-3 minimized conformers. Shows how far
# the apo-sampled ensemble drifts from the crystallographic binding mode.
#
# Usage:
#   scripts/viz/render_ensemble.sh <tf_name> [n_states] [seed] [out_dir]
#     n_states  number of Stage-3 states to sample (default 10)
#     seed      RNG seed for the sample (default 42; change for a different draw)
#
# Requires cluster conda envs 'bioemu' (mdtraj+matplotlib) and 'pymol' (PyMOL 3.x
# with the edu license). Run from the repo root; ray-tracing needs a compute node
# (srun ...), not the login node.
set -euo pipefail

TF="${1:?Usage: $0 <tf_name> [n_states] [seed] [out_dir]}"
NSTATES="${2:-10}"
SEED="${3:-42}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$REPO_ROOT"

source lib/common.sh
source "config/pilots/${TF}.sh"
load_pilot_config "$TF" >/dev/null 2>&1 || true
: "${PDB_ID:?PDB_ID not set}" "${PROTEIN_CHAIN:?}" "${DNA_CHAINS:?}"

OUT_DIR="${4:-output/viz/${TF}}"
mkdir -p "$OUT_DIR"
BIOEMU_PY="/home1/shewchuk/.conda/envs/bioemu/bin/python"
PYMOL_BIN="/home1/shewchuk/.conda/envs/pymol/bin/pymol"
export PYMOL_LICENSE_FILE="${PYMOL_LICENSE_FILE:-/project2/rohs_102/shewchuk/pymol/pymol-edu-license.lic}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1

REF_CIF="structures/source_chains/${PDB_ID}_chains/${PDB_ID}.cif"
S2_DIR="output/stage2_docked/${TF}"
S3_DIR="output/stage3_min/${TF}"
XTAL_PDB="${OUT_DIR}/${PDB_ID}_crystal_ref.pdb"
MANIFEST="${OUT_DIR}/.${TF}_ensemble_manifest.txt"

# 1) extract crystal protein (config chain) + bound DNA into one PDB (crystal frame)
"$BIOEMU_PY" - "$REF_CIF" "$PROTEIN_CHAIN" "$DNA_CHAINS" "$XTAL_PDB" <<'PYEOF'
import sys, mdtraj as md
ref, prot, dna, out = sys.argv[1:5]
t = md.load(ref)
sel = f"chainid {prot} or " + " or ".join(f"chainid {c}" for c in dna.split(","))
t.atom_slice(t.topology.select(sel)).save_pdb(out)
print(f"crystal ref: chain {prot} + DNA {dna} -> {out}", file=sys.stderr)
PYEOF

# 2) sample N Stage-3 states (seeded), pair each with its Stage-2 partner, build manifest
"$BIOEMU_PY" - "$S2_DIR" "$S3_DIR" "$PDB_ID" "$NSTATES" "$SEED" "$XTAL_PDB" "$MANIFEST" <<'PYEOF'
import sys, os, glob, random
s2d, s3d, pdb, n, seed, xtal, man = sys.argv[1:8]
s3 = sorted(glob.glob(f"{s3d}/{pdb}_state_*.pdb"))
random.seed(int(seed))
pick = sorted(random.sample(s3, min(int(n), len(s3))))
lines = [f"XTAL {xtal}"]
kept = 0
for p in pick:
    base = os.path.basename(p)
    s2 = os.path.join(s2d, base)
    if os.path.exists(s2):
        lines.append(f"S2 {s2}"); lines.append(f"S3 {p}"); kept += 1
open(man, "w").write("\n".join(lines) + "\n")
print(f"sampled {kept}/{len(s3)} Stage-3 states (seed {seed})", file=sys.stderr)
PYEOF

# 3) render (structure only)
RAW_PNG="${OUT_DIR}/.${TF}_ensemble_raw.png"
OUT_PNG="${OUT_DIR}/${TF}_${PDB_ID}_ensemble_${NSTATES}states.png"
RENDER_LINE=$("$PYMOL_BIN" -cq "$SCRIPT_DIR/_render_ensemble.py" -- "$MANIFEST" "$RAW_PNG" 2>/dev/null | grep RENDER_OK)
echo "$RENDER_LINE"
NS2=$(echo "$RENDER_LINE"  | sed -n 's/.*n_s2=\([0-9]*\).*/\1/p')
DEVS3=$(echo "$RENDER_LINE" | sed -n 's/.*mean_dev_s3=\([0-9.]*\).*/\1/p')
DMIN=$(echo "$RENDER_LINE"  | sed -n 's/.*dev_s3_min=\([0-9.]*\).*/\1/p')
DMAX=$(echo "$RENDER_LINE"  | sed -n 's/.*dev_s3_max=\([0-9.]*\).*/\1/p')

# 4) composite title + legend
"$BIOEMU_PY" - "$RAW_PNG" "$OUT_PNG" "$TF $PDB_ID" "$NS2" "$DEVS3" "$DMIN" "$DMAX" <<'PYEOF'
import sys, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.image as mpimg
from matplotlib.patches import Patch
raw, out, label, ns2, devs3, dmin, dmax = sys.argv[1:8]
img = mpimg.imread(raw); h, w = img.shape[:2]
fig, ax = plt.subplots(figsize=(w/200, h/200 + 1.2), dpi=200)
ax.imshow(img); ax.axis("off")
ax.set_title(f"{label}: {ns2} BioEmu conformers vs crystal bound pose\n"
             f"Stage 3 (minimized) C\u03b1 deviation from crystal: "
             f"mean {devs3} \u00c5  (range {dmin}\u2013{dmax} \u00c5)",
             fontsize=12, pad=10)
handles = [Patch(color="#b22222", label="crystal bound (PDB)"),
           Patch(color="#4169e1", label="Stage 2 docked (pre-min)"),
           Patch(color="#228b22", label="Stage 3 minimized"),
           Patch(color="#ffa500", label="crystal DNA")]
ax.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
          bbox_to_anchor=(0.5, -0.06), fontsize=10)
fig.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
print("COMPOSITE_OK", out)
PYEOF

rm -f "$RAW_PNG"
echo "WROTE $OUT_PNG"
