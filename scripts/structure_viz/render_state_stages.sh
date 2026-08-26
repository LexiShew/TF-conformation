#!/bin/bash
# render_state_stages.sh — visualize how one TF-DNA state changes across the
# pipeline stages (Stage 1 apo -> Stage 2 docked -> Stage 3 minimized).
#
# Produces a 4-panel PNG: apo conformer, docked complex, minimized complex, and
# a Stage2-vs-Stage3 overlay annotated with all-atom RMSD.
#
# Usage:
#   scripts/structure_viz/render_state_stages.sh <tf_name> [state_number] [out_dir]
#
#   tf_name       one of the pilot configs (tbp, egr1, foxa, lef1, engrailed, ets1, dux4)
#   state_number  1-based state index (default 1). state N == Stage-1 frame N-1.
#   out_dir       output directory (default: output/viz/<tf>)
#
# Requires (cluster conda envs): 'bioemu' (mdtraj) to slice the Stage-1 frame,
# 'pymol' (PyMOL 3.x) to render. Run from the repo root on Endeavour.
set -euo pipefail

TF="${1:?Usage: $0 <tf_name> [state] [out_dir]}"
STATE="${2:-1}"
STATE3=$(printf "%03d" "$STATE")

# resolve repo root (this script lives in scripts/structure_viz/)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$REPO_ROOT"

# load pilot config for PDB_ID + the BINDING_CHAIN stage-1 dir
source lib/common.sh
source "config/pilots/${TF}.sh"
load_pilot_config "$TF" >/dev/null 2>&1 || true
: "${PDB_ID:?PDB_ID not set by config}"

OUT_DIR="${3:-output/viz/${TF}}"
mkdir -p "$OUT_DIR"

BIOEMU_PY="/home1/shewchuk/.conda/envs/bioemu/bin/python"
PYMOL_BIN="/home1/shewchuk/.conda/envs/pymol/bin/pymol"
# Point the licensed PyMOL build at the edu license so renders have no
# "evaluation only" watermark. Override PYMOL_LICENSE_FILE in the environment
# if your license lives elsewhere.
export PYMOL_LICENSE_FILE="${PYMOL_LICENSE_FILE:-/project2/rohs_102/shewchuk/pymol/pymol-edu-license.lic}"

# Stage-1 conformers live under the BINDING_CHAIN dir. Find it (handles chainA/chainC/etc).
S1_DIR=$(ls -d output/stage1_bioemu/${PDB_ID}_chain*_conformations 2>/dev/null | head -1)
[ -z "$S1_DIR" ] && { echo "ERROR: no Stage-1 dir for $PDB_ID" >&2; exit 1; }

S2_PDB="output/stage2_docked/${TF}/${PDB_ID}_state_${STATE3}.pdb"
S3_PDB="output/stage3_min/${TF}/${PDB_ID}_state_${STATE3}.pdb"
[ -f "$S2_PDB" ] || { echo "ERROR: missing $S2_PDB" >&2; exit 1; }
[ -f "$S3_PDB" ] || { echo "ERROR: missing $S3_PDB" >&2; exit 1; }

# extract the exact Stage-1 apo frame (state N -> frame N-1)
S1_PDB="${OUT_DIR}/${PDB_ID}_state_${STATE3}_stage1.pdb"
FRAME=$((STATE - 1))
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
"$BIOEMU_PY" - "$S1_DIR" "$FRAME" "$S1_PDB" <<'PYEOF'
import sys, mdtraj
s1_dir, frame, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
# samples_sidechain_rec.xtc is the all-atom (HPacker-rebuilt) trajectory; its
# matching topology is samples_sidechain_rec.pdb (NOT topology.pdb, which is the
# backbone-only BioEmu topology with a different atom count).
t = mdtraj.load(f"{s1_dir}/samples_sidechain_rec.xtc",
                top=f"{s1_dir}/samples_sidechain_rec.pdb")
if frame >= t.n_frames:
    sys.exit(f"frame {frame} out of range ({t.n_frames} frames)")
t[frame].save_pdb(out)
print(f"stage1 frame {frame}/{t.n_frames} -> {out}", file=sys.stderr)
PYEOF

RAW_PNG="${OUT_DIR}/.${TF}_${PDB_ID}_state_${STATE3}_raw.png"
OUT_PNG="${OUT_DIR}/${TF}_${PDB_ID}_state_${STATE3}_stages.png"
LABEL="${TF} ${PDB_ID} state ${STATE3}"

RENDER_LINE=$("$PYMOL_BIN" -cq "$SCRIPT_DIR/_render_stages.py" -- \
    "$S1_PDB" "$S2_PDB" "$S3_PDB" "$RAW_PNG" "$LABEL" 2>/dev/null | grep RENDER_OK)
echo "$RENDER_LINE"
RMS_MIN=$(echo "$RENDER_LINE" | sed -n 's/.*rms_min_S2vsS3=\([0-9.]*\).*/\1/p')

# composite a clean title + color legend below the structure (matplotlib in bioemu)
"$BIOEMU_PY" - "$RAW_PNG" "$OUT_PNG" "$LABEL" "${RMS_MIN:-nan}" <<'PYEOF'
import sys, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.patches import Patch
raw, out, label, rms = sys.argv[1:5]
img = mpimg.imread(raw)
h, w = img.shape[:2]
fig, ax = plt.subplots(figsize=(w/200, h/200 + 1.1), dpi=200)
ax.imshow(img); ax.axis("off")
ax.set_title(f"{label}    Stage 2\u2192Stage 3 minimization RMSD = {rms} \u00c5",
             fontsize=13, pad=8)
handles = [Patch(color="#808080", label="Stage 1 apo (BioEmu)"),
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
