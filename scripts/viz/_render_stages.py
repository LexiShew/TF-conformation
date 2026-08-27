"""
_render_stages.py — overlay one TF-DNA state across pipeline stages (PyMOL).

Invoked by render_state_stages.sh. Produces a SINGLE aligned/superimposed view
with all three stages overlaid on a common frame so the conformational change is
directly visible:

  Stage 1 apo    (grey, thin)   — BioEmu ensemble member that seeded this state
  Stage 2 docked (marine)       — that conformer Kabsch-docked onto crystal DNA
  Stage 3 min    (forest)       — after OpenMM minimization
  DNA            (orange)       — crystal duplex, shown once (identical in S2/S3)

Alignment: Stage 3 is superposed onto Stage 2 by DNA (identical crystal duplex,
so this is an exact common frame); the Stage 1 apo protein is `super`-fit onto
the Stage 2 protein (it is the same conformer up to the rigid docking transform,
so it lands on Stage 2 — confirming docking is rigid-body and isolating the
minimization shift as the S2->S3 difference). RMSDs are printed to stdout and
annotated on the image.

Usage (from the wrapper):
  pymol -cq _render_stages.py -- <s1_apo.pdb> <s2.pdb> <s3.pdb> <out.png> <label>
"""
import sys
from pymol import cmd

s1, s2, s3, out_png, title = sys.argv[1:6]

cmd.set("ray_opaque_background", 0)
cmd.set("depth_cue", 0)
cmd.set("ray_shadows", 0)
cmd.set("antialias", 2)
cmd.set("cartoon_transparency", 0.0)
cmd.bg_color("white")

cmd.load(s1, "s1_apo")
cmd.load(s2, "s2_dock")
cmd.load(s3, "s3_min")

DNA = "polymer.nucleic"
PROT = "polymer.protein"
has_dna = cmd.count_atoms(f"s2_dock and {DNA}") > 0

# ---- bring everything into the Stage-2 frame ----
if has_dna:
    rms_dna = cmd.align(f"s3_min and {DNA}", f"s2_dock and {DNA}")[0]
else:
    rms_dna = float("nan")
# apo protein onto docked protein (sequence-based superposition, tolerant of
# the apo/docked atom-set differences)
rms_apo = cmd.super("s1_apo", f"s2_dock and {PROT}")[0]
# minimization shift: S2 vs S3 protein, all heavy atoms, in the shared DNA frame
rms_min = cmd.rms_cur(f"s3_min and {PROT} and not hydro",
                      f"s2_dock and {PROT} and not hydro")

# ---- representation: single overlaid scene ----
cmd.hide("everything")

# DNA once (from Stage 2), orange sticks+cartoon rings
if has_dna:
    cmd.show("cartoon", f"s2_dock and {DNA}")
    cmd.set("cartoon_ring_mode", 3, "s2_dock")
    cmd.set("cartoon_ring_finder", 1, "s2_dock")
    cmd.color("orange", f"s2_dock and {DNA}")

# proteins overlaid
cmd.show("cartoon", f"s1_apo and {PROT}")
cmd.color("grey70", f"s1_apo and {PROT}")
cmd.set("cartoon_transparency", 0.55, "s1_apo")

cmd.show("cartoon", f"s2_dock and {PROT}")
cmd.color("marine", f"s2_dock and {PROT}")
cmd.set("cartoon_transparency", 0.35, "s2_dock")

cmd.show("cartoon", f"s3_min and {PROT}")
cmd.color("forest", f"s3_min and {PROT}")

# Title/legend are composited afterward by the wrapper (matplotlib) — PyMOL 3D
# label placement is unreliable. Emit the RMSD so the wrapper can annotate.
cmd.orient(f"s2_dock")
cmd.turn("y", 8)
cmd.zoom("all", buffer=6, complete=1)
cmd.set("ray_trace_mode", 1)   # outline for crisp cartoons

W, H = 1800, 1350
cmd.ray(W, H)
cmd.png(out_png, dpi=200, ray=0)
print("RENDER_OK %s | rms_apo_vs_dock=%.3f rms_dna_align=%.3f rms_min_S2vsS3=%.3f"
      % (out_png, rms_apo, rms_dna, rms_min))
