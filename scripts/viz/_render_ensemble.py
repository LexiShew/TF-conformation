"""
_render_ensemble.py — overlay the crystal bound structure against a docked +
minimized BioEmu ensemble for one pilot (PyMOL).

Invoked by render_ensemble.sh. Reads a manifest and produces a single scene:

  crystal (firebrick, opaque)   — the original PDB bound protein pose (ground truth)
  Stage 2 (marine, transparent) — N docked BioEmu conformers (pre-minimization)
  Stage 3 (forest, transparent) — the same N conformers after OpenMM minimization
  DNA     (orange)              — crystal duplex, shown once

Everything is registered on the crystal DNA frame (Stage 2 carries the crystal
DNA rigidly; Stage 3 is aligned back onto it), so the protein spread relative to
the bound pose is read directly. The figure visualizes how far the apo-sampled
ensemble drifts from the crystallographic binding mode.

Manifest lines:  "XTAL <path>" | "S2 <path>" | "S3 <path>"
Usage: pymol -cq _render_ensemble.py -- <manifest> <out.png>
"""
import sys
from pymol import cmd

manifest, out_png = sys.argv[1:3]

cmd.set("ray_opaque_background", 0)
cmd.set("depth_cue", 0)
cmd.set("ray_shadows", 0)
cmd.set("antialias", 2)
cmd.bg_color("white")

xtal = None
s2_objs, s3_objs = [], []
with open(manifest) as fh:
    for i, line in enumerate(fh):
        parts = line.split()
        if len(parts) != 2:
            continue
        kind, path = parts
        if kind == "XTAL":
            xtal = "xtal"; cmd.load(path, xtal)
        elif kind == "S2":
            o = f"s2_{len(s2_objs):02d}"; cmd.load(path, o); s2_objs.append(o)
        elif kind == "S3":
            o = f"s3_{len(s3_objs):02d}"; cmd.load(path, o); s3_objs.append(o)

DNA = "polymer.nucleic"
PROT = "polymer.protein"
has_dna = xtal and cmd.count_atoms(f"{xtal} and {DNA}") > 0

# ---- register every state onto the crystal DNA frame ----
if has_dna:
    for o in s2_objs + s3_objs:
        if cmd.count_atoms(f"{o} and {DNA}") > 0:
            cmd.align(f"{o} and {DNA}", f"{xtal} and {DNA}")

# ---- deviation of the minimized ensemble from the bound pose ----
# CA RMSD in the DNA-registered frame (transform=0: measure, don't move).
dev_s2, dev_s3 = [], []
for o in s2_objs:
    try:
        dev_s2.append(cmd.align(f"{o} and {PROT} and name CA",
                                f"{xtal} and {PROT} and name CA",
                                transform=0, cycles=0)[0])
    except Exception:
        pass
for o in s3_objs:
    try:
        dev_s3.append(cmd.align(f"{o} and {PROT} and name CA",
                                f"{xtal} and {PROT} and name CA",
                                transform=0, cycles=0)[0])
    except Exception:
        pass
def mean(x): return sum(x)/len(x) if x else float("nan")

# ---- representation ----
cmd.hide("everything")
# DNA once, from crystal
if has_dna:
    cmd.show("cartoon", f"{xtal} and {DNA}")
    cmd.set("cartoon_ring_mode", 3, xtal)
    cmd.set("cartoon_ring_finder", 1, xtal)
    cmd.color("orange", f"{xtal} and {DNA}")
# ensembles: transparent clouds
for o in s2_objs:
    cmd.show("cartoon", f"{o} and {PROT}"); cmd.color("marine", f"{o} and {PROT}")
    cmd.set("cartoon_transparency", 0.55, o)
for o in s3_objs:
    cmd.show("cartoon", f"{o} and {PROT}"); cmd.color("forest", f"{o} and {PROT}")
    cmd.set("cartoon_transparency", 0.55, o)
# crystal bound pose: opaque reference, drawn last
cmd.show("cartoon", f"{xtal} and {PROT}")
cmd.color("firebrick", f"{xtal} and {PROT}")
cmd.set("cartoon_transparency", 0.0, xtal)

cmd.orient(xtal)
cmd.turn("y", 8)
cmd.zoom("all", buffer=5, complete=1)
cmd.set("ray_trace_mode", 1)

W, H = 1800, 1350
cmd.ray(W, H)
cmd.png(out_png, dpi=200, ray=0)
print("RENDER_OK %s | n_s2=%d n_s3=%d mean_dev_s2=%.2f mean_dev_s3=%.2f dev_s3_min=%.2f dev_s3_max=%.2f"
      % (out_png, len(s2_objs), len(s3_objs), mean(dev_s2), mean(dev_s3),
         min(dev_s3) if dev_s3 else float("nan"), max(dev_s3) if dev_s3 else float("nan")))
