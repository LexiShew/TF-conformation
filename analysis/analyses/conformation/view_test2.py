import os, sys, glob, re, math
from pymol import cmd

BASE = "/project2/rohs_102/shewchuk/TF-conformation"
STAGE2 = os.path.join(BASE, "output", "stage2_docked")
OUT = os.path.join(BASE, "analysis", "figures", "pymol", "_viewtest")
os.makedirs(OUT, exist_ok=True)

def base_quality():
    cmd.set("cartoon_sampling", 10); cmd.set("cartoon_fancy_helices", 0)
    cmd.set("cartoon_flat_sheets", 1); cmd.set("transparency_mode", 1)
    cmd.set("ray_opaque_background", 0); cmd.set("antialias", 2)
    cmd.set("ray_shadows", 0); cmd.set("stick_quality", 15); cmd.bg_color("white")

def style(obj, dna_t=0.5):
    prot = "(%s) and polymer.protein" % obj
    dna = "(%s) and polymer.nucleic" % obj
    cmd.hide("everything", "(%s)" % obj)
    cmd.set("cartoon_side_chain_helper", 1)
    cmd.show("cartoon", prot); cmd.cartoon("automatic", prot)
    cmd.color("0x4CA83C", prot)
    iface = ("%s and (sidechain or name CA) and byres ( %s within 5 of %s )") % (prot, prot, dna)
    cmd.show("sticks", iface)
    cmd.show("surface", dna); cmd.color("white", dna); cmd.set("transparency", dna_t, dna)
    cmd.set("cartoon_ring_finder", 1); cmd.set("cartoon_ring_mode", 3)
    return prot, dna, iface

def norm(v):
    n = math.sqrt(sum(x*x for x in v)) or 1.0
    return [x/n for x in v]
def cross(a, b):
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
def sub(a, b):
    return [a[0]-b[0], a[1]-b[1], a[2]-b[2]]

def interface_view(prot, dna, iface):
    """Look along DNA->protein so the protein (and its interface sticks) is
    nearest the camera and the DNA surface sits behind it. Roll so the
    interface contact region is centred."""
    P = cmd.centerofmass(prot)
    D = cmd.centerofmass(dna)
    # camera toward-viewer axis (z) points from DNA to protein -> protein in front
    zc = norm(sub(P, D))
    up0 = [0.0, 1.0, 0.0]
    if abs(sum(a*b for a, b in zip(up0, zc))) > 0.9:
        up0 = [1.0, 0.0, 0.0]
    xc = norm(cross(up0, zc))
    yc = cross(zc, xc)
    view = list(cmd.get_view())
    view[0:9] = [xc[0], xc[1], xc[2], yc[0], yc[1], yc[2], zc[0], zc[1], zc[2]]
    cmd.set_view(view)
    cmd.zoom("(%s) or (%s)" % (prot, dna), buffer=4)

f = sorted(glob.glob(os.path.join(STAGE2, "ets1", "*_state_*.pdb")))
target = [p for p in f if re.search(r"state_0*42\.pdb", p)][0]

# D: interface-view, transparency 0.5
cmd.reinitialize(); base_quality(); cmd.load(target, "mol")
prot, dna, iface = style("mol", 0.5)
interface_view(prot, dna, iface)
cmd.png(os.path.join(OUT, "ets1_D_ifaceview_t50.png"), width=900, height=900, dpi=150, ray=1)
print("[render] D"); sys.stdout.flush()

# E: interface-view, transparency 0.65
cmd.reinitialize(); base_quality(); cmd.load(target, "mol")
prot, dna, iface = style("mol", 0.65)
interface_view(prot, dna, iface)
cmd.png(os.path.join(OUT, "ets1_E_ifaceview_t65.png"), width=900, height=900, dpi=150, ray=1)
print("[render] E"); sys.stdout.flush()

# F: interface-view but from the DNA side rotated 30deg (a 3/4 view) t0.6
cmd.reinitialize(); base_quality(); cmd.load(target, "mol")
prot, dna, iface = style("mol", 0.6)
interface_view(prot, dna, iface)
cmd.turn("y", 30)
cmd.zoom("(%s) or (%s)" % (prot, dna), buffer=4)
cmd.png(os.path.join(OUT, "ets1_F_ifaceview_3q.png"), width=900, height=900, dpi=150, ray=1)
print("[render] F"); sys.stdout.flush()
print("VT2_DONE")
