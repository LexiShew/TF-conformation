import os, sys, glob, re
from pymol import cmd

BASE = "/project2/rohs_102/shewchuk/TF-conformation"
STAGE2 = os.path.join(BASE, "output", "stage2_docked")
OUT = "/project2/rohs_102/shewchuk/TF-conformation/analysis/figures/pymol/_viewtest"
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
    iface = ("%s and (sidechain or name CA) and byres ( %s within 5 of %s )") % (prot, prot, dna)
    cmd.show("sticks", iface)
    cmd.show("surface", dna); cmd.set("transparency", dna_t, dna)
    cmd.set("cartoon_ring_finder", 1); cmd.set("cartoon_ring_mode", 3)
    return iface, dna

def cam_z(sel, view):
    """camera-space z of the selection centroid (PyMOL: +Z toward viewer)."""
    R = view[0:9]; org = view[12:15]
    try:
        c = cmd.centerofmass(sel)
    except Exception:
        return None
    return R[6]*(c[0]-org[0]) + R[7]*(c[1]-org[1]) + R[8]*(c[2]-org[2])

f = sorted(glob.glob(os.path.join(STAGE2, "ets1", "*_state_*.pdb")))
# state 42
target = [p for p in f if re.search(r"state_0*42\.pdb", p)][0]

def render(tag):
    cmd.png(os.path.join(OUT, "ets1_%s.png" % tag), width=900, height=900, dpi=150, ray=1)
    print("[render]", tag); sys.stdout.flush()

# ---- variant A: orient whole complex (current behavior) ----
cmd.reinitialize(); base_quality(); cmd.load(target, "mol")
iface, dna = style("mol")
cmd.orient("mol"); cmd.zoom("mol", buffer=6)
v = cmd.get_view()
zi, zd = cam_z(iface, v), cam_z(dna, v)
print("A orient-whole: z_iface=%.2f z_dna=%.2f (dna_infront=%s)" % (zi, zd, zd > zi))
render("A_orient_whole")

# ---- variant B: same, but flip 180 about y if DNA in front ----
cmd.reinitialize(); base_quality(); cmd.load(target, "mol")
iface, dna = style("mol")
cmd.orient("mol")
v = cmd.get_view()
if cam_z(dna, v) > cam_z(iface, v):
    cmd.turn("y", 180)
cmd.zoom("mol", buffer=6)
v2 = cmd.get_view()
print("B flip-if-dna-front: z_iface=%.2f z_dna=%.2f" % (cam_z(iface, v2), cam_z(dna, v2)))
render("B_flip")

# ---- variant C: orient on interface, then flip if DNA in front ----
cmd.reinitialize(); base_quality(); cmd.load(target, "mol")
iface, dna = style("mol")
cmd.orient(iface)
v = cmd.get_view()
if cam_z(dna, v) > cam_z(iface, v):
    cmd.turn("y", 180)
cmd.zoom("mol", buffer=6)
render("C_iface_flip")
