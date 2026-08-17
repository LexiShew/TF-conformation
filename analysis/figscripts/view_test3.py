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

def style(obj, dna_t):
    prot = "(%s) and polymer.protein" % obj
    dna = "(%s) and polymer.nucleic" % obj
    cmd.hide("everything", "(%s)" % obj)
    cmd.set("cartoon_side_chain_helper", 1)
    cmd.show("cartoon", prot); cmd.cartoon("automatic", prot); cmd.color("0x2E7FA8", prot)
    iface = ("%s and (sidechain or name CA) and byres ( %s within 5 of %s )") % (prot, prot, dna)
    cmd.show("sticks", iface)
    # emphasise interface sticks with a warm colour so they pop against teal/white
    cmd.color("0xE08A2B", iface)
    cmd.show("surface", dna); cmd.color("white", dna); cmd.set("transparency", dna_t, dna)
    cmd.set("cartoon_ring_finder", 1); cmd.set("cartoon_ring_mode", 3)
    return prot, dna, iface

def _norm(v):
    n=math.sqrt(sum(x*x for x in v)) or 1.0; return [x/n for x in v]
def _cross(a,b):
    return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
def _sub(a,b): return [a[0]-b[0],a[1]-b[1],a[2]-b[2]]

def dna_long_axis(dna_sel):
    """Approx DNA helical axis = vector between the two most distant P atoms."""
    m = cmd.get_model("(%s) and name P" % dna_sel)
    pts = [(a.coord[0],a.coord[1],a.coord[2]) for a in m.atom]
    if len(pts) < 2:
        return [0.0,0.0,1.0]
    best=(0,pts[0],pts[0])
    for i in range(len(pts)):
        for j in range(i+1,len(pts)):
            d=sum((pts[i][k]-pts[j][k])**2 for k in range(3))
            if d>best[0]: best=(d,pts[i],pts[j])
    return _norm(_sub(best[2],best[1]))

def set_axes(x,y,z):
    v=list(cmd.get_view()); v[0:9]=[x[0],x[1],x[2],y[0],y[1],y[2],z[0],z[1],z[2]]; cmd.set_view(v)

f = sorted(glob.glob(os.path.join(STAGE2, "ets1", "*_state_*.pdb")))
target=[p for p in f if re.search(r"state_0*42\.pdb",p)][0]

def contact_region(prot,dna):
    return "((%s) within 6 of (%s)) or ((%s) within 6 of (%s))" % (prot,dna,dna,prot)

# ---- G: orient on the interface region (face-on to the contact patch), t=0.6 ----
cmd.reinitialize(); base_quality(); cmd.load(target,"mol")
prot,dna,iface=style("mol",0.6)
cmd.orient(contact_region(prot,dna))
cmd.zoom(contact_region(prot,dna), buffer=8)
cmd.png(os.path.join(OUT,"ets1_G_orientiface_t60.png"),width=900,height=900,dpi=150,ray=1)
print("[render] G"); sys.stdout.flush()

# ---- H: side-on (x=P->D, z perp to x and DNA-long-axis), t=0.65 ----
cmd.reinitialize(); base_quality(); cmd.load(target,"mol")
prot,dna,iface=style("mol",0.65)
P=cmd.centerofmass(prot); D=cmd.centerofmass(dna)
xc=_norm(_sub(P,D))                 # interface axis -> horizontal
L=dna_long_axis(dna)
zc=_norm(_cross(xc,L))              # look perp to interface axis AND DNA helix
yc=_cross(zc,xc)
set_axes(xc,yc,zc)
cmd.zoom(contact_region(prot,dna), buffer=8)
cmd.png(os.path.join(OUT,"ets1_H_sideon_t65.png"),width=900,height=900,dpi=150,ray=1)
print("[render] H"); sys.stdout.flush()

# ---- I: side-on, higher transparency t=0.8 ----
cmd.reinitialize(); base_quality(); cmd.load(target,"mol")
prot,dna,iface=style("mol",0.8)
P=cmd.centerofmass(prot); D=cmd.centerofmass(dna)
xc=_norm(_sub(P,D)); L=dna_long_axis(dna); zc=_norm(_cross(xc,L)); yc=_cross(zc,xc)
set_axes(xc,yc,zc)
cmd.zoom(contact_region(prot,dna), buffer=8)
cmd.png(os.path.join(OUT,"ets1_I_sideon_t80.png"),width=900,height=900,dpi=150,ray=1)
print("[render] I"); sys.stdout.flush()

# ---- J: orient on interface region, high transparency t=0.8 ----
cmd.reinitialize(); base_quality(); cmd.load(target,"mol")
prot,dna,iface=style("mol",0.8)
cmd.orient(contact_region(prot,dna))
cmd.zoom(contact_region(prot,dna), buffer=8)
cmd.png(os.path.join(OUT,"ets1_J_orientiface_t80.png"),width=900,height=900,dpi=150,ray=1)
print("[render] J"); sys.stdout.flush()
print("VT3_DONE")
