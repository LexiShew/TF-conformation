#!/usr/bin/env python3
"""
render_prepost_min.py — before/after-minimization renders for all pilots, using
the representation recipe in generic_tf_dna_settings.pml
(sidechains within 5 A as sticks, backbone as cartoon, DNA as semi-transparent
white surface).

Per pilot, for 2 sample states (best-fnat and median-fnat, from
perstate_metrics.csv):
  * <tf>_state<NNN>_before.png          — Stage-2 docked pose alone (TEAL)
  * <tf>_state<NNN>_after.png           — Stage-3 minimized pose alone (GREEN)
  * <tf>_state<NNN>_overlay.png         — docked+minimized overlaid
  * <tf>_state<NNN>_overlay_crystal.png — docked+minimized+crystal (dark grey)
  * <tf>_crystal.png                    — crystal reference alone (dark grey), once/pilot
=> 9 PNGs per pilot (4 single + 2 overlay + 2 overlay_crystal + 1 crystal).

Camera: interface_view() gives a SIDE-ON view — perpendicular to both the
protein->DNA contact axis and the DNA helical axis — so the interface sticks
face the viewer edge-on instead of being sandwiched behind the DNA surface. DNA
surface transparency is 0.8 so groove sticks read through. All poses in an
overlay are aligned on protein Ca (crystal Ca when present) so they share one
frame; interface sticks inherit each pose's color (teal=docked, green=min,
grey=crystal).

Non-destructive: SKIP_EXISTING (default on) never overwrites an existing PNG;
pass --no-skip to force re-render.

Runs headless under the `pymol` env on a COMPUTE NODE (login node has no GL).
Usage:  pymol -cq render_prepost_min.py -- <pilot|all> [more pilots ...] [--no-skip]
Reuses make_S.py's pdb_of / read_fnat / state_file helpers (same directory).
"""
import os
import sys
import csv
import re
import glob
from pymol import cmd

BASE = "/project2/rohs_102/shewchuk/TF-conformation"
OUTDIR = os.path.join(BASE, "analysis", "figures", "pymol", "prepost_min")
PERSTATE = os.path.join(BASE, "analysis", "data", "perstate_metrics.csv")
PILOTS_DIR = os.path.join(BASE, "config", "pilots")
STAGE2 = os.path.join(BASE, "output", "stage2_docked")
STAGE3 = os.path.join(BASE, "output", "stage3_min")
SRC_CHAINS = os.path.join(BASE, "structures", "source_chains")

C_S2   = "0x2E7FA8"   # TEAL       (stage2 docked = before)
C_S3   = "0x4CA83C"   # GREEN      (stage3 min    = after)
C_CRYS = "0x3A3A3A"   # dark grey  (crystal reference)

# set by main(): if True, never overwrite an existing PNG
SKIP_EXISTING = True

ALL_PILOTS = ["csl", "dux4", "egr1", "engrailed", "err", "ets1", "foxa",
              "hsf", "irf", "lef1", "nfat", "runx", "tbp"]


# ---- helpers (mirrors make_S.py) -------------------------------------------
def pdb_of(tf):
    path = os.path.join(PILOTS_DIR, "%s.sh" % tf)
    if not os.path.exists(path):
        return None
    m = re.search(r'PDB_ID=["\']?([^"\'\n#]+)', open(path).read())
    return m.group(1).strip().lower() if m else None


def read_fnat(tf):
    out = {}
    with open(PERSTATE) as fh:
        for row in csv.DictReader(fh):
            if row.get("pilot") != tf:
                continue
            if "3" not in str(row.get("stage", "")):
                continue
            try:
                fn = float(row["fnat"])
            except (ValueError, KeyError, TypeError):
                continue
            st = re.search(r"(\d+)", str(row.get("state", "")))
            if st:
                out[int(st.group(1))] = fn
    return out


def state_file(dir_tf, state_int):
    for p in sorted(glob.glob(os.path.join(dir_tf, "*.pdb"))):
        b = os.path.basename(p)
        m = re.search(r"state[_-]?0*(\d+)", b) or re.search(r"0*(\d+)\.pdb$", b)
        if m and int(m.group(1)) == state_int:
            return p
    return None


def crystal_paths(pdb):
    """Return (protein_pdb, dna_pdb) for the crystal reference, or (None, None)."""
    prot = sorted(glob.glob("%s/%s_chains/%s_chain*_protein.pdb" % (SRC_CHAINS, pdb, pdb)))
    dna = "%s/%s_chains/%s_dna.pdb" % (SRC_CHAINS, pdb, pdb)
    return (prot[0] if prot else None), (dna if os.path.isfile(dna) else None)


def _skip(out):
    """True if the file exists and SKIP_EXISTING is on (never overwrite)."""
    if SKIP_EXISTING and os.path.exists(out):
        print("[skip] exists: %s" % out)
        sys.stdout.flush()
        return True
    return False


def sample_states(tf):
    """Return up to 2 sample state indices: best-fnat and median-fnat."""
    fnat = read_fnat(tf)
    if not fnat:
        # fall back to whatever docked states exist
        cand = sorted(glob.glob(os.path.join(STAGE2, tf, "*_state_*.pdb")))
        idx = []
        for p in cand:
            m = re.search(r"state[_-]?0*(\d+)", os.path.basename(p))
            if m:
                idx.append(int(m.group(1)))
        idx.sort()
        return idx[:1] + idx[len(idx) // 2:len(idx) // 2 + 1] if idx else []
    ordered = sorted(fnat, key=fnat.get)
    best = ordered[-1]
    med = ordered[len(ordered) // 2]
    states = [best] if best == med else [best, med]
    return states


# ---- rendering -------------------------------------------------------------
import math

TURN_3Q = 30.0   # 3/4-view rotation about the camera y-axis after facing


def base_quality():
    cmd.set("cartoon_sampling", 10)
    cmd.set("cartoon_fancy_helices", 0)
    cmd.set("cartoon_flat_sheets", 1)
    cmd.set("transparency_mode", 1)
    cmd.set("ray_opaque_background", 0)
    cmd.set("antialias", 2)
    cmd.set("ray_shadows", 0)
    cmd.set("stick_quality", 15)
    cmd.bg_color("white")


def _norm(v):
    n = math.sqrt(sum(x * x for x in v)) or 1.0
    return [x / n for x in v]


def _cross(a, b):
    return [a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]]


def _dna_long_axis(dna_sel):
    """Approx DNA helical axis = unit vector between the two most distant P atoms."""
    m = cmd.get_model("(%s) and name P" % dna_sel)
    pts = [(a.coord[0], a.coord[1], a.coord[2]) for a in m.atom]
    if len(pts) < 2:
        return [0.0, 0.0, 1.0]
    best = (0.0, pts[0], pts[0])
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            d = sum((pts[i][k] - pts[j][k]) ** 2 for k in range(3))
            if d > best[0]:
                best = (d, pts[i], pts[j])
    return _norm([best[2][k] - best[1][k] for k in range(3)])


def interface_view(prot_sel, dna_sel, zoom_sel=None):
    """Side-on view of the protein-DNA interface so the interface sticks are
    NOT sandwiched behind the DNA surface.

    The camera looks perpendicular to BOTH the protein->DNA contact axis and the
    DNA helical axis: the protein sits on one side, the DNA on the other, and the
    contact patch faces the viewer edge-on (so sticks reaching into the groove
    read clearly through the semi-transparent surface). Anchored on the CONTACTED
    DNA (within 8 A of protein) so a long crystal duplex doesn't drag the view."""
    contact = "(%s) within 8 of (%s)" % (dna_sel, prot_sel)
    if cmd.count_atoms(contact) < 3:
        contact = dna_sel
    P = cmd.centerofmass(prot_sel)
    D = cmd.centerofmass(contact)
    xc = _norm([P[i] - D[i] for i in range(3)])       # interface axis -> horizontal
    L = _dna_long_axis(contact)                        # DNA helix -> vertical-ish
    zc = _norm(_cross(xc, L))                           # camera looks perp to both
    if sum(z * z for z in zc) < 1e-6:                   # xc parallel to L: fallback
        up0 = [0.0, 1.0, 0.0] if abs(xc[1]) < 0.9 else [1.0, 0.0, 0.0]
        zc = _norm(_cross(xc, up0))
    yc = _cross(zc, xc)
    view = list(cmd.get_view())
    view[0:9] = [xc[0], xc[1], xc[2], yc[0], yc[1], yc[2], zc[0], zc[1], zc[2]]
    cmd.set_view(view)
    if zoom_sel is None:
        zoom_sel = "(%s) or (%s)" % (prot_sel, contact)
    cmd.zoom(zoom_sel, buffer=6)


def style_object(obj, dna_transparency=0.8):
    """The generic_tf_dna_settings recipe, scoped to one object."""
    prot = "(%s) and polymer.protein" % obj
    dna = "(%s) and polymer.nucleic" % obj
    cmd.hide("everything", "(%s)" % obj)
    cmd.set("cartoon_side_chain_helper", 1)
    cmd.show("cartoon", prot)
    cmd.cartoon("automatic", prot)  # standard cartoon by secondary structure
    iface = ("%s and (sidechain or name CA) "
             "and byres ( %s within 5 of %s )") % (prot, prot, dna)
    cmd.show("sticks", iface)
    cmd.show("surface", dna)
    cmd.set("transparency", float(dna_transparency), dna)
    cmd.set("cartoon_ring_finder", 1)
    cmd.set("cartoon_ring_mode", 3)


def render_single(tf, state, stage_dir, tag, color):
    """Render one pose (before or after) in isolation."""
    f = state_file(os.path.join(stage_dir, tf), state)
    if not f:
        print("[%s] SKIP %s state %d: no file in %s" % (tag, tf, state, stage_dir))
        return None
    os.makedirs(OUTDIR, exist_ok=True)
    out = os.path.join(OUTDIR, "%s_state%03d_%s.png" % (tf, state, tag))
    if _skip(out):
        return out
    cmd.reinitialize()
    base_quality()
    cmd.load(f, "mol")
    style_object("mol", dna_transparency=0.8)
    cmd.color(color, "mol and polymer.protein")
    cmd.color("white", "mol and polymer.nucleic")
    interface_view("mol and polymer.protein", "mol and polymer.nucleic")
    cmd.png(out, width=1200, height=1200, dpi=200, ray=1)
    print("[render] %s" % out)
    sys.stdout.flush()
    return out


def render_overlay(tf, state):
    """Render docked (before, TEAL) + minimized (after, GREEN) overlaid."""
    s2 = state_file(os.path.join(STAGE2, tf), state)
    s3 = state_file(os.path.join(STAGE3, tf), state)
    if not s2 or not s3:
        print("[overlay] SKIP %s state %d: s2=%s s3=%s" % (tf, state, s2, s3))
        return None
    os.makedirs(OUTDIR, exist_ok=True)
    out = os.path.join(OUTDIR, "%s_state%03d_overlay.png" % (tf, state))
    if _skip(out):
        return out
    cmd.reinitialize()
    base_quality()
    cmd.load(s2, "before")
    cmd.load(s3, "after")
    # align minimized onto docked by protein Ca so the DNA context matches
    try:
        cmd.align("after and polymer.protein and name CA",
                  "before and polymer.protein and name CA")
    except Exception:
        pass
    # style both proteins (cartoon + interface sticks); one shared DNA surface
    for o, col in (("before", C_S2), ("after", C_S3)):
        prot = "%s and polymer.protein" % o
        dna = "%s and polymer.nucleic" % o
        cmd.hide("everything", o)
        cmd.show("cartoon", prot)
        cmd.cartoon("automatic", prot)
        iface = ("%s and (sidechain or name CA) "
                 "and byres ( %s within 5 of %s )") % (prot, prot, dna)
        cmd.show("sticks", iface)
        cmd.color(col, prot)
    # DNA surface: show the 'before' DNA only (they overlap), semi-transparent white
    cmd.show("surface", "before and polymer.nucleic")
    cmd.color("white", "before and polymer.nucleic")
    cmd.set("transparency", 0.8, "before and polymer.nucleic")
    cmd.set("cartoon_side_chain_helper", 1)
    interface_view("(before or after) and polymer.protein",
                   "before and polymer.nucleic")
    cmd.png(out, width=1200, height=1200, dpi=200, ray=1)
    print("[render] %s" % out)
    sys.stdout.flush()
    return out


def render_crystal(tf):
    """Crystal reference alone: protein cartoon (dark grey) + interface sticks,
    DNA semi-transparent white surface. -> <tf>_crystal.png"""
    pdb = pdb_of(tf)
    if not pdb:
        print("[crystal] SKIP %s: no PDB id" % tf)
        return None
    prot_f, dna_f = crystal_paths(pdb)
    if not prot_f or not dna_f:
        print("[crystal] SKIP %s (%s): missing crystal files (%s, %s)"
              % (tf, pdb, prot_f, dna_f))
        return None
    os.makedirs(OUTDIR, exist_ok=True)
    out = os.path.join(OUTDIR, "%s_crystal.png" % tf)
    if _skip(out):
        return out
    cmd.reinitialize()
    base_quality()
    cmd.load(prot_f, "cryP")
    cmd.load(dna_f, "cryD")
    # style: protein cartoon+interface sticks (dark grey), DNA white surface
    cmd.set("cartoon_side_chain_helper", 1)
    cmd.show("cartoon", "cryP")
    cmd.cartoon("automatic", "cryP")
    cmd.color(C_CRYS, "cryP")
    iface = ("cryP and (sidechain or name CA) "
             "and byres ( cryP within 5 of cryD )")
    cmd.show("sticks", iface)
    cmd.show("surface", "cryD")
    cmd.color("white", "cryD")
    cmd.set("transparency", 0.8, "cryD")
    cmd.set("cartoon_ring_finder", 1)
    cmd.set("cartoon_ring_mode", 3)
    # tight zoom: protein + only the DNA contact patch (within 5 A) so the long
    # crystal duplex doesn't crowd the protein/interface sticks
    interface_view("cryP", "cryD",
                   zoom_sel="cryP or (cryD within 5 of cryP)")
    cmd.png(out, width=1200, height=1200, dpi=200, ray=1)
    print("[render] %s" % out)
    sys.stdout.flush()
    return out


def render_overlay_crystal(tf, state):
    """docked (TEAL) + minimized (GREEN) + crystal (dark grey) overlaid,
    all aligned on protein Ca; one shared DNA surface (crystal DNA).
    -> <tf>_state<NNN>_overlay_crystal.png"""
    s2 = state_file(os.path.join(STAGE2, tf), state)
    s3 = state_file(os.path.join(STAGE3, tf), state)
    pdb = pdb_of(tf)
    prot_f, dna_f = crystal_paths(pdb) if pdb else (None, None)
    if not s2 or not s3 or not prot_f or not dna_f:
        print("[overlay_crystal] SKIP %s state %d: s2=%s s3=%s cryP=%s cryD=%s"
              % (tf, state, s2, s3, prot_f, dna_f))
        return None
    os.makedirs(OUTDIR, exist_ok=True)
    out = os.path.join(OUTDIR, "%s_state%03d_overlay_crystal.png" % (tf, state))
    if _skip(out):
        return out
    cmd.reinitialize()
    base_quality()
    cmd.load(s2, "before")
    cmd.load(s3, "after")
    cmd.load(prot_f, "cryP")
    cmd.load(dna_f, "cryD")
    # align docked + minimized onto the CRYSTAL protein Ca (shared frame)
    for mob in ("before", "after"):
        try:
            cmd.align("%s and polymer.protein and name CA" % mob,
                      "cryP and name CA")
        except Exception:
            pass
    # style the three proteins (cartoon + interface sticks)
    for o, col, dsel in (("before", C_S2, "before and polymer.nucleic"),
                         ("after", C_S3, "after and polymer.nucleic"),
                         ("cryP", C_CRYS, "cryD")):
        prot = "%s and polymer.protein" % o if o != "cryP" else "cryP"
        cmd.hide("everything", o)
        cmd.show("cartoon", prot)
        cmd.cartoon("automatic", prot)
        iface = ("%s and (sidechain or name CA) "
                 "and byres ( %s within 5 of %s )") % (prot, prot, dsel)
        cmd.show("sticks", iface)
        cmd.color(col, prot)
    # one shared DNA surface: the crystal DNA (dark-grey reference frame)
    cmd.show("surface", "cryD")
    cmd.color("white", "cryD")
    cmd.set("transparency", 0.8, "cryD")
    cmd.set("cartoon_side_chain_helper", 1)
    cmd.set("cartoon_ring_finder", 1)
    cmd.set("cartoon_ring_mode", 3)
    prot_all = "(before or after or cryP) and polymer.protein"
    interface_view(prot_all, "cryD",
                   zoom_sel="(%s) or (cryD within 5 of (%s))" % (prot_all, prot_all))
    cmd.png(out, width=1200, height=1200, dpi=200, ray=1)
    print("[render] %s" % out)
    sys.stdout.flush()
    return out


def do_pilot(tf):
    states = sample_states(tf)
    if not states:
        print("[pilot] SKIP %s: no sample states" % tf)
        return
    print("[pilot] %s sample states: %s" % (tf, states))
    render_crystal(tf)                       # <tf>_crystal.png (once per pilot)
    for st in states:
        render_single(tf, st, STAGE2, "before", C_S2)
        render_single(tf, st, STAGE3, "after", C_S3)
        render_overlay(tf, st)
        render_overlay_crystal(tf, st)       # <tf>_state<NNN>_overlay_crystal.png


def main():
    global SKIP_EXISTING
    argv = [a for a in sys.argv[1:] if a != "--" and not a.endswith(".py")]
    if "--no-skip" in argv:
        SKIP_EXISTING = False
        argv = [a for a in argv if a != "--no-skip"]
    pilots = argv if argv else ["all"]
    if pilots == ["all"] or "all" in pilots:
        pilots = ALL_PILOTS
    print("[render_prepost_min] pilots=%s skip_existing=%s" % (pilots, SKIP_EXISTING))
    for tf in pilots:
        do_pilot(tf)


main()
