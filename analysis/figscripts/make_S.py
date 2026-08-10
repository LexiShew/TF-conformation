#!/usr/bin/env python3
"""
make_S.py — PyMOL structure renders for the TF-conformation book.

Runs headless under the `pymol` env on a COMPUTE NODE (login node segfaults —
no GL). Uses cmd.png(..., ray=1). Palette matches fig_common:
  crystal protein = GREY  0x7C8698
  best-fnat frame = TEAL  0x2E7FA8
  worst-fnat frame= ALARM 0xD1495B
  stage2 docked   = TEAL  0x2E7FA8
  stage3 min      = GREEN 0x4CA83C
  DNA context     = sand  0xCBB994  (neutral, not a palette entity)

Two figure types (per pilot):
  S1_bestworst_<tf>.png       best-fnat frame (TEAL) + worst-fnat frame (ALARM)
                              overlaid on crystal (GREY) + DNA (sand).
  S2_stage_progression_<tf>.png  crystal (GREY) vs stage2-docked (TEAL) vs
                              stage3-min (GREEN) for ONE representative state
                              (the median-fnat state) — tight overlap for a
                              rigid pilot (ets1).

Best/worst/representative states are chosen from analysis/data/perstate_metrics.csv
(stage3 rows, fnat column). State→file mapping is by the zero-padded state index
embedded in the docked/minimized PDB filenames.

Stdlib only (csv/glob/re) + pymol.cmd — the pymol env has no pandas/numpy-guaranteed.

Usage (via render_S.sbatch, on a compute node):
    pymol -cq make_S.py -- <pilot|all> [more pilots ...] [--mode S1,S2]

Defaults: pilots = ets1 lef1 ; mode = S1,S2  (S1 for both; S2 for ets1).
"""
import os
import sys
import csv
import re
import glob
from pymol import cmd

BASE = "/project2/rohs_102/shewchuk/TF-conformation"
OUTDIR = os.path.join(BASE, "analysis", "figures", "pymol")
PERSTATE = os.path.join(BASE, "analysis", "data", "perstate_metrics.csv")
PILOTS_DIR = os.path.join(BASE, "config", "pilots")
SRC_CHAINS = os.path.join(BASE, "structures", "source_chains")
STAGE2 = os.path.join(BASE, "output", "stage2_docked")
STAGE3 = os.path.join(BASE, "output", "stage3_min")

C_CRYST = "0x7C8698"   # GREY
C_BEST  = "0x2E7FA8"   # TEAL
C_WORST = "0xD1495B"   # ALARM
C_S2    = "0x2E7FA8"   # TEAL   (stage2)
C_S3    = "0x4CA83C"   # GREEN  (stage3)
C_DNA   = "0xCBB994"   # sand


# ---------------------------------------------------------------------------
def pdb_of(tf):
    """Parse PDB_ID from config/pilots/<tf>.sh (stdlib regex)."""
    path = os.path.join(PILOTS_DIR, f"{tf}.sh")
    if not os.path.exists(path):
        return None
    txt = open(path).read()
    m = re.search(r'PDB_ID=["\']?([^"\'\n#]+)', txt)
    return m.group(1).strip().lower() if m else None


def read_fnat(tf):
    """Return {state_int: fnat} for stage3 rows of this pilot from perstate_metrics."""
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
    """Find the PDB in dir_tf whose embedded state index == state_int.
    Handles both zero-padded (…_state_007…) and bare integer file indices."""
    cand = sorted(glob.glob(os.path.join(dir_tf, "*.pdb")))
    for p in cand:
        b = os.path.basename(p)
        m = re.search(r"state[_-]?0*(\d+)", b) or re.search(r"0*(\d+)\.pdb$", b)
        if m and int(m.group(1)) == state_int:
            return p
    return None


def crystal_paths(pdb):
    prot = sorted(glob.glob(f"{SRC_CHAINS}/{pdb}_chains/{pdb}_chain*_protein.pdb"))
    dna = f"{SRC_CHAINS}/{pdb}_chains/{pdb}_dna.pdb"
    return (prot[0] if prot else None), (dna if os.path.isfile(dna) else None)


def base_quality():
    cmd.set("cartoon_sampling", 8)
    cmd.set("cartoon_fancy_helices", 0)
    cmd.set("cartoon_flat_sheets", 1)
    cmd.set("transparency_mode", 1)
    cmd.set("ray_opaque_background", 0)
    cmd.set("antialias", 2)
    cmd.set("ray_shadows", 0)
    cmd.bg_color("white")


def setup_scene(pdb):
    """Load crystal protein (ref) + DNA context. Returns (ref_sel, has_dna)."""
    cmd.reinitialize()
    base_quality()
    prot, dna = crystal_paths(pdb)
    if not prot:
        raise RuntimeError(f"no crystal protein for {pdb}")
    cmd.load(prot, "cryst")
    ref_sel = "cryst and name CA"
    has_dna = False
    if dna:
        cmd.load(dna, "dna_all")
        contact = sorted({at.chain for at in cmd.get_model(
            "dna_all within 6 of cryst").atom})
        if contact:
            sel = " or ".join(f"chain {ch}" for ch in contact)
            cmd.create("dna", f"dna_all and ({sel})")
        else:
            cmd.create("dna", "dna_all")
        cmd.delete("dna_all")
        has_dna = True
    return ref_sel, has_dna


def add_frame(name, path, ref_sel):
    """Load one frame, align its Cα onto the crystal, keep protein chain A."""
    cmd.load(path, "tmp")
    try:
        cmd.align("tmp and chain A and name CA", ref_sel)
    except Exception:
        cmd.align("tmp and polymer.protein and name CA", ref_sel)
    cmd.create(name, "tmp and chain A and polymer.protein")
    cmd.delete("tmp")


def finalize_and_view(objs, has_dna):
    cmd.hide("everything")
    for o in objs:
        cmd.show("cartoon", o)
    cmd.set("cartoon_transparency", 0.0, "cryst")
    cmd.set("cartoon_tube_radius", 0.3, "cryst")
    view_sel = "cryst or dna" if has_dna else "cryst"
    if has_dna:
        cmd.show("surface", "dna")
        cmd.color("white", "dna")
        cmd.color(C_DNA, "dna")
        cmd.set("transparency", 0.55, "dna")
    cmd.orient(view_sel)
    cmd.zoom(view_sel, buffer=8)


def render_png(tag, tf):
    os.makedirs(OUTDIR, exist_ok=True)
    out = os.path.join(OUTDIR, f"{tag}_{tf}.png")
    cmd.png(out, width=1200, height=1200, dpi=200, ray=1)
    sys.stdout.flush()
    print(f"[render] {out}")
    sys.stdout.flush()
    return out


# ---------------------------------------------------------------------------
def make_S1(tf):
    pdb = pdb_of(tf)
    if not pdb:
        print(f"[S1] SKIP {tf}: no PDB_ID"); return
    fnat = read_fnat(tf)
    if not fnat:
        print(f"[S1] SKIP {tf}: no stage3 fnat rows"); return
    best_state = max(fnat, key=fnat.get)
    worst_state = min(fnat, key=fnat.get)
    best_f = state_file(os.path.join(STAGE3, tf), best_state)
    worst_f = state_file(os.path.join(STAGE3, tf), worst_state)
    if not best_f or not worst_f:
        print(f"[S1] SKIP {tf}: cannot map states "
              f"best={best_state}->{best_f} worst={worst_state}->{worst_f}")
        return
    ref_sel, has_dna = setup_scene(pdb)
    add_frame("best", best_f, ref_sel)
    add_frame("worst", worst_f, ref_sel)
    cmd.color(C_CRYST, "cryst")
    cmd.color(C_BEST, "best")
    cmd.color(C_WORST, "worst")
    cmd.set("cartoon_transparency", 0.25, "best")
    cmd.set("cartoon_transparency", 0.25, "worst")
    finalize_and_view(["cryst", "best", "worst"], has_dna)
    print(f"[S1] {tf}: best state {best_state} (fnat={fnat[best_state]:.3f}) "
          f"worst state {worst_state} (fnat={fnat[worst_state]:.3f})")
    render_png("S1_bestworst", tf)


def make_S2(tf):
    pdb = pdb_of(tf)
    if not pdb:
        print(f"[S2] SKIP {tf}: no PDB_ID"); return
    fnat = read_fnat(tf)
    if not fnat:
        print(f"[S2] SKIP {tf}: no stage3 fnat rows"); return
    # representative = median-fnat state
    ordered = sorted(fnat, key=fnat.get)
    rep = ordered[len(ordered) // 2]
    s2f = state_file(os.path.join(STAGE2, tf), rep)
    s3f = state_file(os.path.join(STAGE3, tf), rep)
    if not s2f or not s3f:
        print(f"[S2] SKIP {tf}: cannot map rep state {rep} "
              f"s2={s2f} s3={s3f}")
        return
    ref_sel, has_dna = setup_scene(pdb)
    add_frame("stage2", s2f, ref_sel)
    add_frame("stage3", s3f, ref_sel)
    cmd.color(C_CRYST, "cryst")
    cmd.color(C_S2, "stage2")
    cmd.color(C_S3, "stage3")
    cmd.set("cartoon_transparency", 0.35, "stage2")
    cmd.set("cartoon_transparency", 0.35, "stage3")
    finalize_and_view(["cryst", "stage2", "stage3"], has_dna)
    print(f"[S2] {tf}: representative state {rep} (fnat={fnat[rep]:.3f})")
    render_png("S2_stage_progression", tf)


# ---------------------------------------------------------------------------
def main():
    args = [a for a in sys.argv[1:] if a != "--" and not a.endswith(".py")]
    mode = "S1,S2"
    pilots = []
    for a in args:
        if a.startswith("--mode"):
            mode = a.split("=", 1)[1] if "=" in a else mode
        else:
            pilots.append(a)
    # allow "--mode S1,S2" as two tokens
    if "--mode" in args:
        i = args.index("--mode")
        if i + 1 < len(args):
            mode = args[i + 1]
            pilots = [p for p in pilots if p != mode]
    if not pilots:
        pilots = ["ets1", "lef1"]
    modes = set(mode.split(","))
    print(f"[make_S] pilots={pilots} modes={sorted(modes)}")
    for tf in pilots:
        if "S1" in modes:
            make_S1(tf)
        if "S2" in modes:
            make_S2(tf)


main()
