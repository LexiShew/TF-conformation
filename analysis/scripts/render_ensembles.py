#!/usr/bin/env python3
"""Render stacked crystal + AF3 + BioEmu ensembles per pilot (PyMOL, headless).

For each pilot, all protein conformations are superposed on the crystal
binding-chain Cα:
  - crystal protein : opaque dark-grey cartoon (the reference bound pose)
  - crystal DNA     : gold cartoon, shown once for binding context
  - AF3 (10 models) : lavender cartoon, semi-transparent (tight bundle)
  - BioEmu (~90)    : teal cartoon, high transparency (broad fan)

Colours match the D1 diversity figure. Writes analysis/figures/pymol/<pilot>_stack.png.

Usage:  pymol -cq render_ensembles.py -- <pilot|all>
Run in the `pymol` conda env; cap BLAS threads on the login node.
"""
import os, sys, glob
from pymol import cmd

BASE = "/project2/rohs_102/shewchuk/TF-conformation"
OUTDIR = os.path.join(BASE, "analysis", "figures", "pymol")

PILOTS = {   # pilot -> (pdb, af3 subdir, docked subdir)
    "ets1":      ("1k79", "ets1_1k79",      "ets1"),
    "tbp":       ("1tgh", "tbp_1tgh",       "tbp"),
    "egr1":      ("1aay", "egr1_1aay",      "egr1"),
    "engrailed": ("3hdd", "engrailed_3hdd", "engrailed"),
    "foxa":      ("1vtn", "foxa_1vtn",      "foxa"),
    "lef1":      ("2lef", "lef1_2lef",      "lef1"),
}
ORDER = ["ets1", "tbp", "egr1", "engrailed", "foxa", "lef1"]

C_CRYST = "0x404040"   # dark grey
C_DNA   = "0xcbb994"   # muted warm sand (soft context, not competing)
C_AF3   = "0xb9a7d6"   # lavender  (matches D1)
C_BIO   = "0x4c9aa8"   # teal      (matches D1)


def build_stack(name, paths, ref_sel, max_frames=None):
    """Load each protein-chain-A structure, align its Cα onto ref_sel, and
    accumulate the aligned protein into a single multi-state object `name`.
    If max_frames is set, evenly subsample the input paths (a conformational
    cloud reads the same with ~30 frames as with ~90, and cartoon ray-tracing
    cost scales with frame count)."""
    if max_frames and len(paths) > max_frames:
        step = len(paths) / float(max_frames)
        paths = [paths[int(i * step)] for i in range(max_frames)]
    st = 0
    for p in paths:
        cmd.load(p, "tmp")
        # align protein chain A Cα onto the crystal reference
        try:
            cmd.align("tmp and chain A and name CA", ref_sel)
        except Exception:
            cmd.align("tmp and polymer.protein and name CA", ref_sel)
        st += 1
        cmd.create(name, "tmp and chain A and polymer.protein", 1, st)
        cmd.delete("tmp")
    return st


def render_pilot(pilot):
    pdb, af3dir, dockdir = PILOTS[pilot]
    cmd.reinitialize()
    cmd.bg_color("white")

    # crystal reference protein (single binding chain)
    cryst_prot = sorted(glob.glob(f"{BASE}/structures/source_chains/{pdb}_chains/{pdb}_chain*_protein.pdb"))[0]
    cmd.load(cryst_prot, "cryst")
    ref_sel = "cryst and name CA"

    # crystal DNA (context). Some crystals (e.g. 1k79) contain a second DNA
    # duplex from a packing copy; keep only chains contacting the binding
    # protein so the view isn't dominated by a distant duplex.
    dna = f"{BASE}/structures/source_chains/{pdb}_chains/{pdb}_dna.pdb"
    has_dna = os.path.isfile(dna)
    if has_dna:
        cmd.load(dna, "dna_all")
        contact_chains = sorted({at.chain for at in cmd.get_model(
            "dna_all within 5 of cryst").atom})
        if contact_chains:
            sel = " or ".join(f"chain {ch}" for ch in contact_chains)
            cmd.create("dna", f"dna_all and ({sel})")
        else:
            cmd.create("dna", "dna_all")
        cmd.delete("dna_all")

    # AF3 bundle + BioEmu sample (each aligned onto crystal)
    af3_paths = sorted(glob.glob(f"{BASE}/af3/output/{af3dir}/seed-*/*_model.cif"))
    bio_paths = sorted(glob.glob(f"{BASE}/output/stage2_docked/{dockdir}/*.pdb"))
    n_af3 = build_stack("af3", af3_paths, ref_sel)                # all 10
    n_bio = build_stack("bio", bio_paths, ref_sel, max_frames=10) # sample 10

    # global cartoon / quality settings (once)
    cmd.set("cartoon_sampling", 6)            # coarser spline = cheaper
    cmd.set("cartoon_fancy_helices", 0)
    cmd.set("cartoon_flat_sheets", 1)
    cmd.set("all_states", 1)                  # show every state of a stack
    cmd.set("transparency_mode", 1)           # order-independent (no depth sort)
    cmd.set("ray_opaque_background", 0)
    cmd.set("antialias", 2)
    cmd.set("ray_shadows", 0)
    cmd.set("ray_trace_mode", 0)
    cmd.color(C_CRYST, "cryst")
    cmd.color(C_AF3, "af3")
    cmd.color(C_BIO, "bio")
    if has_dna:
        cmd.set("cartoon_ring_mode", 3, "dna")
        cmd.set("cartoon_ring_finder", 1, "dna")
        cmd.color(C_DNA, "dna")

    # fixed view for both panels (frame on the crystal complex)
    view_sel = "cryst or dna" if has_dna else "cryst"
    cmd.orient(view_sel)
    cmd.zoom(view_sel, buffer=8)

    os.makedirs(OUTDIR, exist_ok=True)

    def render(ensemble_obj, ens_transp, tag):
        """Show crystal (opaque dark reference) + one ensemble (transparent)
        + DNA; hide the other ensemble; ray-trace to <pilot>_<tag>.png."""
        cmd.hide("everything")
        cmd.show("cartoon", f"cryst or {ensemble_obj}")
        cmd.set("cartoon_transparency", 0.0, "cryst")
        cmd.set("cartoon_tube_radius", 0.3, "cryst")
        cmd.set("cartoon_transparency", ens_transp, ensemble_obj)
        if has_dna:
            cmd.show("cartoon", "dna")
            cmd.set("cartoon_transparency", 0.10, "dna")
        out = os.path.join(OUTDIR, f"{pilot}_{tag}.png")
        cmd.png(out, width=1200, height=1200, dpi=200, ray=1)
        sys.stdout.flush()
        print(f"{pilot}_{tag}: -> {out}")
        sys.stdout.flush()

    # Panel 1: AF3 (transparent purple).  Panel 2: 10 BioEmu (transparent teal).
    # AF3's 10 models nearly coincide, so 10 stacked cartoon layers compound
    # toward opaque; use high transparency so the crystal shows through.
    render("af3", 0.75, "af3")
    render("bio", 0.55, "bioemu")
    print(f"{pilot}: AF3 states={n_af3} BioEmu states={n_bio} DNA={has_dna}")
    sys.stdout.flush()


def main():
    # PyMOL puts args after `--` into sys.argv; be tolerant of both launch modes
    args = [a for a in sys.argv[1:] if a != "--"]
    target = args[-1] if args else "all"
    if target.endswith(".py"):        # only the script path present
        target = "all"
    pilots = ORDER if target == "all" else [target]
    for p in pilots:
        render_pilot(p)


# PyMOL executes this file with __name__ != "__main__"; run unconditionally.
main()
