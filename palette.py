"""
palette.py — the ONE canonical color scheme for the TF-conformation project.

Import this everywhere (new figures AND the existing analysis/ + figure_scripts/
scripts) so every figure in the book uses the same hue for the same entity.

------------------------------------------------------------------------------
SEMANTIC MAPPING (color IS the cross-reference — §4.1)
------------------------------------------------------------------------------
Entity                                     role / where it appears            hue
- baseline model  /  crystal reference     grey   "the fixed reference"       GREY
- BioEmu ensemble  /  augmented·FROZEN DNA  teal   "the free-state ensemble"   TEAL
- augmented·RELAXED DNA                     green  "DNA allowed to move"       GREEN
- AlphaFold3 (comparator, motivation only)  lavender "the other tool"         AF3
- negative / "augmentation hurts" / alarm   rose   annotations & signs ONLY   ALARM

Why teal threads BioEmu -> frozen-DNA augmentation: BioEmu is the source of the
frozen-DNA augmented ensemble, so they are ONE narrative thread and share one
hue. Relaxed-DNA is the intervention on that thread, so it gets its own (green).
AF3 is a different tool shown only to motivate the method (Ch. 3), so it gets a
clearly distinct hue reserved for that comparison. The rose ALARM hue is NEVER a
data-series color — it marks "hurts"/negative/excluded only (§4.5).

Supersedes the two inconsistencies in the old scripts:
  * figure_scripts/_common.py had relaxed-DNA as BOTH orange (C_RELAX='#E07B39')
    and a green ramp -> relaxed is GREEN, always.
  * plot_diversity.py used a different teal (#4c9aa8) for BioEmu than the
    frozen-DNA teal (#3B7EA1) -> one teal, TEAL below.
"""

# ---- core entity hues (base shade) ----
GREY  = "#7C8698"   # baseline / crystal reference
TEAL  = "#2E7FA8"   # BioEmu / augmented·frozen DNA  (focal)
GREEN = "#4CA83C"   # augmented·relaxed DNA
AF3   = "#9B7FB8"   # AlphaFold3 comparator (motivation only)
ALARM = "#D1495B"   # negative / "hurts" / excluded  — annotations only, never a series

# ---- 3-shade lightness ramps for nested categories (all / other / same) (§4.3) ----
# outer level = entity hue family, inner level = light->dark within it.
GREY_R  = ["#C3C8D1", "#7C8698", "#4C5464"]   # baseline: all / other / same
TEAL_R  = ["#A9CFE0", "#5BA0C2", "#2E7FA8"]   # frozen-aug
GREEN_R = ["#BDE0A8", "#7DC162", "#4CA83C"]   # relaxed-aug

# ---- convenience maps ----
ARM_COLOR = {                # model-arm palette (results / benchmark figures)
    "baseline": GREY,
    "augmented_frozen": TEAL, "frozen": TEAL, "augmented": TEAL,
    "augmented_relaxed": GREEN, "relax": GREEN, "relaxed": GREEN,
}
SOURCE_COLOR = {"bioemu": TEAL, "af3": AF3}          # ensemble-source palette (motivation)
STRUCT_COLOR = {"crystal": GREY, "frozen": TEAL, "relaxed": GREEN}  # structural renders / dna_relax

# ---------------------------------------------------------------------------
# SEPARATE categorical axes (NOT the entity palette above).
# Some figures key color to TF-identity or to minimization-variant rather than
# to model-arm. Those are different semantic axes and get their own qualitative
# palettes — forcing them onto the 5 entity hues would break §4.1 (one hue per
# entity). Kept visually distinct from the entity hues so no cross-figure
# confusion arises.
# ---------------------------------------------------------------------------
# per-TF qualitative palette (12 pilots) — used by rmsd_analysis TF-keyed figures.
TF_PALETTE = {
    "ets1":"#2E7FA8", "tbp":"#4CA83C", "egr1":"#E1812C", "engrailed":"#8E6FB5",
    "foxa":"#D1495B", "lef1":"#C6A015", "csl":"#5AA6A0", "err":"#B5651D",
    "nfat":"#5B6FB0", "runx":"#C55FA8", "hsf":"#7A8B3A", "dux4":"#4A4A4A",
}
# NOTE: dux4 (the excluded dimer) uses a dark charcoal, deliberately distinct
# from the baseline/crystal GREY (#7C8698) — an earlier version reused GREY here,
# which collided with the baseline entity color. All 3 pilots ever drawn together
# in one figure (egr1/tbp/dux4, rmsd_analysis) are mutually distinct; the full
# 12-way set is NOT guaranteed CVD-distinct if all appear in one figure at once.
# minimization-variant axis (cage story): cage ON is the default/good arm (teal
# thread), cage OFF ("legacy") a muted warm contrast; never red (reserved alarm).
VARIANT_COLOR = {"metal_cage": TEAL, "legacy": "#C98A3A"}


# ---- shared rc (mirrors figure-style ladder; sizes 8/7/6) ----
def apply_style():
    import matplotlib as mpl
    mpl.rcParams.update({
        "figure.dpi": 120, "savefig.dpi": 300, "savefig.bbox": "tight",
        "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
        "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 7,
        "axes.grid": False, "figure.facecolor": "white", "axes.facecolor": "white",
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.7, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    })
