"""
fig_common.py — shared foundation for the F / I / R / S / P / M figure suite.

Everything the plotting scripts need that must stay consistent as pilots are
added: pilot auto-discovery, the pilot→PDB / pilot→family maps, the canonical
fnat-pass-rate ordering, labels, and the color palette (imported from the
repo-root palette.py). NOTHING here hardcodes the original six pilots — the
pilot list is discovered from disk at run time, so a newly-run TF appears in
every figure automatically once its pipeline output and eval JSON exist.

Import in every make_*.py as:  from fig_common import *
"""
import os, glob, re, json
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths — override with env vars; defaults assume this file sits in the repo
# under analysis/analyses/conformation/ (…/TF-conformation/analysis/analyses/conformation/fig_common.py).
# ---------------------------------------------------------------------------
_HERE = os.path.dirname(os.path.abspath(__file__))
TFCONF   = os.environ.get("TFCONF_DIR",   os.path.abspath(os.path.join(_HERE, "..", "..")))
DATA_DIR = os.environ.get("TFCONF_DATA",  os.path.join(TFCONF, "analysis", "data"))
FIG_DIR  = os.environ.get("TFCONF_FIGS",  os.path.join(TFCONF, "analysis", "figures"))
PILOTS_DIR = os.environ.get("TFCONF_PILOTS", os.path.join(TFCONF, "config", "pilots"))
STAGE3_DIR = os.path.join(TFCONF, "output", "stage3_min")
EVAL_DIR   = os.path.join(TFCONF, "output", "stage7_eval")

# palette lives at repo root (a copy sits beside this file for local/dev use)
import sys as _sys
_sys.path.insert(0, TFCONF)
_sys.path.insert(0, _HERE)
from palette import (GREY, TEAL, GREEN, AF3, ALARM, GREY_R, TEAL_R, GREEN_R,
                     TF_PALETTE, apply_style)

# ---------------------------------------------------------------------------
# Curated biology — the ONE place family/label knowledge lives. Add a row when
# you add a pilot; everything else is discovered. `family` is the motif-level
# DBD family; `dna_deform` is the qualitative literature call used only by M1.
# ---------------------------------------------------------------------------
PILOT_META = {
    # pilot        family                         short         dna_deform
    "ets1":     ("ETS",                           "ETS",        "minimal"),
    "tbp":      ("TBP / β-saddle",                "TBP",        "extreme (~80° kink)"),
    "egr1":     ("C2H2 zinc finger",              "C2H2-ZF",    "minimal"),
    "engrailed":("Homeodomain",                   "Homeodom.",  "modest"),
    "foxa":     ("Forkhead",                       "Forkhead",   "modest"),
    "lef1":     ("HMG-box",                        "HMG-box",    "severe (~110° bend)"),
    "csl":      ("CSL / RBPJ",                     "CSL/RBPJ",   "modest"),
    "err":      ("Nuclear receptor",              "Nuc. rec.",  "modest"),
    "nfat":     ("Rel / NFAT",                     "Rel/NFAT",   "modest"),
    "runx":     ("Runt",                           "Runt",       "modest"),
    "hsf":      ("HSF",                            "HSF",        "modest"),
    "irf":      ("IRF",                            "IRF",        "minimal"),
    "dux4":     ("Homeodomain (dimer)",           "Homeodom.²", "modest"),
}
# pilots excluded from the benchmark but still valid structural pilots (dux4 = dimer)
NONMONOMER = {"dux4"}


def parse_pilot_config(tf):
    """Read PDB_ID / PWM_LABEL / BINDING_CHAIN / FOLD from a pilot .sh config."""
    path = os.path.join(PILOTS_DIR, f"{tf}.sh")
    if not os.path.exists(path):
        return {}
    txt = open(path).read()
    out = {}
    for key in ("PDB_ID", "PWM_LABEL", "BINDING_CHAIN", "FOLD"):
        m = re.search(rf'{key}=["\']?([^"\'\n#]+)', txt)
        if m:
            out[key] = m.group(1).strip()
    return out


def discover_pilots(require=("stage3",)):
    """Return the pilot list discovered from disk (never hardcoded).

    require: which evidence a pilot must have to be included —
      'stage3' : an output/stage3_min/<pilot>/ dir (structurally processed)
      'eval'   : an output/stage7_eval/id_benchmark_<pilot>.json (benchmarked)
      'config' : a config/pilots/<pilot>.sh
    A pilot is returned only if it satisfies ALL requested evidence types.
    Base pilots only (…_dnarelax / _pass / _legacy suffixes are collapsed/ignored).
    """
    cand = set()
    for cf in glob.glob(os.path.join(PILOTS_DIR, "*.sh")):
        name = os.path.basename(cf)[:-3]
        if name.endswith("_dnarelax"):
            continue
        cand.add(name)
    for d in glob.glob(os.path.join(STAGE3_DIR, "*")):
        if os.path.isdir(d):
            name = os.path.basename(d)
            if name.endswith("_pass"):
                continue
            cand.add(name)
    def has(tf, kind):
        if kind == "stage3":
            return os.path.isdir(os.path.join(STAGE3_DIR, tf))
        if kind == "eval":
            return os.path.exists(os.path.join(EVAL_DIR, f"id_benchmark_{tf}.json"))
        if kind == "config":
            return os.path.exists(os.path.join(PILOTS_DIR, f"{tf}.sh"))
        return True
    return [tf for tf in sorted(cand) if all(has(tf, k) for k in require)]


def family_of(tf):    return PILOT_META.get(tf, (tf, tf, "unknown"))[0]
def short_of(tf):     return PILOT_META.get(tf, (tf, tf, "unknown"))[1]
def dna_deform_of(tf):return PILOT_META.get(tf, (tf, tf, "unknown"))[2]
def label_of(tf):     return f"{tf.upper()} · {short_of(tf)}"


def order_by_passrate(pilots, passrate):
    """Canonical repo ordering: descending fnat-gate pass-rate.
    passrate: dict pilot -> fraction (0..1). Missing values sort last."""
    return sorted(pilots, key=lambda t: (-(passrate.get(t, -1)), t))


def pilot_color_ordered(pilots):
    """Assign each pilot a distinct hue from TF_PALETTE (perceptually-spaced
    tab20 fallback for a pilot not yet in the palette)."""
    import matplotlib.cm as cm
    out = {}; extra = 0
    for tf in pilots:
        if tf in TF_PALETTE:
            out[tf] = TF_PALETTE[tf]
        else:
            out[tf] = cm.tab20(extra % 20); extra += 1
    return out


def savefig(fig, name, subdir=None):
    """Save into FIG_DIR (or a subdir), 300 dpi, tight. Returns the path."""
    d = FIG_DIR if subdir is None else os.path.join(FIG_DIR, subdir)
    os.makedirs(d, exist_ok=True)
    path = os.path.join(d, name)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    return path
