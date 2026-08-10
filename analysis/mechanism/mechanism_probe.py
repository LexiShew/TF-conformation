#!/usr/bin/env python
"""
mechanism_probe.py -- why the measured induced-fit index fails, and what replaces it.

mechanism_analysis.py showed that a STATIC induced-fit index built from crystal DNA
geometry (bend + minor-groove deviation from canonical B-DNA) does NOT predict the sign
of the augmentation effect, and that a median split on it separates the pilots in the
wrong direction. This script establishes why, and tests the alternative that survives.

Three probes
------------
P1  Static vs dynamic DNA axes.
    Static  = how deformed the BOUND duplex is (crystal bend, crystal groove width).
    Dynamic = how much DNA shape FLUCTUATES across the ensemble (MGW-FL, bend IQR).
    The hypothesis as originally written is static. The data favour dynamic.

P2  Self-family effect.
    The cross-benchmark dPearson mixes 130 entries from every family, most unrelated to
    the pilot. The mechanistically meaningful quantity is whether augmenting with TF X's
    ensemble improves prediction on TF X's OWN family. Computed seed-level from
    perseed_perentry.csv via the selffam flag.

P3  Family-level DNA shape.
    Each benchmark family is characterised by the measured DNA geometry of the pilots
    that belong to it, then correlated with that family's seed-level augmentation effect.
    This is the "differences in performance across families" question asked directly.

Statistics
----------
Seed is the unit of replication throughout. Cross-pilot correlations are Spearman at
n = 11-13; they are reported as direction and rank, not as significance claims. Where a
median split is used it is reported with the split value so the reader can see how
unbalanced it is.

Usage
-----
    source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh && conda activate deeppbs
    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
    cd /project2/rohs_102/shewchuk/TF-conformation
    python analysis/mechanism/mechanism_probe.py

Outputs (analysis/mechanism/data/)
    static_vs_dynamic.csv     both axis families vs dPearson, side by side
    selffam_effects.csv       seed-level own-family augmentation effect per pilot
    selffam_vs_shape.csv      own-family effect vs measured DNA shape
    family_shape_effects.csv  per-family effect joined to family DNA geometry
"""
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent / "data"
OUT.mkdir(parents=True, exist_ok=True)

MECH = pd.read_csv(OUT / "mechanism_table.csv")


# ---------------------------------------------------------------- P1
STATIC = [
    ("cry_bend_uu",   "Crystal DNA bend (deg)"),
    ("bend_dev",      "Crystal bend deviation from B-DNA"),
    ("cry_minor_w",   "Crystal minor-groove width (A)"),
    ("groove_dev",    "Crystal minor-groove deviation from B"),
    ("cry_shortening","Crystal duplex shortening"),
    ("induced_fit_index", "Induced-fit index (static, measured)"),
]
DYNAMIC = [
    ("mgwfl_froz",        "MGW fluctuation, frozen (A)"),
    ("mgwfl_relax",       "MGW fluctuation, relaxed (A)"),
    ("iface_mgwfl_froz",  "Interface MGW fluctuation, frozen (A)"),
    ("iface_mgwfl_relax", "Interface MGW fluctuation, relaxed (A)"),
    ("froz_bend_iqr",     "Ensemble DNA-bend IQR, frozen (deg)"),
]
PROTEIN = [
    ("d_min",      "Reachability d_min (A)"),
    ("spread",     "Free-state ensemble spread (A)"),
    ("rmsf_mean",  "Ensemble RMSF (A)"),
]


def p1_static_vs_dynamic(outcome="froz_dP_mean"):
    rows = []
    for group, axes in (("static DNA", STATIC), ("dynamic DNA", DYNAMIC), ("protein", PROTEIN)):
        for col, label in axes:
            sub = MECH[[col, outcome, "tf"]].dropna()
            if len(sub) < 4:
                continue
            rho, p = stats.spearmanr(sub[col], sub[outcome])
            # Leave-one-out stability: how much does one pilot drive the correlation?
            loo = [stats.spearmanr(sub.drop(i)[col], sub.drop(i)[outcome])[0]
                   for i in sub.index]
            rows.append(dict(group=group, axis=col, label=label, n=len(sub),
                             rho=rho, p=p,
                             loo_min=np.min(loo), loo_max=np.max(loo),
                             loo_sign_stable=bool(np.all(np.sign(loo) == np.sign(rho))),
                             most_influential=sub.loc[
                                 sub.index[int(np.argmax(np.abs(np.array(loo) - rho)))], "tf"]))
    return (pd.DataFrame(rows)
            .sort_values("rho", key=abs, ascending=False)
            .reset_index(drop=True))


# ---------------------------------------------------------------- P2
def p2_selffam():
    """
    Own-family augmentation effect, seed level.

    perseed_perentry.csv marks each benchmark entry with selffam = True when the entry's
    family matches the augmenting pilot's family. Averaging within (pilot, dna, seed,
    selffam) and differencing arms gives a seed-paired own-family delta. LEF1 has no
    other HMG-box benchmark entries and drops out.
    """
    pe = pd.read_csv(ROOT / "analysis/figure_scripts/perseed_perentry.csv")
    g = (pe.groupby(["tf", "arm", "dna", "seed", "selffam"])
           .agg(m_pearsonr=("m_pearsonr", "mean"), m_mae=("m_mae", "mean"),
                n=("entry", "size")).reset_index())
    w = g.pivot_table(index=["tf", "dna", "seed", "selffam", "n"],
                      columns="arm", values=["m_pearsonr", "m_mae"])
    w.columns = [f"{a}_{b}" for a, b in w.columns]
    w = w.reset_index()
    w["dP"] = w["m_pearsonr_augmented"] - w["m_pearsonr_baseline"]
    w["dMAE"] = w["m_mae_augmented"] - w["m_mae_baseline"]
    w = w.dropna(subset=["dP"])

    rows = []
    for (tf, dna, sf), sub in w.groupby(["tf", "dna", "selffam"]):
        d = sub["dP"].to_numpy()
        n = len(d)
        if n > 1:
            se = d.std(ddof=1) / np.sqrt(n)
            tc = stats.t.ppf(0.975, n - 1)
            lo, hi = d.mean() - tc * se, d.mean() + tc * se
            p = stats.ttest_1samp(d, 0).pvalue
        else:
            lo = hi = p = np.nan
        rows.append(dict(tf=tf, dna=dna, selffam=bool(sf), n_seeds=n,
                         n_entries=int(sub.n.median()),
                         dP_mean=d.mean(), dP_lo=lo, dP_hi=hi, dP_p=p,
                         dMAE_mean=sub.dMAE.mean(),
                         n_neg_seeds=int((d < 0).sum())))
    return pd.DataFrame(rows)


def p2_join_shape(self_df, dna_cond="frozen"):
    """Own-family effect for one DNA condition, joined to the measured DNA-shape axes."""
    s = (self_df[(self_df.dna == dna_cond) & (self_df.selffam)]
         .set_index("tf")[["dP_mean", "dP_lo", "dP_hi", "dP_p", "n_seeds",
                           "n_entries", "n_neg_seeds"]]
         .add_prefix("selffam_").reset_index())
    keep = ["tf", "d_min", "spread", "rmsf_mean", "cry_bend_uu", "cry_minor_w",
            "cry_shortening", "bend_dev", "groove_dev", "induced_fit_index",
            "conf_selection_index", "mgwfl_froz", "mgwfl_relax", "iface_mgwfl_froz",
            "iface_mgwfl_relax", "froz_bend_iqr", "froz_dP_mean"]
    j = MECH[keep].merge(s, on="tf", how="inner")

    rows = []
    for col, label in STATIC + DYNAMIC + PROTEIN:
        sub = j[[col, "selffam_dP_mean"]].dropna()
        if len(sub) < 4:
            continue
        rho, p = stats.spearmanr(sub[col], sub["selffam_dP_mean"])
        rows.append(dict(axis=col, label=label, n=len(sub), rho=rho, p=p))
    corr = pd.DataFrame(rows).sort_values("rho", key=abs, ascending=False).reset_index(drop=True)
    return j, corr


# ---------------------------------------------------------------- P3
def p3_family_shape():
    """
    Join each benchmark family's seed-level augmentation effect to the measured DNA
    geometry of the pilot(s) representing that family. Only families with a pilot in
    the panel get a shape value.
    """
    fam = pd.read_csv(OUT / "family_effects_seedlevel.csv")
    fam = fam[fam.dna == "frozen"].copy()

    # pilot -> benchmark family label, taken from the annotation actually used by the
    # figure scripts (motif-level assignment).
    pe = pd.read_csv(ROOT / "analysis/figure_scripts/perseed_perentry.csv",
                     usecols=["tf", "family", "selffam"])
    own = (pe[pe.selffam].groupby("tf")["family"]
           .agg(lambda s: s.mode().iat[0]).rename("family").reset_index())

    shape = MECH[["tf", "cry_bend_uu", "cry_minor_w", "bend_dev", "groove_dev",
                  "induced_fit_index", "d_min", "spread",
                  "mgwfl_froz", "iface_mgwfl_froz", "froz_bend_iqr"]]
    fam_shape = (own.merge(shape, on="tf", how="left")
                    .groupby("family").agg(
                        pilots=("tf", lambda s: " ".join(sorted(s))),
                        n_pilots=("tf", "nunique"),
                        cry_bend=("cry_bend_uu", "mean"),
                        cry_minor_w=("cry_minor_w", "mean"),
                        induced_fit_index=("induced_fit_index", "mean"),
                        d_min=("d_min", "mean"),
                        mgwfl_froz=("mgwfl_froz", "mean"),
                        iface_mgwfl_froz=("iface_mgwfl_froz", "mean"),
                        bend_iqr=("froz_bend_iqr", "mean"))
                    .reset_index())

    j = fam.merge(fam_shape, on="family", how="left")
    rows = []
    for col in ["cry_bend", "cry_minor_w", "induced_fit_index", "d_min",
                "mgwfl_froz", "iface_mgwfl_froz", "bend_iqr"]:
        sub = j[[col, "dP_mean"]].dropna()
        if len(sub) < 4:
            continue
        rho, p = stats.spearmanr(sub[col], sub["dP_mean"])
        rows.append(dict(axis=col, n_families=len(sub), rho=rho, p=p))
    return j, pd.DataFrame(rows).sort_values("rho", key=abs, ascending=False).reset_index(drop=True)


def main():
    p1 = p1_static_vs_dynamic()
    p1.to_csv(OUT / "static_vs_dynamic.csv", index=False)
    print("=== P1  static vs dynamic DNA axes vs cross-benchmark dPearson (frozen) ===")
    print(p1[["group", "label", "n", "rho", "p", "loo_min", "loo_max",
              "loo_sign_stable", "most_influential"]].round(3).to_string(index=False))

    sf = p2_selffam()
    sf.to_csv(OUT / "selffam_effects.csv", index=False)
    tab = sf[(sf.dna == "frozen") & sf.selffam].sort_values("dP_mean", ascending=False)
    print("\n=== P2  OWN-FAMILY augmentation effect, seed level (frozen) ===")
    print(tab[["tf", "n_seeds", "n_entries", "dP_mean", "dP_lo", "dP_hi",
               "dP_p", "n_neg_seeds"]].round(4).to_string(index=False))

    print("\n=== P2a  OWN-FAMILY effect, FROZEN vs RELAX side by side ===")
    both = (sf[sf.selffam]
            .pivot_table(index="tf", columns="dna",
                         values=["dP_mean", "dP_p", "n_seeds", "n_neg_seeds"]))
    both.columns = [f"{a}_{b}" for a, b in both.columns]
    both = both.sort_values("dP_mean_frozen", ascending=False)
    print(both[["dP_mean_frozen", "dP_p_frozen", "n_neg_seeds_frozen",
                "dP_mean_relax", "dP_p_relax", "n_neg_seeds_relax"]]
          .round(4).to_string())
    both.to_csv(OUT / "selffam_frozen_vs_relax.csv")

    j, corr = p2_join_shape(sf)
    j.to_csv(OUT / "selffam_vs_shape.csv", index=False)
    print("\n=== P2b  own-family effect vs measured shape (frozen) ===")
    print(corr.round(3).to_string(index=False))

    fj, fcorr = p3_family_shape()
    fj.to_csv(OUT / "family_shape_effects.csv", index=False)
    print("\n=== P3  family effect joined to family DNA shape ===")
    print(fj[["family", "pilots", "n_seeds", "dP_mean", "dP_p", "cry_bend",
              "cry_minor_w", "d_min", "iface_mgwfl_froz"]]
          .round(3).to_string(index=False))
    print("\n=== P3b  family shape vs family effect ===")
    print(fcorr.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
