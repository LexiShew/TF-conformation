#!/usr/bin/env python
"""
mechanism_analysis.py -- Induced fit vs conformational selection, tested quantitatively.

Question
--------
The project's mechanism hypothesis (analysis/mechanism_and_roadmap.md) says
conformational augmentation of DeepPBS helps in the CONFORMATIONAL-SELECTION regime
(free protein already reaches its bound pose; recognition needs little DNA deformation)
and hurts under INDUCED FIT ON THE DNA (the model holds DNA geometry fixed, so
free-protein frames carry no information about the distortion recognition requires).

That hypothesis was originally supported with a CURATED, literature-derived
DNA-deformation axis ("minimal / modest / severe"). Every DNA-shape quantity it
needs is now measured on disk. This script replaces the curated axis with measured
pyCurves geometry and tests the hypothesis as a statistical claim.

What it does
------------
1. Builds one per-pilot mechanism table from:
     protein side  : analysis/data/reachability.csv          (d_min, spread, rmsf, reach_ratio)
     DNA side      : analysis/dna_relax/data/pycurves_all_perstructure.csv  (crystal + ensemble geometry)
                     analysis/dna_relax/data/mgw_fl_summary.csv             (minor-groove-width fluctuation)
                     analysis/dna_relax/data/iface_mgwfl_vs_accuracy.csv    (interface-restricted MGW-FL)
     outcome       : analysis/figure_scripts/perseed_summary.csv            (seed-paired dPearson)
2. Defines a MEASURED induced-fit index from crystal DNA geometry (deviation of the
   bound duplex from canonical B-DNA), replacing the curated dna_deform labels.
3. Correlates every protein-side and DNA-side axis against the seed-paired
   augmentation effect, and reports which side of the interface actually predicts it.
4. Runs the family-level analysis at the correct unit of replication (seed, not entry).

Statistics
----------
The unit of replication for an augmentation effect is the SEED (augmentation is applied
once per retraining), not the benchmark entry. All per-pilot effects are seed-level
means with seed-level CIs; family effects use seed-level paired tests. This follows the
correction in analysis/figure_scripts/REVIEW_figure_scripts.md, which found that treating
entry x seed rows as independent replicates shrinks the standard error ~2x.

Correlations across pilots are Spearman with n = 11-13 pilots. At that n almost nothing
reaches significance; the script reports rho with p and CI and the text below is written
to describe direction and rank, not to claim significance.

Usage
-----
    source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh && conda activate deeppbs
    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
    cd /project2/rohs_102/shewchuk/TF-conformation
    python analysis/mechanism/mechanism_analysis.py

Outputs (analysis/mechanism/)
    data/mechanism_table.csv        one row per pilot, every axis + outcome
    data/axis_correlations.csv      each axis vs dPearson: rho, p, n
    data/family_effects_seedlevel.csv  per-family seed-level augmentation effect
    data/regime_contrast.csv        conformational-selection vs induced-fit split
"""
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
(OUT / "data").mkdir(parents=True, exist_ok=True)

# Canonical B-DNA reference values (Curves+ convention), used to express how far a
# bound duplex is deformed from ideal B-form. Minor-groove width 5.7 A is the value the
# project already uses in analysis/dna_relax/README.md; bend 0 deg is ideal straight B-DNA.
B_DNA_MINOR_W = 5.7
B_DNA_BEND = 0.0

# DNA-condition labels as they appear in the benchmark CSVs. These are asserted against
# the data at load time rather than assumed: a mismatch here drops a whole arm silently.
DNA_COND_FROZEN = "frozen"
DNA_COND_RELAX = "relax"


# ---------------------------------------------------------------- loaders
def load_protein_axes():
    """Free-state protein ensemble geometry: how far the apo ensemble is from bound."""
    df = pd.read_csv(ROOT / "analysis/data/reachability.csv")
    return df.rename(columns={"pilot": "tf", "n": "n_states_reach"})


def load_dna_geometry():
    """Measured DNA shape, per pilot, split by condition (crystal / frozen / relaxed)."""
    ps = pd.read_csv(ROOT / "analysis/dna_relax/data/pycurves_all_perstructure.csv")

    cry = (ps[ps.cond == "crystal"]
           .set_index("tf")[["bend_uu", "bend_pp", "minor_w", "major_w", "shortening"]]
           .add_prefix("cry_"))

    # Ensemble spread of DNA shape under the frozen pipeline: how much the DNA geometry
    # varies across docked/minimized frames even when it is tethered.
    froz = ps[ps.cond == "frozen"].groupby("tf").agg(
        froz_bend_med=("bend_uu", "median"),
        froz_bend_iqr=("bend_uu", lambda s: s.quantile(.75) - s.quantile(.25)),
        froz_minor_med=("minor_w", "median"),
        n_froz=("state", "size"),
    )
    rel = ps[ps.cond == "relaxed"].groupby("tf").agg(
        rel_bend_med=("bend_uu", "median"),
        rel_minor_med=("minor_w", "median"),
        n_rel=("state", "size"),
    )
    return cry.join(froz, how="outer").join(rel, how="outer").reset_index()


def load_mgwfl():
    """Minor-groove-width fluctuation: whole-molecule and interface-restricted."""
    whole = (pd.read_csv(ROOT / "analysis/dna_relax/data/mgw_fl_summary.csv")
             .rename(columns={"pilot": "tf",
                              "mgw_fl_frozen_mean": "mgwfl_froz",
                              "mgw_fl_relaxed_mean": "mgwfl_relax"})
             [["tf", "mgwfl_froz", "mgwfl_relax"]])
    iface = (pd.read_csv(ROOT / "analysis/dna_relax/data/iface_mgwfl_vs_accuracy.csv")
             .rename(columns={"pilot": "tf"})
             [["tf", "iface_mgwfl_af3", "iface_mgwfl_froz", "iface_mgwfl_relax", "n_iface"]])
    return whole.merge(iface, on="tf", how="outer")


def load_seed_effects():
    """Seed-paired augmentation effect. The seed is the unit of replication."""
    ss = pd.read_csv(ROOT / "analysis/figure_scripts/perseed_summary.csv")
    wide = ss.pivot_table(index=["tf", "dna", "seed"], columns="arm",
                          values=["mean_pearsonr", "mean_mae"])
    wide.columns = [f"{a}_{b}" for a, b in wide.columns]
    wide = wide.reset_index()
    wide["dP"] = wide["mean_pearsonr_augmented"] - wide["mean_pearsonr_baseline"]
    wide["dMAE"] = wide["mean_mae_augmented"] - wide["mean_mae_baseline"]
    wide = wide.dropna(subset=["dP"])

    found = set(wide.dna.unique())
    expected = {DNA_COND_FROZEN, DNA_COND_RELAX}
    if not expected <= found:
        raise ValueError(
            f"DNA condition labels in perseed_summary.csv are {sorted(found)}, "
            f"expected to contain {sorted(expected)}. Filtering on a wrong label "
            f"silently drops an entire arm -- fix DNA_COND_* rather than the filter.")
    return wide


def summarize_seed_effects(seed_df):
    """Per pilot x DNA-condition: seed-level mean, CI, and paired t-test (n = n_seeds)."""
    rows = []
    for (tf, dna), g in seed_df.groupby(["tf", "dna"]):
        d = g["dP"].to_numpy()
        n = len(d)
        mean = d.mean()
        if n > 1:
            se = d.std(ddof=1) / np.sqrt(n)
            tcrit = stats.t.ppf(0.975, n - 1)
            lo, hi = mean - tcrit * se, mean + tcrit * se
            p = stats.ttest_1samp(d, 0).pvalue
        else:
            lo = hi = p = np.nan
        rows.append(dict(tf=tf, dna=dna, n_seeds=n, dP_mean=mean,
                         dP_lo=lo, dP_hi=hi, dP_p=p,
                         dMAE_mean=g["dMAE"].mean(),
                         n_neg_seeds=int((d < 0).sum())))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- mechanism axes
def build_mechanism_table():
    prot = load_protein_axes()
    dna = load_dna_geometry()
    mgw = load_mgwfl()
    eff = summarize_seed_effects(load_seed_effects())

    froz = (eff[eff.dna == "frozen"]
            .set_index("tf")[["dP_mean", "dP_lo", "dP_hi", "dP_p", "n_seeds", "n_neg_seeds", "dMAE_mean"]]
            .add_prefix("froz_"))
    # NOTE: the benchmark CSVs label the DNA conditions "frozen" and "relax"
    # (not "relaxed"). An earlier version of this script filtered on "relaxed" and
    # silently produced an empty relaxed arm. DNA_COND_RELAX is asserted below.
    rel = (eff[eff.dna == DNA_COND_RELAX]
           .set_index("tf")[["dP_mean", "dP_lo", "dP_hi", "dP_p", "n_seeds",
                             "n_neg_seeds", "dMAE_mean"]]
           .add_prefix("rel_"))

    t = (prot.merge(dna, on="tf", how="outer")
              .merge(mgw, on="tf", how="outer")
              .merge(froz.reset_index(), on="tf", how="left")
              .merge(rel.reset_index(), on="tf", how="left"))

    # --- MEASURED induced-fit index -------------------------------------------------
    # How far is the BOUND (crystal) duplex deformed from canonical B-DNA? This is the
    # quantity the curated "minimal / modest / severe" labels were standing in for.
    # Two components, z-scored across pilots and averaged so neither dominates:
    #   bend deviation       : |crystal bend - 0 deg|   (TBP kink, LEF1 bend)
    #   groove deviation     : |crystal minor width - 5.7 A|  (TBP's widened minor groove)
    t["bend_dev"] = (t["cry_bend_uu"] - B_DNA_BEND).abs()
    t["groove_dev"] = (t["cry_minor_w"] - B_DNA_MINOR_W).abs()
    for c in ("bend_dev", "groove_dev"):
        t[f"z_{c}"] = (t[c] - t[c].mean()) / t[c].std(ddof=1)
    t["induced_fit_index"] = t[["z_bend_dev", "z_groove_dev"]].mean(axis=1)

    # --- Conformational-selection index ---------------------------------------------
    # How well does the FREE protein ensemble already reach the bound pose?
    # Small d_min = reaches it. Sign-flipped so larger = more selection-like.
    t["z_reach"] = -(t["d_min"] - t["d_min"].mean()) / t["d_min"].std(ddof=1)
    t["conf_selection_index"] = t["z_reach"]

    # The hypothesis in one number: selection-like minus induced-fit-like.
    t["mechanism_score"] = t["conf_selection_index"] - t["induced_fit_index"]
    return t


# ---------------------------------------------------------------- correlations
AXES = [
    # (column, side, human label, direction hypothesis predicts vs dPearson)
    ("d_min",               "protein", "Reachability d_min (A)",                 "neg"),
    ("spread",              "protein", "Free-state ensemble spread (A)",         "neg"),
    ("rmsf_mean",           "protein", "Ensemble RMSF (A)",                      "neg"),
    ("reach_ratio",         "protein", "Reach ratio",                            "neg"),
    ("cry_bend_uu",         "DNA",     "Crystal DNA bend (deg)",                 "neg"),
    ("bend_dev",            "DNA",     "Crystal bend deviation from B-DNA",      "neg"),
    ("groove_dev",          "DNA",     "Crystal minor-groove deviation from B",  "neg"),
    ("cry_minor_w",         "DNA",     "Crystal minor-groove width (A)",         "neg"),
    ("cry_shortening",      "DNA",     "Crystal duplex shortening",              "neg"),
    ("mgwfl_froz",          "DNA",     "MGW fluctuation, frozen (A)",            "neg"),
    ("mgwfl_relax",         "DNA",     "MGW fluctuation, relaxed (A)",           "neg"),
    ("iface_mgwfl_froz",    "DNA",     "Interface MGW fluctuation, frozen (A)",  "neg"),
    ("iface_mgwfl_relax",   "DNA",     "Interface MGW fluctuation, relaxed (A)", "neg"),
    ("froz_bend_iqr",       "DNA",     "Ensemble DNA-bend IQR, frozen (deg)",    "neg"),
    ("induced_fit_index",   "DNA",     "Induced-fit index (measured)",           "neg"),
    ("conf_selection_index","protein", "Conformational-selection index",         "pos"),
    ("mechanism_score",     "combined","Mechanism score (selection - induced)",  "pos"),
]


def correlate_axes(t, outcome="froz_dP_mean"):
    rows = []
    for col, side, label, pred in AXES:
        sub = t[[col, outcome]].dropna()
        if len(sub) < 4:
            continue
        rho, p = stats.spearmanr(sub[col], sub[outcome])
        r_p, p_p = stats.pearsonr(sub[col], sub[outcome])
        rows.append(dict(axis=col, side=side, label=label, predicted_sign=pred,
                         n=len(sub), spearman_rho=rho, spearman_p=p,
                         pearson_r=r_p, pearson_p=p_p,
                         sign_matches_hypothesis=(rho < 0) == (pred == "neg")))
    return (pd.DataFrame(rows)
            .sort_values("spearman_rho", key=abs, ascending=False)
            .reset_index(drop=True))


# ---------------------------------------------------------------- family analysis
def family_effects_seedlevel():
    """
    Per-family augmentation effect, computed at the seed level.

    perseed_perentry.csv carries a family label per benchmark entry (motif-level
    assignment). For each (pilot, dna, seed, family) we average the per-entry metric
    within the family, then form the seed-paired augmented-minus-baseline delta. The
    test across seeds within a family is therefore on n = seeds, not n = entries.
    """
    pe = pd.read_csv(ROOT / "analysis/figure_scripts/perseed_perentry.csv")
    g = (pe.groupby(["tf", "arm", "dna", "seed", "family"])
           .agg(m_pearsonr=("m_pearsonr", "mean"),
                m_mae=("m_mae", "mean"),
                n_entries=("entry", "size"))
           .reset_index())
    w = g.pivot_table(index=["tf", "dna", "seed", "family", "n_entries"],
                      columns="arm", values=["m_pearsonr", "m_mae"])
    w.columns = [f"{a}_{b}" for a, b in w.columns]
    w = w.reset_index()
    w["dP"] = w["m_pearsonr_augmented"] - w["m_pearsonr_baseline"]
    w["dMAE"] = w["m_mae_augmented"] - w["m_mae_baseline"]
    w = w.dropna(subset=["dP"])

    rows = []
    for (fam, dna), sub in w.groupby(["family", "dna"]):
        # Average across pilots within a seed first, so pilots do not inflate n.
        per_seed = sub.groupby("seed")["dP"].mean()
        d = per_seed.to_numpy()
        n = len(d)
        if n < 2:
            lo = hi = p = np.nan
        else:
            se = d.std(ddof=1) / np.sqrt(n)
            tc = stats.t.ppf(0.975, n - 1)
            lo, hi = d.mean() - tc * se, d.mean() + tc * se
            p = stats.ttest_1samp(d, 0).pvalue
        rows.append(dict(family=fam, dna=dna, n_seeds=n,
                         n_pilots=sub.tf.nunique(),
                         n_entries_median=int(sub.n_entries.median()),
                         dP_mean=d.mean(), dP_lo=lo, dP_hi=hi, dP_p=p,
                         dMAE_mean=sub.groupby("seed")["dMAE"].mean().mean()))
    out = pd.DataFrame(rows)
    # BH-FDR within each DNA condition
    for dna, idx in out.groupby("dna").groups.items():
        p = out.loc[idx, "dP_p"]
        ok = p.notna()
        if ok.sum():
            ranked = p[ok].rank(method="first")
            out.loc[p[ok].index, "dP_q"] = (p[ok] * ok.sum() / ranked).clip(upper=1)
    return out.sort_values(["dna", "dP_mean"], ascending=[True, False]).reset_index(drop=True)


# ---------------------------------------------------------------- regime contrast
def regime_contrast(t):
    """
    Split pilots by the MEASURED induced-fit index at its median and compare the
    augmentation effect between regimes. This is the hypothesis stated as a
    two-group comparison rather than a correlation.
    """
    sub = t.dropna(subset=["induced_fit_index", "froz_dP_mean"]).copy()
    med = sub["induced_fit_index"].median()
    sub["regime"] = np.where(sub["induced_fit_index"] > med,
                             "induced-fit (measured)", "conformational-selection (measured)")
    rows = []
    for reg, g in sub.groupby("regime"):
        rows.append(dict(regime=reg, n_pilots=len(g),
                         pilots=" ".join(sorted(g.tf)),
                         dP_mean=g.froz_dP_mean.mean(),
                         dP_median=g.froz_dP_mean.median(),
                         dMAE_mean=g.froz_dMAE_mean.mean()))
    res = pd.DataFrame(rows)
    a = sub[sub.regime.str.startswith("induced")]["froz_dP_mean"]
    b = sub[sub.regime.str.startswith("conformational")]["froz_dP_mean"]
    u, pu = stats.mannwhitneyu(a, b, alternative="two-sided")
    tt, pt = stats.ttest_ind(a, b, equal_var=False)
    res.attrs["test"] = dict(mannwhitney_u=float(u), mannwhitney_p=float(pu),
                             welch_t=float(tt), welch_p=float(pt),
                             split_value=float(med))
    return res, sub


# ---------------------------------------------------------------- main
def main():
    t = build_mechanism_table()
    t.to_csv(OUT / "data/mechanism_table.csv", index=False)

    corr_froz = correlate_axes(t, outcome="froz_dP_mean").assign(dna="frozen")
    corr_rel = correlate_axes(t, outcome="rel_dP_mean").assign(dna="relax")
    corr = pd.concat([corr_froz, corr_rel], ignore_index=True)
    corr.to_csv(OUT / "data/axis_correlations.csv", index=False)

    fam = family_effects_seedlevel()
    fam.to_csv(OUT / "data/family_effects_seedlevel.csv", index=False)

    reg, tagged = regime_contrast(t)
    reg.to_csv(OUT / "data/regime_contrast.csv", index=False)
    with open(OUT / "data/regime_contrast_test.json", "w") as fh:
        json.dump(reg.attrs["test"], fh, indent=2)

    cols = ["tf", "d_min", "spread", "cry_bend_uu", "cry_minor_w",
            "induced_fit_index", "conf_selection_index", "mechanism_score",
            "froz_dP_mean", "froz_n_seeds"]
    print("=== mechanism table ===")
    print(t[cols].round(3).to_string(index=False))
    for dna_cond in ("frozen", "relax"):
        print(f"\n=== axis vs seed-paired dPearson ({dna_cond}) ===")
        sub = corr[corr.dna == dna_cond]
        print(sub[["label", "side", "n", "spearman_rho", "spearman_p",
                   "sign_matches_hypothesis"]].round(3).to_string(index=False))
    print("\n=== frozen vs relaxed augmentation effect, per pilot ===")
    print(t[["tf", "froz_dP_mean", "froz_dP_p", "rel_dP_mean", "rel_dP_p",
             "froz_n_seeds", "rel_n_seeds"]].dropna(subset=["froz_dP_mean"])
          .round(4).to_string(index=False))
    print("\n=== family effects, seed-level ===")
    print(fam[fam.dna == "frozen"][["family", "n_pilots", "n_seeds",
                                    "n_entries_median", "dP_mean", "dP_lo",
                                    "dP_hi", "dP_p"]].round(4).to_string(index=False))
    print("\n=== regime contrast ===")
    print(reg.round(4).to_string(index=False))
    print(json.dumps(reg.attrs["test"], indent=2))


if __name__ == "__main__":
    main()
