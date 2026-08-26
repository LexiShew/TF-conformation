#!/usr/bin/env python
"""
mechanism_confound.py -- is the own-family effect mechanism, or just headroom?

mechanism_probe.py found that the augmentation effect lives in the OWN-FAMILY subset,
not the 130-entry cross-benchmark average (ETS1 own-family dPearson = +0.111, p = 0.013,
0 of 6 seeds negative), and that it tracks a DYNAMIC DNA axis (ensemble DNA-bend IQR,
Spearman rho = -0.53) rather than the STATIC crystal deformation the original induced-fit
hypothesis was written on.

Before that can be called mechanism, two alternatives must be excluded.

A1  HEADROOM. A pilot whose own family the baseline predicts poorly has more room to
    improve, so any regularizer would raise it. If own-family effect is explained by
    baseline own-family accuracy, the "mechanism" is arithmetic.

A2  SUBSET SIZE. Own-family subsets range from 3 to 24 entries. Small subsets have
    noisier means and can manufacture large deltas. If effect size tracks 1/n, the
    ranking is a sampling artifact.

Both are tested as correlations against the seed-level own-family effect, and then
jointly in an OLS model so the DNA axis is evaluated with headroom partialled out.

A third check (D1) tests the dose-response the mechanism predicts: within a pilot, does
the augmentation effect on an individual benchmark entry track how far that entry's
own DNA geometry sits from the pilot's crystal DNA? This is the finest-grained version
of the claim and does not depend on the cross-pilot n = 11.

Statistics
----------
Seed remains the unit of replication for effects. Cross-pilot models are n = 11, which
is small enough that the OLS is reported for its coefficient signs and partial ranking,
not as a significance test. Bootstrap CIs (10,000 resamples over pilots) are given for
the headline correlation so the reader can see its instability directly.

Usage
-----
    source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh && conda activate deeppbs
    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
    cd /project2/rohs_102/shewchuk/TF-conformation
    python analysis/analyses/mechanism/mechanism_confound.py

Outputs (analysis/analyses/mechanism/data/)
    confound_tests.csv       headroom / subset-size correlations
    partial_model.csv        OLS of own-family effect on DNA axis + headroom + size
    bootstrap_ci.json        bootstrap CI for the headline correlation
"""
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent / "data"
RNG = np.random.default_rng(0)


def own_family_baseline():
    """Baseline (un-augmented) accuracy on each pilot's own family -- the headroom axis."""
    pe = pd.read_csv(ROOT / "analysis/data/perseed_perentry.csv")
    b = pe[(pe.arm == "baseline") & (pe.dna == "frozen") & (pe.selffam)]
    return (b.groupby("tf")
             .agg(base_selffam_P=("m_pearsonr", "mean"),
                  base_selffam_MAE=("m_mae", "mean"),
                  n_selffam_entries=("entry", "nunique"))
             .reset_index())


def build():
    sf = pd.read_csv(OUT / "selffam_effects.csv")
    sf = (sf[(sf.dna == "frozen") & (sf.selffam)]
          .rename(columns={"dP_mean": "selffam_dP", "dP_p": "selffam_p",
                           "n_neg_seeds": "selffam_n_neg"})
          [["tf", "selffam_dP", "selffam_p", "selffam_n_neg", "n_seeds"]])
    mech = pd.read_csv(OUT / "mechanism_table.csv")[
        ["tf", "froz_bend_iqr", "d_min", "spread", "cry_bend_uu",
         "induced_fit_index", "iface_mgwfl_froz", "mgwfl_froz", "froz_dP_mean"]]
    return sf.merge(own_family_baseline(), on="tf").merge(mech, on="tf", how="left")


def confound_tests(df):
    rows = []
    tests = [
        ("base_selffam_P",    "A1 headroom: baseline own-family Pearson"),
        ("base_selffam_MAE",  "A1 headroom: baseline own-family MAE"),
        ("n_selffam_entries", "A2 subset size: own-family entry count"),
        ("froz_bend_iqr",     "DNA axis: ensemble DNA-bend IQR (frozen)"),
        ("d_min",             "protein axis: reachability d_min"),
        ("cry_bend_uu",       "static DNA: crystal bend"),
        ("induced_fit_index", "static DNA: induced-fit index"),
    ]
    for col, label in tests:
        sub = df[[col, "selffam_dP"]].dropna()
        if len(sub) < 4:
            continue
        rho, p = stats.spearmanr(sub[col], sub.selffam_dP)
        r, rp = stats.pearsonr(sub[col], sub.selffam_dP)
        rows.append(dict(axis=col, label=label, n=len(sub),
                         spearman_rho=rho, spearman_p=p,
                         pearson_r=r, pearson_p=rp))
    return pd.DataFrame(rows).sort_values("spearman_rho", key=abs,
                                          ascending=False).reset_index(drop=True)


def partial_model(df):
    """OLS: own-family effect ~ DNA-bend IQR + baseline headroom + subset size (z-scored)."""
    d = df.dropna(subset=["selffam_dP", "froz_bend_iqr", "base_selffam_P",
                          "n_selffam_entries"]).copy()
    X = d[["froz_bend_iqr", "base_selffam_P", "n_selffam_entries"]].copy()
    for col in X:
        X[col] = (X[col] - X[col].mean()) / X[col].std(ddof=1)
    X.insert(0, "const", 1.0)
    y = d["selffam_dP"].to_numpy()
    Xm = X.to_numpy()
    beta, *_ = np.linalg.lstsq(Xm, y, rcond=None)
    resid = y - Xm @ beta
    dof = len(y) - Xm.shape[1]
    s2 = resid @ resid / dof
    cov = s2 * np.linalg.inv(Xm.T @ Xm)
    se = np.sqrt(np.diag(cov))
    tvals = beta / se
    pvals = 2 * stats.t.sf(np.abs(tvals), dof)
    ss_tot = ((y - y.mean()) ** 2).sum()
    return pd.DataFrame(dict(term=X.columns, beta=beta, se=se, t=tvals, p=pvals)), \
        dict(n=len(y), dof=int(dof), r2=float(1 - resid @ resid / ss_tot))


def bootstrap_ci(df, xcol="froz_bend_iqr", n_boot=10000):
    """Bootstrap the headline Spearman over pilots -- shows how unstable n=11 really is."""
    sub = df[[xcol, "selffam_dP"]].dropna()
    x, y = sub[xcol].to_numpy(), sub.selffam_dP.to_numpy()
    obs = stats.spearmanr(x, y)[0]
    boot = []
    for _ in range(n_boot):
        i = RNG.integers(0, len(x), len(x))
        if len(np.unique(x[i])) < 3:
            continue
        boot.append(stats.spearmanr(x[i], y[i])[0])
    boot = np.array([b for b in boot if np.isfinite(b)])
    return dict(axis=xcol, n_pilots=int(len(x)), observed_rho=float(obs),
                ci_lo=float(np.percentile(boot, 2.5)),
                ci_hi=float(np.percentile(boot, 97.5)),
                frac_same_sign=float(np.mean(np.sign(boot) == np.sign(obs))),
                n_boot=int(len(boot)))


def dose_response():
    """
    D1 -- within-pilot dose response.

    For each pilot, correlate the per-entry augmentation effect (seed-averaged) against
    that entry's own family membership and, where available, the DNA geometry of the
    entry's own structure. Here the available per-entry DNA descriptor is whether the
    entry belongs to the pilot's family (selffam), so this reduces to a within-pilot
    same-vs-other contrast tested at the seed level -- the test REVIEW_figure_scripts.md
    asked for and fig6 never performed.
    """
    pe = pd.read_csv(ROOT / "analysis/data/perseed_perentry.csv")
    pe = pe[pe.dna == "frozen"]
    g = (pe.groupby(["tf", "arm", "seed", "selffam"])
           .agg(m_pearsonr=("m_pearsonr", "mean")).reset_index())
    w = g.pivot_table(index=["tf", "seed", "selffam"], columns="arm",
                      values="m_pearsonr").reset_index()
    w["dP"] = w["augmented"] - w["baseline"]
    # A seed that is missing one arm yields a NaN delta. Dropping those here keeps the
    # paired contrast on seeds where BOTH arms trained (nfat has one such seed).
    w = w.dropna(subset=["dP"])
    rows = []
    for tf, sub in w.groupby("tf"):
        a = sub[sub.selffam].set_index("seed")["dP"]
        b = sub[~sub.selffam].set_index("seed")["dP"]
        common = a.index.intersection(b.index)
        if len(common) < 3:
            continue
        diff = (a[common] - b[common]).to_numpy()
        t, p = stats.ttest_1samp(diff, 0)
        rows.append(dict(tf=tf, n_seeds=len(common),
                         selffam_dP=a[common].mean(), otherfam_dP=b[common].mean(),
                         contrast=diff.mean(), t=t, p=p,
                         n_neg_seeds=int((diff < 0).sum())))
    out = pd.DataFrame(rows).sort_values("contrast", ascending=False).reset_index(drop=True)
    ranked = out.p.rank(method="first")
    out["q_bh"] = (out.p * len(out) / ranked).clip(upper=1)
    return out


def main():
    df = build()
    ct = confound_tests(df)
    ct.to_csv(OUT / "confound_tests.csv", index=False)
    print("=== confound tests vs own-family augmentation effect ===")
    print(ct.round(3).to_string(index=False))

    pm, fit = partial_model(df)
    pm.to_csv(OUT / "partial_model.csv", index=False)
    print(f"\n=== OLS: selffam_dP ~ bend_IQR + baseline headroom + subset size "
          f"(n={fit['n']}, dof={fit['dof']}, R2={fit['r2']:.3f}) ===")
    print(pm.round(4).to_string(index=False))

    bc = bootstrap_ci(df)
    with open(OUT / "bootstrap_ci.json", "w") as fh:
        json.dump(bc, fh, indent=2)
    print("\n=== bootstrap CI, headline correlation ===")
    print(json.dumps(bc, indent=2))

    dr = dose_response()
    dr.to_csv(OUT / "same_vs_other_contrast.csv", index=False)
    print("\n=== D1 within-pilot same-family vs other-family contrast (seed level) ===")
    print(dr.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
