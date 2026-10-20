"""
stats.py — the project's one correct implementation of its recurring statistics.

WHY THIS EXISTS
---------------
Paired-delta testing and FDR correction were implemented twice, in
stage7_eval/evaluate_id_benchmark.py and analysis/analyses/mechanism/
mechanism_analysis.py, with *different* bugs in each:

  * evaluate_id_benchmark.py converted its t-statistic through the NORMAL CDF
    (`2*(1-0.5*(1+erf(|t|/sqrt(2))))`). At the n=5 that most on-target subsets
    run, that understates p by 2.6x (t=2.0) to 14.8x (t=3.0). Every
    "significant" subset p-value in the stage-7 logs is overstated.
  * mechanism_analysis.py used the correct t-distribution but computed BH-FDR
    without the step-up cumulative minimum, so q-values can be non-monotonic in
    p — a more-significant test can be rejected while a weaker one passes.
  * evaluate_id_benchmark.py also re-seeded its bootstrap RNG *inside* the
    per-metric loop, so all six metrics resampled on identical indices and
    their CIs were perfectly correlated rather than independent.

All three are fixed here. Import this rather than reimplementing:

    import sys, os
    sys.path.insert(0, os.path.join(TFCONF, "analysis", "common"))
    from stats import paired_delta, bh_fdr, effective_n

Findings F9 and F12 of docs/codebase_review_2026-10-15.md.
"""
import numpy as np

try:
    from scipy import stats as _st
    _HAVE_SCIPY = True
except ImportError:  # pragma: no cover
    _HAVE_SCIPY = False


def _t_sf_two_sided(t, df):
    """Two-sided p from the t-distribution. Falls back loudly, never silently."""
    if df <= 0:
        return float("nan")
    if _HAVE_SCIPY:
        return float(2 * _st.t.sf(abs(t), df))
    raise ImportError(
        "scipy is required for correct t-distribution p-values. Do NOT "
        "substitute the normal approximation: at n=5 it understates p by "
        "2.6-15x, which is exactly the bug this module exists to fix.")


def paired_delta(baseline, augmented, n_boot=10000, seed=42, ci=95.0):
    """Paired delta (augmented - baseline) with a correct p-value and bootstrap CI.

    Parameters
    ----------
    baseline, augmented : 1-D array-likes of equal length, paired elementwise.
        NaN in either member drops that pair.
    n_boot : bootstrap resamples for the CI of the mean delta.
    seed : RNG seed. Seeded ONCE per call, so separate metrics get independent
        resamples (the old code re-seeded per metric inside the loop).
    ci : central CI width in percent.

    Returns
    -------
    dict with n_pairs, mean_delta, median_delta, sd, ci_lo, ci_hi, t_stat,
    p_value (two-sided, t-distribution with df = n-1), and cohens_dz.
    Returns n_pairs and NaNs when fewer than 2 usable pairs remain.
    """
    b = np.asarray(baseline, dtype=float)
    a = np.asarray(augmented, dtype=float)
    if b.shape != a.shape:
        raise ValueError(f"baseline {b.shape} and augmented {a.shape} must pair "
                         f"elementwise")
    ok = ~(np.isnan(b) | np.isnan(a))
    b, a = b[ok], a[ok]
    n = b.size
    nan = float("nan")
    if n < 2:
        return dict(n_pairs=int(n), mean_delta=nan, median_delta=nan, sd=nan,
                    ci_lo=nan, ci_hi=nan, t_stat=nan, p_value=nan, cohens_dz=nan)

    d = a - b
    mean = float(d.mean())
    sd = float(d.std(ddof=1))

    # Guard against degenerate spread. If the deltas are constant to within
    # float noise, the t-statistic explodes and reports a meaningless p like
    # 1e-46 with dz ~ 1e15. That is not a finding, it is a numerically
    # degenerate input (e.g. every entry improved by exactly the same amount,
    # which in practice means the "replicates" are not independent).
    scale = max(abs(mean), np.abs(d).max(), 1.0)
    degenerate = sd <= 1e-12 * scale

    if degenerate:
        t = float("nan")
        p = float("nan")
        dz = float("nan")
    else:
        t = mean / (sd / np.sqrt(n))
        p = _t_sf_two_sided(t, n - 1)
        dz = mean / sd

    rng = np.random.default_rng(seed)
    idx = rng.integers(0, n, size=(n_boot, n))
    boot = d[idx].mean(axis=1)
    lo, hi = np.percentile(boot, [(100 - ci) / 2, 100 - (100 - ci) / 2])

    return dict(n_pairs=int(n), mean_delta=mean, median_delta=float(np.median(d)),
                sd=sd, ci_lo=float(lo), ci_hi=float(hi), t_stat=float(t),
                p_value=float(p), cohens_dz=float(dz))


def bh_fdr(pvals):
    """Benjamini-Hochberg q-values, WITH the step-up monotonicity enforcement.

    The enforcement is the part the previous implementation omitted:
    q must be non-decreasing in p, which requires a reverse cumulative minimum
    over the p-ordered values. Without it:

        p=0.001 -> q=0.0050     p=0.001 -> q=0.0050
        p=0.030 -> q=0.0750     p=0.030 -> q=0.0667   <- correct
        p=0.040 -> q=0.0667     p=0.040 -> q=0.0667
                   (as coded)              (BH)

    NaN p-values are passed through as NaN and excluded from m.
    """
    p = np.asarray(pvals, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    m = int(ok.sum())
    if m == 0:
        return out
    pv = p[ok]
    order = np.argsort(pv, kind="stable")
    ranked = pv[order] * m / np.arange(1, m + 1)
    # step-up: enforce monotonicity from the largest p downward
    q_sorted = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(m)
    q[order] = np.clip(q_sorted, 0.0, 1.0)
    out[ok] = q
    return out


def effective_n(entries):
    """Number of independent structures behind a list of benchmark entries.

    Entries are per-CHAIN (`<pdb>_<chain>_<pwm>.npz`), so four chains of one
    crystal are four entries but ONE structure. Treating them as independent is
    pseudoreplication: it inflates n and deflates p. The retired IRF pilot's
    "n=4" on-target subset was four chains of a single PDB (7jm4) — effective
    n=1. See docs/retired_pilots.md.

    Returns (n_entries, n_structures). Report both whenever you report a
    paired test over benchmark entries.
    """
    entries = list(entries)
    pdbs = {str(e).split("_")[0].lower() for e in entries if e}
    return len(entries), len(pdbs)


def describe_paired(res, n_struct=None, label=""):
    """One-line human summary of a paired_delta result, with honest caveats."""
    if res["n_pairs"] < 2:
        return f"{label}: n={res['n_pairs']} — no paired statistic possible"
    s = (f"{label}: Δ={res['mean_delta']:+.4f} "
         f"[{res['ci_lo']:+.4f}, {res['ci_hi']:+.4f}] "
         f"p={res['p_value']:.3g} (n={res['n_pairs']}, dz={res['cohens_dz']:+.2f})")
    if n_struct is not None and n_struct < res["n_pairs"]:
        s += (f"  ⚠ pseudoreplicated: {res['n_pairs']} entries span only "
              f"{n_struct} structure(s); effective n={n_struct}")
    if res["n_pairs"] < 5:
        s += "  ⚠ n<5"
    return s
