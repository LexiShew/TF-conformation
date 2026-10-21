#!/usr/bin/env python
"""Decompose, localize and rebuild the DNA bend fluctuation from local
base-pair-step geometry -- the three analyses that withdrew the August
bend-IQR mechanism claim.

Background. The August claim was that the augmentation effect tracks DNA bend
*fluctuation* (IQR of pyCurves `bend_uu` across an ensemble) rather than static
crystal bend, at rho = -0.87. `bend_uu` is a terminal axis-vector angle: it
looks at the two ends of the duplex and ignores everything between them. These
three analyses ask whether the quantity it summarizes behaves the same way.

  1. DECOMPOSE   which local degree of freedom carries the bend fluctuation.
                 Re-accumulate the duplex bend with roll, tilt or twist held at
                 its per-step mean across the ensemble, and compare the
                 resulting IQR with the all-varying case.
  2. LOCALIZE    where along the duplex the fluctuation sits: per-step IQR of
                 roll and tilt as a function of normalized position.
  3. REBUILD     accumulate the local steps into a duplex bend four ways
                 (complete, central window, each end, end-to-end frame
                 orientation) and correlate each against both `bend_uu` IQR and
                 the augmentation effect.

Two independent accumulation methods are used, because agreeing on a null is
worth more than either alone:
  * twist-aware vector sum -- each step contributes a bending vector
    (roll, tilt) rotated by the twist accumulated before it, then summed. This
    is the standard twist-corrected curvature accumulation.
  * rotation-chain orientation -- compose each step's rotation
    (tilt, roll, twist as a rotation vector in the local frame) and measure the
    angle between the first and last helical-axis direction. No small-angle or
    planarity assumption.

Inputs (all already on disk; nothing is recomputed)
  analysis/analyses/dna_relax/pycurves/<pilot>/<cond>_state_NNN_curvesplus.json
      per-state pyCurves output, `dataframes.local_inter_base_pair` giving
      per-step shift/slide/rise/tilt/roll/twist. Gitignored as regenerable
      (~331 MB) but retained.
  analysis/analyses/dna_relax/data/pycurves_all_perstructure.csv
      whole-duplex `bend_uu` per state -- the quantity under scrutiny.
  analysis/analyses/mechanism/data/mechanism_summary_table.csv
      per-pilot augmentation effect (`ownfam_dPearson_frozen`) and the original
      `bend_IQR_deg` axis, so the August correlation can be reproduced here and
      the reconstructions tested against the same effect.

Outputs
  data/bend_recon_perstate.csv     one row per pilot x cond x state (checkpoint)
  data/bend_decomp_perpilot.csv    IQRs, component decomposition, per pilot x cond
  data/bend_localization.csv       per-step roll/tilt/twist IQR along the duplex
  data/bend_axis_correlations.csv  every axis vs bend_uu IQR and vs the effect
  figures/B1_bend_decomposition.png
  figures/B2_bend_localization.png
  figures/B3_bend_reconstruction.png

CAVEAT ON PROVENANCE. The August numbers (rho = -0.87, and rho = +0.05 between
`bend_uu` IQR and the complete reconstruction) were computed outside the repo
and no script or table for them survives. This is a RE-DERIVATION, not a
reproduction: where it disagrees, these values are the ones with a committed
script behind them.

Usage (parse is cached -- pass --reparse to redo it)
    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
    cd analysis/analyses/mechanism
    /project2/rohs_102/shewchuk/conda/envs/deeppbs/bin/python \
        make_bend_decomposition.py [--reparse] [--cond frozen]

Parse only -- no JAX, no pyCurves re-run, so a login node is fine (the thread
pin is still required for numpy there). ~4800 JSON files, a few minutes.
"""
import argparse
import glob
import json
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.transform import Rotation as Rot

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
RELAX = os.path.normpath(os.path.join(HERE, "..", "dna_relax"))
PYC = os.path.join(RELAX, "pycurves")
DATA = os.path.join(HERE, "data")
FIGS = os.path.join(HERE, "figures")
sys.path.insert(0, REPO)
from palette import GREY, GREY_R, TEAL, TEAL_R, ALARM  # noqa: E402

MG = GREY
STEP_COLS = ["tilt", "roll", "twist"]


# --------------------------------------------------------------------------- #
# accumulation
# --------------------------------------------------------------------------- #
def bend_vector_sum(tilt, roll, twist):
    """Twist-aware accumulation of the local bending vectors, in degrees.

    Step i bends by (roll_i, tilt_i) in its own frame; that frame is rotated
    about the helical axis by the twist accumulated over steps 0..i-1, so the
    contributions are rotated into a common frame before being summed. Ignoring
    the twist is what makes a naive sum of |roll| overestimate curvature.
    """
    phi = np.radians(np.concatenate([[0.0], np.cumsum(twist)[:-1]]))
    bx = np.sum(roll * np.cos(phi) - tilt * np.sin(phi))
    by = np.sum(roll * np.sin(phi) + tilt * np.cos(phi))
    return float(np.hypot(bx, by))


def frame_orientation(tilt, roll, twist):
    """Angle between the first and last helical-axis direction, in degrees.

    Composes each step's rotation exactly (rotation vector (tilt, roll, twist)
    in the local frame, x = short axis, y = long axis, z = helical axis), so
    unlike the vector sum this makes no small-angle or planarity assumption.
    """
    R = Rot.identity()
    for t, r, w in zip(tilt, roll, twist):
        R = R * Rot.from_rotvec(np.radians([t, r, w]))
    z = R.apply([0.0, 0.0, 1.0])
    return float(np.degrees(np.arccos(np.clip(z[2], -1.0, 1.0))))


def windows(n):
    """Central half, 5' quarter and 3' quarter of an n-step duplex."""
    q = max(1, n // 4)
    return {"central": slice(q, n - q), "end5": slice(0, q), "end3": slice(n - q, n)}


# --------------------------------------------------------------------------- #
# parse
# --------------------------------------------------------------------------- #
def read_steps(path):
    """Per-step tilt/roll/twist for one state, ordered 5'->3'.

    Returns None when the state has no base-pair-level step table. That happens
    for irf, whose palindromic ISRE duplex pyCurves does not pair: those files
    carry only `local_inter_base` (per strand). Averaging the two strands would
    need the antiparallel tilt sign convention to be right, so irf is excluded
    and reported rather than guessed at.
    """
    try:
        with open(path) as fh:
            d = json.load(fh)
    except (ValueError, OSError):
        return None
    rows = d.get("dataframes", {}).get("local_inter_base_pair")
    if not rows:
        return None
    df = pd.DataFrame(rows)
    key = "level" if "level" in df.columns else "sequence_index"
    df = df.sort_values(key)
    if not set(STEP_COLS) <= set(df.columns):
        return None
    df["_lvl"] = pd.to_numeric(df[key], errors="coerce")
    out = df[["_lvl"] + STEP_COLS].apply(pd.to_numeric, errors="coerce").dropna()
    if not len(out):
        return None
    return out["_lvl"].to_numpy(float), out[STEP_COLS].to_numpy(float)


def parse_all():
    """Walk every retained per-state JSON once; return per-state reconstructions
    plus the raw per-step arrays needed for the component decomposition."""
    recs, raw, skipped = [], {}, {}
    pilots = sorted(p for p in os.listdir(PYC) if os.path.isdir(os.path.join(PYC, p)))
    for pilot in pilots:
        for path in sorted(glob.glob(os.path.join(PYC, pilot, "*_curvesplus.json"))):
            m = re.match(r"(frozen|relaxed)_state_(\d+)_curvesplus\.json",
                         os.path.basename(path))
            if not m:
                continue
            cond, state = m.group(1), int(m.group(2))
            got = read_steps(path)
            if got is None or len(got[1]) < 6:
                skipped[pilot] = skipped.get(pilot, 0) + 1
                continue
            lvl, arr = got
            tilt, roll, twist = arr[:, 0], arr[:, 1], arr[:, 2]
            n = len(tilt)
            w = windows(n)
            recs.append(dict(
                tf=pilot, cond=cond, state=state, n_steps=n,
                recon_complete=bend_vector_sum(tilt, roll, twist),
                recon_central=bend_vector_sum(*[x[w["central"]] for x in (tilt, roll, twist)]),
                recon_end5=bend_vector_sum(*[x[w["end5"]] for x in (tilt, roll, twist)]),
                recon_end3=bend_vector_sum(*[x[w["end3"]] for x in (tilt, roll, twist)]),
                recon_orient=frame_orientation(tilt, roll, twist),
                mean_abs_roll=float(np.mean(np.abs(roll))),
                mean_abs_tilt=float(np.mean(np.abs(tilt))),
                mean_twist=float(np.mean(twist))))
            raw.setdefault((pilot, cond), []).append((state, lvl, arr))
        print("  parsed %-10s %d states (%d without a base-pair step table)"
              % (pilot, sum(1 for r in recs if r["tf"] == pilot),
                 skipped.get(pilot, 0)), flush=True)
    dropped = {k: v for k, v in skipped.items()
               if not any(r["tf"] == k for r in recs)}
    if dropped:
        print("EXCLUDED (no base-pair-level steps in any state): %s"
              % ", ".join("%s (%d files)" % kv for kv in sorted(dropped.items())),
              flush=True)
    return pd.DataFrame(recs), raw


def align(mats, min_frac=.90):
    """Put a pilot x cond ensemble on ONE set of base-pair levels.

    pyCurves does not detect the same number of base pairs in every state, so
    positional index i is not the same step across states. Align on the level
    label instead: keep the levels present in >= min_frac of states, then keep
    only the states that have all of them. Returns (states, levels, M) with M
    shaped states x levels x (tilt, roll, twist).
    """
    from collections import Counter
    seen = Counter()
    for _, lvl, _ in mats:
        seen.update(set(lvl.tolist()))
    keep = sorted(l for l, c in seen.items() if c >= min_frac * len(mats))
    if len(keep) < 6:
        return None
    kept_states, rows = [], []
    for state, lvl, arr in mats:
        idx = {l: j for j, l in enumerate(lvl.tolist())}
        if not all(l in idx for l in keep):
            continue
        kept_states.append(state)
        rows.append(np.stack([arr[idx[l]] for l in keep]))
    if len(rows) < 10:
        return None
    return kept_states, keep, np.stack(rows)


def decompose(raw):
    """Which local degree of freedom carries the bend fluctuation?

    Re-accumulate each state with one component varying and the other two held
    at their per-step ensemble mean. If the duplex bend fluctuation is carried
    by roll, the roll-only IQR recovers the all-varying IQR.
    """
    out = []
    for (pilot, cond), mats in raw.items():
        got = align(mats)
        if got is None:
            continue
        _, keep, M = got                     # states x levels x (tilt, roll, twist)
        mean = M.mean(axis=0)
        full = np.array([bend_vector_sum(s[:, 0], s[:, 1], s[:, 2]) for s in M])
        row = dict(tf=pilot, cond=cond, n_states=len(M), n_steps=M.shape[1],
                   n_states_raw=len(mats), iqr_all=float(stats.iqr(full)))
        for i, name in enumerate(STEP_COLS):
            V = np.repeat(mean[None], len(M), axis=0).copy()
            V[:, :, i] = M[:, :, i]          # only this component varies
            vals = np.array([bend_vector_sum(s[:, 0], s[:, 1], s[:, 2]) for s in V])
            row[f"iqr_{name}_only"] = float(stats.iqr(vals))
        out.append(row)
    D = pd.DataFrame(out)
    for name in STEP_COLS:                   # share of the all-varying IQR
        D[f"frac_{name}"] = D[f"iqr_{name}_only"] / D.iqr_all
    return D


def localize(raw):
    """Per-step IQR of each component, at normalized position along the duplex."""
    out = []
    for (pilot, cond), mats in raw.items():
        got = align(mats)
        if got is None:
            continue
        _, keep, M = got
        n = M.shape[1]
        for j in range(n):
            out.append(dict(tf=pilot, cond=cond, step=j, n_steps=n,
                            pos=(j + .5) / n,
                            **{f"iqr_{c}": float(stats.iqr(M[:, j, i]))
                               for i, c in enumerate(STEP_COLS)},
                            **{f"med_{c}": float(np.median(M[:, j, i]))
                               for i, c in enumerate(STEP_COLS)}))
    return pd.DataFrame(out)


# --------------------------------------------------------------------------- #
# figures
# --------------------------------------------------------------------------- #
def style():
    plt.rcParams.update({
        "font.size": 7, "axes.titlesize": 8.5, "axes.labelsize": 7.4,
        "xtick.labelsize": 6.6, "ytick.labelsize": 6.6, "legend.fontsize": 6.4,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": .7, "figure.dpi": 300, "savefig.dpi": 300})


def fig_B1(D, cond, out):
    S = D[D.cond == cond].sort_values("iqr_all", ascending=False)
    fig, ax = plt.subplots(figsize=(6.6, 3.7))
    x = np.arange(len(S))
    w = .26
    for k, (name, col) in enumerate(zip(STEP_COLS, [TEAL_R[0], TEAL, GREY_R[1]])):
        ax.bar(x + (k - 1) * w, S[f"iqr_{name}_only"], width=w, color=col,
               label=f"{name} only", zorder=3)
    ax.plot(x, S.iqr_all, "k_", ms=13, mew=1.6, label="all three varying", zorder=4)
    ax.set_xticks(x)
    ax.set_xticklabels(S.tf, rotation=38, ha="right")
    ax.set_ylabel("IQR of reconstructed duplex bend (°)")
    ax.set_title("Which local degree of freedom carries the bend fluctuation?\n"
                 "one component varies, the other two held at their per-step mean")
    ax.legend(frameon=False, ncol=4, loc="upper right")
    ax.text(.5, -.42, "%s ensembles · twist-aware accumulation of every base-pair step"
            % cond, transform=ax.transAxes, ha="center", fontsize=6, color=MG)
    fig.savefig(os.path.join(out, "B1_bend_decomposition.png"), bbox_inches="tight")
    plt.close(fig)


def fig_B2(L, cond, out):
    S = L[L.cond == cond]
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.3))
    for ax, comp, lab in zip(axes, ["roll", "tilt"], ["roll", "tilt"]):
        for tf, g in S.groupby("tf"):
            g = g.sort_values("pos")
            ax.plot(g.pos, g[f"iqr_{comp}"], lw=.8, color=TEAL, alpha=.42, zorder=2)
        med = S.groupby(pd.cut(S.pos, np.linspace(0, 1, 11)), observed=True)[
            f"iqr_{comp}"].median()
        ctr = [iv.mid for iv in med.index]
        ax.plot(ctr, med.values, lw=2.1, color=GREY_R[2], zorder=4,
                label="median over pilots")
        ax.axvspan(.25, .75, color=GREY_R[0], alpha=.22, zorder=1)
        ax.text(.5, .02, "central half", transform=ax.transAxes, ha="center",
                fontsize=6, color=MG)
        ax.set_xlabel("normalized position along the duplex (5' → 3')")
        ax.set_ylabel(f"per-step IQR of {lab} (°)")
        ax.set_title(f"{lab}")
    axes[0].legend(frameon=False, loc="upper center")
    fig.suptitle("The fluctuation is centrally weighted — and the ends, which are "
                 "all bend_uu looks at, are the quietest part of the duplex",
                 fontsize=9, y=1.03)
    fig.text(.5, -.08, "%s ensembles · one line per pilot. Median per-step roll IQR runs "
             "~1.5° at the termini against ~4.5° mid-duplex." % cond,
             ha="center", fontsize=6, color=MG)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "B2_bend_localization.png"), bbox_inches="tight")
    plt.close(fig)


def fig_B3(P, out, cond, eff="effect_relaxed"):
    """Beats 1, 3 and 4 -- but with BOTH reconstructions, which is the finding.

    The complete twist-aware accumulation of every local step tracks neither the
    effect nor bend_uu. The rotation-chain end-to-end orientation -- built from
    the same local steps, by an independent method -- tracks both. So bend_uu is
    specifically an END-TO-END orientation measure, not accumulated local
    curvature. Panels 2 and 3 plot the axis as a RANK so two quantities on
    different degree scales can share an x-axis; Spearman is rank-based, so
    nothing is lost.
    """
    rank = lambda v: stats.rankdata(v)
    series = [("recon_complete_iqr", "complete accumulation\nof every step", TEAL),
              ("recon_orient_iqr", "end-to-end frame\norientation", GREY_R[2])]
    elab = "augmentation effect (own-family Δr, %s arm)" % eff.split("_")[1]
    fig, axes = plt.subplots(1, 3, figsize=(9.3, 3.5))

    S = P.dropna(subset=["bend_uu_iqr", eff])
    rho, pv = stats.spearmanr(S.bend_uu_iqr, S[eff])
    axes[0].scatter(S.bend_uu_iqr, S[eff], s=40, color=TEAL, alpha=.88, zorder=3)
    for _, r in S.iterrows():
        axes[0].annotate(r.tf, (r.bend_uu_iqr, r[eff]), textcoords="offset points",
                         xytext=(5, 3), fontsize=5.4, color=MG)
    axes[0].axhline(0, color=MG, lw=.8, ls="--", zorder=1)
    axes[0].set_title("the August axis reproduces\nρ = %+.2f, p = %.3f  (n = %d)"
                      % (rho, pv, len(S)))
    axes[0].set_xlabel("IQR of bend_uu (°)\nterminal axis-vector angle")
    axes[0].set_ylabel(elab)

    for col, lab, c in series:
        T = P.dropna(subset=[col, eff])
        rho, pv = stats.spearmanr(T[col], T[eff])
        axes[1].scatter(rank(T[col]), T[eff], s=40, color=c, alpha=.88, zorder=3,
                        label="%s\nρ = %+.2f, p = %.3f" % (lab, rho, pv))
    axes[1].axhline(0, color=MG, lw=.8, ls="--", zorder=1)
    axes[1].set_title("rebuilt from the local steps,\ntwo ways")
    axes[1].set_xlabel("rank of the reconstruction's IQR")
    axes[1].set_ylabel(elab)
    axes[1].legend(frameon=False, loc="upper right", fontsize=5.6, labelspacing=.9)

    for col, lab, c in series:
        T = P.dropna(subset=[col, "bend_uu_iqr"])
        rho, pv = stats.spearmanr(T.bend_uu_iqr, T[col])
        axes[2].scatter(rank(T.bend_uu_iqr), rank(T[col]), s=40, color=c, alpha=.88,
                        zorder=3, label="%s\nρ = %+.2f, p = %.4f" % (lab, rho, pv))
    axes[2].set_title("which of them IS bend_uu\na summary of?")
    axes[2].set_xlabel("rank of bend_uu IQR")
    axes[2].set_ylabel("rank of the reconstruction's IQR")
    axes[2].legend(frameon=False, loc="upper left", fontsize=5.6, labelspacing=.9)

    fig.text(.5, -.13, "Axis from the %s ensembles; effect from the %s arm — the pairing FINDINGS.md "
             "records for the committed claim (ρ = −0.76, p = 0.006, n = 11).\n"
             "Re-derived %s. Presentation v4's ρ = −0.87 and ρ = +0.05 have no surviving script; these do."
             % (cond, eff.split("_")[1], pd.Timestamp.today().date()),
             ha="center", fontsize=6, color=ALARM, linespacing=1.5)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "B3_bend_reconstruction.png"), bbox_inches="tight")
    plt.close(fig)


HOLO_OR_GREY = GREY_R[2]


# --------------------------------------------------------------------------- #
def per_pilot_table(R, D, cond):
    """Join every candidate axis onto the augmentation effect, per pilot."""
    S = R[R.cond == cond]
    axes = {f"{c}_iqr": stats.iqr for c in
            ["recon_complete", "recon_central", "recon_end5", "recon_end3", "recon_orient"]}
    P = pd.DataFrame({name: S.groupby("tf")[name.replace("_iqr", "")].apply(fn)
                      for name, fn in axes.items()})
    bu = pd.read_csv(os.path.join(RELAX, "data", "pycurves_all_perstructure.csv"))
    bu = bu[bu.cond == cond]
    P["bend_uu_iqr"] = bu.groupby("tf").bend_uu.apply(stats.iqr)
    P["n_states"] = S.groupby("tf").state.count()
    M = pd.read_csv(os.path.join(DATA, "mechanism_summary_table.csv"))
    M["tf"] = M.TF.str.lower().str.strip()
    # The axis arm and the effect arm are DIFFERENT choices, and FINDINGS.md
    # (lines 266-267) is explicit that the committed claim is the frozen-arm
    # bend axis against the RELAXED-arm augmentation effect: rho = -0.76,
    # bootstrap CI [-0.96, -0.28], p = 0.006, n = 11; the frozen-arm effect
    # gives only -0.53, p = 0.096. So both effects are carried through and the
    # correlation table reports each separately rather than silently pairing
    # each axis with its own arm.
    P = P.join(M.set_index("tf")[["ownfam_dPearson_frozen",
                                  "ownfam_dPearson_relaxed", "bend_IQR_deg"]]
               .rename(columns={"ownfam_dPearson_frozen": "effect_frozen",
                                "ownfam_dPearson_relaxed": "effect_relaxed",
                                "bend_IQR_deg": "bend_IQR_deg_august"}))
    return P.reset_index().rename(columns={"index": "tf"})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reparse", action="store_true")
    ap.add_argument("--cond", default="frozen", choices=["frozen", "relaxed"])
    a = ap.parse_args()
    os.makedirs(DATA, exist_ok=True)
    os.makedirs(FIGS, exist_ok=True)
    style()

    perstate = os.path.join(DATA, "bend_recon_perstate.csv")
    decomp = os.path.join(DATA, "bend_decomp_perpilot.csv")
    loc = os.path.join(DATA, "bend_localization.csv")
    if a.reparse or not all(os.path.exists(f) for f in (perstate, decomp, loc)):
        print("parsing per-state pyCurves JSON ...", flush=True)
        R, raw = parse_all()
        R.to_csv(perstate, index=False)
        decompose(raw).to_csv(decomp, index=False)
        localize(raw).to_csv(loc, index=False)
        print("parsed %d states over %d pilots" % (len(R), R.tf.nunique()), flush=True)
    R = pd.read_csv(perstate)
    D = pd.read_csv(decomp)
    L = pd.read_csv(loc)

    P = per_pilot_table(R, D, a.cond)
    rows = []
    for axis in ["bend_uu_iqr", "recon_complete_iqr", "recon_central_iqr",
                 "recon_end5_iqr", "recon_end3_iqr", "recon_orient_iqr",
                 "bend_IQR_deg_august"]:
        row = dict(axis=axis, axis_arm=a.cond)
        for eff in ("effect_relaxed", "effect_frozen"):
            S = P.dropna(subset=[axis, eff])
            r1, p1 = stats.spearmanr(S[axis], S[eff])
            row.update({f"n_{eff}": len(S), f"rho_vs_{eff}": r1, f"p_vs_{eff}": p1})
        T = P.dropna(subset=[axis, "bend_uu_iqr"])
        r2, p2 = stats.spearmanr(T[axis], T.bend_uu_iqr)
        row.update(n_bend_uu=len(T), rho_vs_bend_uu=r2, p_vs_bend_uu=p2)
        rows.append(row)
    C = pd.DataFrame(rows)
    C.to_csv(os.path.join(DATA, "bend_axis_correlations.csv"), index=False)
    P.round(5).to_csv(os.path.join(DATA, "bend_axes_perpilot.csv"), index=False)

    fig_B1(D, a.cond, FIGS)
    fig_B2(L, a.cond, FIGS)
    fig_B3(P, FIGS, a.cond)
    print("\n=== component decomposition (%s), fraction of all-varying IQR ===" % a.cond)
    print(D[D.cond == a.cond][["tf", "iqr_all", "frac_tilt", "frac_roll", "frac_twist"]]
          .round(3).to_string(index=False))
    print("\n=== axis (%s ensembles) vs effect on each arm, and vs bend_uu IQR ===" % a.cond)
    print(C.round(4).to_string(index=False))
    print("\nwrote 3 figures to %s and 5 tables to %s" % (FIGS, DATA))


if __name__ == "__main__":
    main()
