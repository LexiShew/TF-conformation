#!/usr/bin/env python3
"""
make_R.py — Cα/backbone-RMSD stage-evolution figures for the TF-conformation book.

Pilot-agnostic: pilots, ordering and colours all come from fig_common (discovered
from disk). Renders three panels into analysis/figures/:

  R1_ca_rmsd_stages.png     per-pilot whole-protein backbone RMSD-to-crystal at
                            stage2 (docked) vs stage3 (minimized), paired.
                            Minimization is local -> the two arms sit nearly on
                            top of each other.
  R2_per_residue_profiles.png  per-pilot mean per-residue Cα-RMSD profile
                            (rigid cores, floppy termini). Coverage = whatever
                            pilots exist in ca_rmsd_perresidue.csv (6 today);
                            noted in the caption.
  R3_minimization_delta.png per-state signed change in backbone RMSD-to-crystal
                            (stage3 - stage2). Negative = minimization nudged the
                            frame toward the crystal. Medians slightly negative.

DATA
  rmsd_analysis/per_state_rmsds.csv   (regenerate with compute_rmsds.py --tfs <all>)
      cols: tf,pdb_id,state,stage,variant,comparison,rmsd_backbone_A,
            rmsd_sidechain_A,rmsd_total_A,rmsd_interface_A,...
      stages 0-2 -> variant 'shared'; stage3 -> 'metal_cage' & 'legacy'.
      comparison 'vs_reference' = RMSD of that stage's frame vs the crystal.
  analysis/data/ca_rmsd_perresidue.csv  cols: pilot,resid_idx,per_res_rmsd,n_states
  analysis/data/perstate_metrics.csv    (fnat -> pass-rate ordering)

There is no pure-Cα column in per_state_rmsds.csv, so R1/R3 use backbone RMSD
(N,CA,C,O) as the whole-protein Cα proxy and label it as such. R2 uses the
genuine per-residue Cα profile from ca_rmsd_perresidue.csv.

Usage:  cd analysis/figscripts && python make_R.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "..", "common"))
from fig_common import (
    TFCONF, DATA_DIR, discover_pilots, order_by_passrate, pilot_color_ordered,
    label_of, short_of, savefig, apply_style,
    GREY, TEAL, GREEN, ALARM,
)

RMSD_CSV = os.path.join(TFCONF, "rmsd_analysis", "per_state_rmsds.csv")
PERRES_CSV = os.path.join(DATA_DIR, "ca_rmsd_perresidue.csv")
PERSTATE_CSV = os.path.join(DATA_DIR, "perstate_metrics.csv")

STAGE3_VARIANT = "metal_cage"   # canonical augmented variant for the book
RMSD_COL = "rmsd_backbone_A"    # whole-protein backbone (Cα proxy)


# ---------------------------------------------------------------------------
def compute_passrate(perstate):
    """fnat pass-rate per pilot = fraction of stage3 states with fnat >= 0.5."""
    s3 = perstate[perstate["stage"].astype(str).str.contains("3")]
    out = {}
    for tf, g in s3.groupby("pilot"):
        fn = pd.to_numeric(g["fnat"], errors="coerce").dropna()
        out[tf] = float((fn >= 0.5).mean()) if len(fn) else np.nan
    return out


def load_vs_reference(rmsd):
    """Return paired stage2/stage3 vs-reference backbone RMSD, long format.
    Columns: tf, state, stage2, stage3 (Å)."""
    vs = rmsd[rmsd["comparison"] == "vs_reference"].copy()
    s2 = vs[(vs["stage"] == 2)][["tf", "state", RMSD_COL]].rename(
        columns={RMSD_COL: "stage2"})
    s3 = vs[(vs["stage"] == 3) & (vs["variant"] == STAGE3_VARIANT)][
        ["tf", "state", RMSD_COL]].rename(columns={RMSD_COL: "stage3"})
    paired = pd.merge(s2, s3, on=["tf", "state"], how="inner")
    return paired


# ---------------------------------------------------------------------------
def fig_R1(paired, order, colors):
    """Per-pilot backbone RMSD-to-crystal, stage2 vs stage3 arms, paired."""
    fig, ax = plt.subplots(figsize=(max(7, 0.75 * len(order) + 2), 4.6))
    xpos = np.arange(len(order))
    dx = 0.16
    for i, tf in enumerate(order):
        g = paired[paired["tf"] == tf]
        if not len(g):
            continue
        for arm, off, col, mk in (("stage2", -dx, GREY, "o"),
                                   ("stage3", +dx, TEAL, "s")):
            v = g[arm].to_numpy()
            med = np.median(v)
            q1, q3 = np.percentile(v, [25, 75])
            ax.errorbar(i + off, med, yerr=[[med - q1], [q3 - med]],
                        fmt=mk, ms=6, color=col, ecolor=col, elinewidth=1.4,
                        capsize=3, zorder=3)
        # faint connector between the two medians to show they nearly coincide
        m2, m3 = np.median(g["stage2"]), np.median(g["stage3"])
        ax.plot([i - dx, i + dx], [m2, m3], color="0.6", lw=0.8, zorder=1)
    ax.set_xticks(xpos)
    ax.set_xticklabels([short_of(t) for t in order], rotation=40, ha="right")
    ax.set_ylabel("backbone RMSD to crystal (Å)")
    ax.set_xlabel("pilot (ordered by fnat pass-rate ➜)")
    ax.set_title("R1 · whole-protein backbone RMSD: docked (stage 2) vs minimized (stage 3)")
    from matplotlib.lines import Line2D
    leg = [Line2D([0], [0], marker="o", color="w", markerfacecolor=GREY,
                  ms=7, label="stage 2 (docked)"),
           Line2D([0], [0], marker="s", color="w", markerfacecolor=TEAL,
                  ms=7, label="stage 3 (minimized)")]
    ax.legend(handles=leg, frameon=False, loc="upper left")
    ax.margins(x=0.03)
    fig.tight_layout()
    return savefig(fig, "R1_ca_rmsd_stages.png")


def fig_R2(perres, colors):
    """Per-pilot mean per-residue Cα RMSD profile."""
    pilots = sorted(perres["pilot"].unique())
    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    for tf in pilots:
        g = perres[perres["pilot"] == tf].sort_values("resid_idx")
        ax.plot(g["resid_idx"], g["per_res_rmsd"], lw=1.6,
                color=colors.get(tf, "0.4"), label=short_of(tf))
    ax.set_xlabel("residue index (Cα)")
    ax.set_ylabel("mean per-residue Cα RMSD (Å)")
    ax.set_title("R2 · per-residue Cα RMSD profile (rigid cores, mobile termini)")
    ax.legend(frameon=False, ncol=2, fontsize=8, loc="upper center")
    fig.tight_layout()
    return savefig(fig, "R2_per_residue_profiles.png")


def fig_R3(paired, order, colors):
    """Signed per-state backbone-RMSD change stage3-stage2 (toward/away crystal)."""
    fig, ax = plt.subplots(figsize=(max(7, 0.75 * len(order) + 2), 4.6))
    ax.axhline(0, color="0.5", lw=1.0, zorder=1)
    data, labels, box_pos = [], [], []
    for i, tf in enumerate(order):
        g = paired[paired["tf"] == tf]
        if not len(g):
            continue
        delta = (g["stage3"] - g["stage2"]).to_numpy()
        data.append(delta)
        labels.append(short_of(tf))
        box_pos.append(i)
    bp = ax.boxplot(data, positions=box_pos, widths=0.6, showfliers=False,
                    patch_artist=True, medianprops=dict(color="black", lw=1.3))
    for patch, tf in zip(bp["boxes"], [order[i] for i in box_pos]):
        med = np.median(paired[paired["tf"] == tf]["stage3"]
                        - paired[paired["tf"] == tf]["stage2"])
        patch.set_facecolor(TEAL if med <= 0 else ALARM)
        patch.set_alpha(0.75)
    ax.set_xticks(box_pos)
    ax.set_xticklabels(labels, rotation=40, ha="right")
    ax.set_ylabel("Δ backbone RMSD to crystal,\nstage3 − stage2 (Å)")
    ax.set_xlabel("pilot (ordered by fnat pass-rate ➜)")
    ax.set_title("R3 · minimization effect: negative = frame nudged toward crystal")
    # headroom so the legend clears both the top spine and the tallest whisker
    y0, y1 = ax.get_ylim()
    ax.set_ylim(y0, y1 + 0.28 * (y1 - y0))
    from matplotlib.patches import Patch
    leg = [Patch(facecolor=TEAL, alpha=0.75, label="median ≤ 0 (toward crystal)"),
           Patch(facecolor=ALARM, alpha=0.75, label="median > 0 (away)")]
    ax.legend(handles=leg, frameon=False, loc="upper left", fontsize=8)
    fig.tight_layout()
    return savefig(fig, "R3_minimization_delta.png")


# ---------------------------------------------------------------------------
def main():
    apply_style()
    rmsd = pd.read_csv(RMSD_CSV)
    perstate = pd.read_csv(PERSTATE_CSV)

    passrate = compute_passrate(perstate)
    pilots = sorted(rmsd["tf"].unique())
    order = order_by_passrate(pilots, passrate)
    colors = pilot_color_ordered(order)

    print("pilots in per_state_rmsds:", pilots)
    print("fnat pass-rate:", {k: round(v, 3) for k, v in
                              sorted(passrate.items(), key=lambda x: -x[1])})

    paired = load_vs_reference(rmsd)
    print("paired stage2/stage3 states per pilot:",
          paired.groupby("tf").size().to_dict())

    p1 = fig_R1(paired, order, colors)
    print("wrote", p1)

    if os.path.exists(PERRES_CSV):
        perres = pd.read_csv(PERRES_CSV)
        p2 = fig_R2(perres, pilot_color_ordered(sorted(perres["pilot"].unique())))
        print("wrote", p2, "| R2 coverage pilots:",
              sorted(perres["pilot"].unique()))
    else:
        print("SKIP R2: no", PERRES_CSV)

    p3 = fig_R3(paired, order, colors)
    print("wrote", p3)


if __name__ == "__main__":
    main()
