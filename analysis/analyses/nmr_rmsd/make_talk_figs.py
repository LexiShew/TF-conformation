#!/usr/bin/env python
"""Rebuild the eight talk figures (N1-N8) for the TF-conformation presentation.

Reads only tracked tables under analysis/analyses/nmr_rmsd/ and
analysis/analyses/conf_vs_crystal_pwm/. No structure files, no GPU, no network.

  python make_talk_figs.py [--out plots/talk] [--tables data/talk]

Conventions fixed here on purpose, because getting them wrong changed what a
slide said during review:

  * fnat summaries use the 13 pairs whose holo control clears the 0.5 gate
    (three_way.ok), pooled over states. Aggregation barely moves the result;
    the pair set moves apo by 0.08. See docs/v5_filled_in.md A1.
  * the cross-pilot augmentation delta is reported for BOTH r_vs_crystal
    (consistency with the crystal's own prediction; the crystal row is 1.0 by
    construction) and r_vs_exp (accuracy against the experimental motif).
    They are not the same quantity. See docs/v5_filled_in.md B1.
  * the across-conformer SD test uses the PILOT as unit (n=12). Seeds within a
    pilot share data and conformations. See docs/v5_filled_in.md A3.
  * pair identity is shown as a protein name, never the internal cNN code,
    except in N2 which carries both for cross-reference to pairs.csv.
"""
import argparse, glob, os, re
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
PWM = os.path.normpath(os.path.join(HERE, "..", "conf_vs_crystal_pwm"))

MG, HOLO, APO, BIO, ALARM = "#7a8288", "#1b6b7a", "#5f8fa6", "#b3452c", "#b3452c"

# Protein names. Taken from pair_core_summary_ca.csv 'family' where the
# deposited entry title is generic (c09 is "PROTEIN (TRANSCRIPTION FACTOR)").
# "Repressor C" is the entry title, not the curated "Mu repressor C" -- the Mu
# assignment is not in the deposited records.
SHORT = {"trf1": "TRF1", "nhp6a": "NHP6A", "vnd": "VND/NK-2", "maze": "MazE",
         "rok": "Rok", "c01": "Repressor C", "c03": "c-Myb", "c04": "MeCP2 MBD",
         "c05": "WRKY4", "c09": "FOXD3", "c10": "SarA", "c11": "Dead ringer",
         "c13": "Lac repressor", "c17": "YB-1", "c20": "THAP1", "c22": "TrpR"}


def style():
    mpl.rcParams.update({
        "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
        "legend.fontsize": 7, "xtick.labelsize": 6, "ytick.labelsize": 6,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.titlelocation": "left", "axes.titlepad": 8,
        "figure.dpi": 110, "savefig.dpi": 300})


def letter(ax, L):
    ax.text(-0.17, 1.15, L, transform=ax.transAxes, fontsize=11,
            fontweight="bold", va="top", ha="left")


def load():
    """Every table the figures need, plus the two long-format fnat tables."""
    d = {}
    for key, rel in [("pca_summary", "pair_pca_summary_ca.csv"),
                     ("bundle", "pair_bundle_summary_ca.csv"),
                     ("core", "pair_core_summary_ca.csv"),
                     ("pca", "pair_pca_ca.csv"),
                     ("trim", "trim_sweep_ca.csv"),
                     ("three_way", "fnat/fnat_three_way.csv"),
                     ("pilots_summary", "fnat_pilots/fnat_pilots_summary.csv")]:
        d[key] = pd.read_csv(os.path.join(HERE, rel))
    d["pairs"] = pd.read_csv(os.path.join(HERE, "pairs.csv"), comment="#")
    d["pwm_state"] = pd.read_csv(os.path.join(PWM, "per_state_all_conditions.csv"))
    d["pwm_tf"] = pd.read_csv(os.path.join(PWM, "per_tf_all_conditions.csv"))
    d["crosspilot"] = pd.read_csv(os.path.join(PWM, "crosspilot_deltas.csv"))

    rows = []
    for f in sorted(glob.glob(os.path.join(HERE, "fnat", "*_fnat_*.csv"))):
        b = os.path.basename(f)[:-4]
        m = re.match(r"(.+?)_fnat_(apo_states|holo_states|bioemu)$", b)
        if not m:
            continue
        t = pd.read_csv(f)
        t["pair"], t["source"] = m.group(1), m.group(2)
        rows.append(t)
    assert rows, "no per-pair fnat files matched fnat/*_fnat_*.csv"
    d["fnat_states"] = pd.concat(rows, ignore_index=True)

    # NOTE the two trees use OPPOSITE suffix order:
    #   fnat/<pair>_fnat_<source>.csv        e.g. c01_fnat_apo_states.csv
    #   fnat_pilots/<pilot>_<source>_fnat.csv  e.g. engrailed_apo_1ENH_fnat.csv
    # A single glob matches only the first and silently yields nothing for the
    # second, so each needs its own pattern.
    rows = []
    for f in sorted(glob.glob(os.path.join(HERE, "fnat_pilots", "*_fnat.csv"))):
        b = os.path.basename(f)[:-len("_fnat.csv")]
        m = re.match(r"([A-Za-z0-9]+)_(.+)$", b)
        if not m:
            continue
        t = pd.read_csv(f)
        t["pilot"] = m.group(1)
        t["entry"] = m.group(2)
        t["source"] = "bioemu" if m.group(2) == "bioemu" else "apo"
        rows.append(t)
    assert rows, "no per-pilot fnat files matched fnat_pilots/*_fnat.csv"
    d["pilot_states"] = pd.concat(rows, ignore_index=True)
    assert (d["pilot_states"].source == "bioemu").any(), "no bioemu rows in pilot_states"
    return d


def derive(d):
    """%across, residual ratio and the gate, joined into one per-pair frame."""
    piv = d["pca_summary"].pivot_table(
        index=["pair_id", "label"], columns="group",
        values=["median_pc1_A", "median_residual_A"])
    out = []
    for (pid, lb), r in piv.iterrows():
        apo, holo = r[("median_pc1_A", "apo_ref")], r[("median_pc1_A", "holo_ref")]
        ens = np.nanmean([r[("median_pc1_A", "apo_ens")], r[("median_pc1_A", "holo_ens")]])
        br = np.nanmean([r[("median_residual_A", "apo_ref")], r[("median_residual_A", "holo_ref")]])
        er = np.nanmean([r[("median_residual_A", "apo_ens")], r[("median_residual_A", "holo_ens")]])
        out.append(dict(pair_id=pid, label=lb,
                        pct_across=100 * (ens - apo) / (holo - apo),
                        resid_ratio=er / br))
    AC = (pd.DataFrame(out)
          .merge(d["pca_summary"].query("group=='apo_ref'")[
              ["pair_id", "pc1_var_frac", "k_components", "apo_holo_sep_pc1_A"]], on="pair_id")
          .merge(d["bundle"][["pair_id", "median_within_A", "median_between_A", "separation_A",
                              "sep_ci_lo_A", "sep_ci_hi_A", "ratio_between_within",
                              "informative"]], on="pair_id"))
    # the gate is reproduced exactly by ratio >= 1.25; assert so a change is caught
    assert (AC.informative == (AC.ratio_between_within >= 1.25).astype(int)).all(), \
        "informative flag no longer matches ratio >= 1.25 -- re-derive the gate"
    return AC


# --------------------------------------------------------------------------- #
# figures
# --------------------------------------------------------------------------- #
NM = {"holo_states": "holo bundle,\nDNA stripped & re-docked",
      "apo_states": "deposited\napo bundle", "bioemu": "BioEmu\nensemble"}
COL = {"holo_states": HOLO, "apo_states": APO, "bioemu": BIO}
SRCS = ["holo_states", "apo_states", "bioemu"]


def fig_N1(d, out):
    ok = set(d["three_way"].query("ok").pair)
    F = d["fnat_states"][d["fnat_states"].pair.isin(ok)]
    fig, ax = plt.subplots(figsize=(5.6, 3.8))
    for i, s in enumerate(SRCS):
        v = F[F.source == s].fnat.values
        for b in ax.violinplot([v], positions=[i], widths=.72,
                               showextrema=False)["bodies"]:
            b.set_facecolor(COL[s]); b.set_alpha(.45)
            b.set_edgecolor(COL[s]); b.set_lw(.9)
        ax.plot([i - .26, i + .26], [np.median(v)] * 2, color="black", lw=1.5, zorder=5)
        ax.text(i + .30, np.median(v), "%.3f\n%.0f%% pass" % (np.median(v), 100 * (v >= .5).mean()),
                ha="left", va="center", fontsize=6.2, zorder=6)
    ax.axhline(.5, color=MG, lw=1.0, ls="--", zorder=1)
    ax.text(-.52, .515, "gate 0.5", fontsize=6, color=MG, ha="left")
    ax.set_xticks(range(3))
    ax.set_xticklabels(["%s\n(n=%d states)" % (NM[s], (F.source == s).sum()) for s in SRCS],
                       fontsize=6.2)
    ax.set_ylabel("fnat — fraction of native\nprotein–DNA contacts recovered")
    ax.set_ylim(-.02, 1.03); ax.set_xlim(-.55, 2.75)
    ax.text(1.015, .5, "higher = better", transform=ax.transAxes, rotation=90,
            va="center", fontsize=6, color=MG)
    ax.set_title("Free-state conformations rebuild most of the interface;\n"
                 "generated ensembles do not")
    fig.savefig(os.path.join(out, "N1_three_way_fnat.png"), bbox_inches="tight")
    plt.close(fig)


def fig_N2(d, AC, out):
    C = d["core"]
    R = (C[["pair_id", "family", "apo_pdb", "holo_pdb", "n_apo_states",
            "n_holo_states", "n_core_positions"]]
         .assign(protein=lambda x: x.pair_id.map(SHORT))
         .merge(AC[["pair_id", "informative", "pct_across"]], on="pair_id")
         .sort_values(["informative", "protein"], ascending=[False, True]))
    cols = [("protein", "Protein", .00), ("family", "Family / fold", .14),
            ("apo_pdb", "apo", .40), ("holo_pdb", "holo", .47),
            ("n_apo_states", "apo\nstates", .555), ("n_holo_states", "holo\nstates", .635),
            ("n_core_positions", "shared\ncore (res)", .725), ("pct_across", "% across", .845)]
    fig, ax = plt.subplots(figsize=(7.4, 4.4)); ax.axis("off")
    for k, h, x in cols:
        ax.text(x, 1.0, h, fontsize=6.4, fontweight="bold", va="top", transform=ax.transAxes)
    for j, (_, r) in enumerate(R.iterrows()):
        y = .93 - j * .0565
        c = "black" if r.informative == 1 else ALARM
        for k, h, x in cols:
            # %across is undefined for a gated-out pair: its apo-holo separation
            # is the denominator and sits inside bundle noise (MazE gives -307%).
            s = ("%.0f%%" % r[k] if r.informative == 1 else "— gated out") if k == "pct_across" else str(r[k])
            ax.text(x, y, s, fontsize=6, va="top", color=c, transform=ax.transAxes)
    ax.text(0, .955 - len(R) * .0565 - .02,
            "14 informative pairs (black) + 2 gated out (red; their apo–holo separation is inside\n"
            "bundle noise, so % across is undefined) · 22 further pairs staged in pairs.csv",
            fontsize=6, color=MG, transform=ax.transAxes, va="top")
    ax.set_title("The pair set: 16 apo/holo solution-NMR pairs, 16 deposited family annotations")
    fig.savefig(os.path.join(out, "N2_pair_roster.png"), bbox_inches="tight")
    plt.close(fig)


def fig_N3(d, out):
    C = d["core"].assign(protein=lambda x: x.pair_id.map(SHORT))
    assert (C.n_intersect - C.n_trimmed - C.n_core_positions == 0).all(), \
        "intersect - trimmed - core != 0"
    fig = plt.figure(figsize=(7.6, 4.2))
    gs = fig.add_gridspec(1, 2, wspace=.08, width_ratios=[1, 1.35])
    a = fig.add_subplot(gs[0, 0]); a.axis("off")
    steps = [("1", "Four-way sequence intersection",
              "positions present in the apo bundle,\nthe holo bundle, and both ensembles"),
             ("2", "Minus expression tags",
              "matched by regex — Rok's two\nmembers both carry LEHHHHHH,\nso the intersection keeps it"),
             ("3", "Minus disordered termini",
              "trimmed inward while either bundle's\nintra-bundle spread exceeds 3.0 Å;\ninterior mobile loops are kept")]
    for i, (n, h, sub) in enumerate(steps):
        y = .90 - i * .265
        a.add_patch(mpl.patches.Circle((.055, y), .032, transform=a.transAxes,
                                       facecolor=HOLO, edgecolor="none", clip_on=False))
        a.text(.055, y, n, transform=a.transAxes, ha="center", va="center",
               fontsize=7, color="white", fontweight="bold")
        a.text(.13, y + .012, h, transform=a.transAxes, fontsize=7, va="center", fontweight="bold")
        a.text(.13, y - .075, sub, transform=a.transAxes, fontsize=6, va="top",
               color=MG, linespacing=1.4)
        if i < 2:
            a.annotate("", xy=(.055, y - .14), xytext=(.055, y - .055),
                       xycoords=a.transAxes, textcoords=a.transAxes,
                       arrowprops=dict(arrowstyle="-|>", lw=.9, color=MG))
    a.text(0, .115, "The core is defined by the two\ndeposited bundles only — the\n"
                    "ensembles are what is being\ntested and must not define\nthe yardstick.",
           transform=a.transAxes, fontsize=6.2, va="top", color=ALARM, linespacing=1.45)
    a.set_title("One honest atom set per pair")
    ax = fig.add_subplot(gs[0, 1])
    Z = C.sort_values("n_core_positions").reset_index(drop=True)
    y = np.arange(len(Z))
    ax.barh(y, Z.n_core_positions, color=HOLO, height=.66, zorder=3, label="shared core — scored")
    ax.barh(y, Z.n_trimmed, left=Z.n_core_positions, color=MG, alpha=.4, height=.66,
            zorder=3, label="trimmed as disordered")
    ax.set_yticks(y); ax.set_yticklabels(Z.protein, fontsize=6.4)
    for i, r in Z.iterrows():
        ax.text(r.n_intersect + 3, i, "%d/%d" % (r.n_core_positions, r.n_intersect),
                fontsize=5.8, va="center", color=MG)
    ax.set_xlabel("residue positions"); ax.set_xlim(0, 152); ax.margins(y=.03)
    ax.legend(frameon=False, fontsize=6, loc="lower right")
    lo, hi = Z.n_core_positions.min(), Z.n_core_positions.max()
    w = Z.loc[Z.n_core_positions.idxmin()]
    ax.set_title("%d–%d positions survive; %s loses %d of %d"
                 % (lo, hi, w.protein, w.n_trimmed, w.n_intersect))
    fig.savefig(os.path.join(out, "N3_shared_core.png"), bbox_inches="tight")
    plt.close(fig)


def fig_N4(AC, out):
    inf = AC.query("informative==1").sort_values("pct_across").reset_index(drop=True)
    y = np.arange(len(inf))
    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    norm = mpl.colors.Normalize(vmin=1.5, vmax=10.5)
    ax.axvline(0, color=MG, lw=1.0); ax.axvline(100, color=MG, lw=1.0)
    for i, r in inf.iterrows():
        ax.plot([0, r.pct_across], [i, i], color=MG, lw=.8, zorder=1)
    sc = ax.scatter(inf.pct_across, y, s=62, c=inf.resid_ratio,
                    cmap=mpl.colormaps["YlOrBr"], norm=norm,
                    edgecolor="black", lw=.5, zorder=3)
    ax.set_yticks(y); ax.set_yticklabels([SHORT[p] for p in inf.pair_id], fontsize=6.4)
    ax.set_xlabel("position along the apo→holo axis (%)")
    ax.set_xlim(-10, 116); ax.set_ylim(-.9, len(inf) + .4)
    ax.text(0, -.75, "apo bundle", fontsize=6, color=MG, ha="center", va="center")
    ax.text(100, -.75, "holo bundle", fontsize=6, color=MG, ha="center", va="center")
    ax.text(inf.pct_across.median(), len(inf) - .05,
            "median %.0f%%" % inf.pct_across.median(), fontsize=6.4,
            color=MG, ha="center", va="bottom")
    cb = fig.colorbar(sc, ax=ax, pad=.02, fraction=.045)
    cb.set_label("ensemble residual ÷ bundle residual", fontsize=6)
    cb.ax.tick_params(labelsize=5.5)
    ax.set_title("Where generated ensembles land spans the whole axis — %.0f%% to %.0f%%,\n"
                 "so there is no systematic apo bias (n = %d informative pairs)"
                 % (inf.pct_across.min(), inf.pct_across.max(), len(inf)))
    fig.savefig(os.path.join(out, "N4_pct_across.png"), bbox_inches="tight")
    plt.close(fig)


def fig_N5(AC, out):
    G = AC.sort_values("separation_A").reset_index(drop=True)
    y = np.arange(len(G))
    fig, ax = plt.subplots(figsize=(6.0, 4.4))
    for i, r in G.iterrows():
        c = HOLO if r.informative == 1 else ALARM
        ax.plot([r.sep_ci_lo_A, r.sep_ci_hi_A], [i, i], color=c, lw=1.6,
                solid_capstyle="round", zorder=2)
        ax.scatter([r.separation_A], [i], s=30, color=c, zorder=3,
                   marker="o" if r.informative == 1 else "D")
    ax.axvline(0, color=MG, lw=.9)
    ax.set_yticks(y); ax.set_yticklabels([SHORT[p] for p in G.pair_id], fontsize=6.4)
    for i, r in G.iterrows():
        if r.informative == 0:
            ax.get_yticklabels()[i].set_color(ALARM)
    ax.set_xlabel("apo–holo separation on the shared core (Å)\n"
                  "median between-bundle minus median within-bundle, 95% CI")
    inf, gat = G.query("informative==1"), G.query("informative==0")
    ax.legend(handles=[
        plt.Line2D([], [], marker="o", ls="", color=HOLO, ms=5,
                   label="informative — ratio ≥ 1.25 (%.1f–%.1f)"
                         % (inf.ratio_between_within.min(), inf.ratio_between_within.max())),
        plt.Line2D([], [], marker="D", ls="", color=ALARM, ms=5,
                   label="gated out — ratio %s"
                         % ", ".join("%.2f" % v for v in sorted(gat.ratio_between_within)))],
        frameon=False, loc="lower right", fontsize=6,
        title="between ÷ within bundle RMSD", title_fontsize=6)
    ax.margins(y=.03)
    ax.set_title("%d of %d pairs separate their own states by more than\n"
                 "their bundles' internal spread" % (len(inf), len(G)))
    fig.savefig(os.path.join(out, "N5_gate.png"), bbox_inches="tight")
    plt.close(fig)


def fig_N6(AC, out):
    R = AC.query("informative==1").sort_values("resid_ratio").reset_index(drop=True)
    x = np.arange(len(R))
    fig, ax = plt.subplots(figsize=(6.0, 3.7))
    ax.bar(x, R.resid_ratio, color=BIO, width=.66, zorder=3)
    ax.axhline(1, color=MG, lw=1.1, ls="--", zorder=2)
    ax.set_xticks(x); ax.set_xticklabels([SHORT[p] for p in R.pair_id],
                                         fontsize=6.2, rotation=40, ha="right")
    ax.set_ylabel("ensemble residual ÷ bundle residual")
    ax.set_ylim(0, R.resid_ratio.max() * 1.12); ax.margins(x=.02)
    ax.text(.03, .95, "every pair above 1× — the dashed line is parity with\n"
                      "the deposited bundles  (%.1f–%.1f×, n = %d)"
            % (R.resid_ratio.min(), R.resid_ratio.max(), len(R)),
            transform=ax.transAxes, fontsize=6.4, va="top", color=BIO)
    ax.set_title("What does reproduce: ensembles sit further outside the\n"
                 "observed subspace than the bundles do — in every pair")
    fig.savefig(os.path.join(out, "N6_residual_ratio.png"), bbox_inches="tight")
    plt.close(fig)


def fig_N7(d, out):
    """Per-pilot BioEmu vs deposited apo.

    The title is COMPUTED, not asserted. The claim is only defensible
    like-for-like against the NMR bundles (the solution-state format BioEmu
    emulates): it holds on 4 of 5 pilots, clearly on 3, and runx INVERTS
    because its two NMR bundles score 0.25/0.37 while its six crystal apo
    entries score 0.67-0.88. Against the crystal-dominated pooled median it
    would read 5 of 5, which is why the pooled median must not be used here.
    """
    PF, PS = d["pilot_states"], d["pilots_summary"]
    pilots = sorted(PS.pilot.unique())
    A = PS[PS.source != "bioemu"].copy()
    A["is_nmr"] = A.method.astype(str).str.contains("NMR")
    cmp_ = []
    for p in pilots:
        bm = PF[(PF.pilot == p) & (PF.source == "bioemu")].fnat.dropna()
        nmr = A[(A.pilot == p) & A.is_nmr].fnat_median
        if not len(bm) or not len(nmr):
            continue
        cmp_.append((p, nmr.median() - bm.max()))
    n_below = sum(1 for _, m in cmp_ if m > 0)
    n_clear = sum(1 for _, m in cmp_ if m > 0.10)
    inverts = [p for p, m in cmp_ if m <= 0]
    rng = np.random.default_rng(1)
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    for i, p in enumerate(pilots):
        bm = PF[(PF.pilot == p) & (PF.source == "bioemu")].fnat.dropna().values
        ax.scatter(np.full(len(bm), i) + rng.uniform(-.13, .13, len(bm)), bm,
                   s=5, color=BIO, alpha=.4, zorder=2)
        ax.plot([i - .22, i + .22], [np.median(bm)] * 2, color=BIO, lw=1.8, zorder=4)
        ax.scatter([i], [bm.max()], s=40, marker="^", facecolor="none",
                   edgecolor=BIO, lw=1.2, zorder=5)
        for _, r in PS[(PS.pilot == p) & (PS.source != "bioemu")].iterrows():
            # crystal entries are single conformers: no bundle spread, so they
            # cannot lose interface coverage to a disorder trim. Mark them apart.
            if "NMR" in str(r.method):
                ax.scatter([i + .30], [r.fnat_median], s=34, marker="s", color=HOLO, zorder=5)
            else:
                ax.scatter([i + .30], [r.fnat_median], s=26, marker="x", color=APO, lw=1.1, zorder=5)
    ax.axhline(.5, color=MG, lw=1.0, ls="--", zorder=1)
    ax.text(-.46, .515, "gate 0.5", fontsize=6, color=MG)
    ax.set_xticks(range(len(pilots))); ax.set_xticklabels(pilots, fontsize=6.4)
    ax.set_ylabel("fnat vs the pilot's own crystal interface")
    ax.set_xlim(-.5, len(pilots) - .3); ax.set_ylim(0, 1.0)
    ax.legend(handles=[
        plt.Line2D([], [], marker="o", ls="", color=BIO, ms=3, label="BioEmu frames (~90 each)"),
        plt.Line2D([], [], marker="^", ls="", mfc="none", mec=BIO, ms=6, label="BioEmu best frame"),
        plt.Line2D([], [], marker="s", ls="", color=HOLO, ms=5, label="deposited apo, NMR bundle (median)"),
        plt.Line2D([], [], marker="x", ls="", color=APO, ms=5, label="deposited apo, crystal (1 conformer)")],
        frameon=False, fontsize=5.8, loc="lower left", ncol=2, columnspacing=.8)
    for p, m in cmp_:
        if m <= 0:
            i = pilots.index(p)
            bm = PF[(PF.pilot == p) & (PF.source == "bioemu")].fnat.dropna()
            nm = sorted(A[(A.pilot == p) & A.is_nmr].fnat_median)
            ax.annotate("%s: its NMR bundles score\n%s, below the best frame"
                        % (p, " / ".join("%.2f" % v for v in nm)),
                        (i, bm.max()), textcoords="offset points", xytext=(-14, 34),
                        fontsize=5.6, color=ALARM, ha="right",
                        arrowprops=dict(arrowstyle="->", lw=.7, color=ALARM))
    ax.set_title("BioEmu's best frame out of ~90 falls short of the deposited apo\n"
                 "median on %d of %d pilots — clearly on %d%s"
                 % (n_below, len(cmp_), n_clear,
                    ", and %s inverts" % ", ".join(inverts) if inverts else ""))
    fig.savefig(os.path.join(out, "N7_pilot_fnat.png"), bbox_inches="tight")
    plt.close(fig)


def fig_N8(d, out, tabdir):
    pv = d["pwm_tf"].pivot_table(index="tf", columns="condition", values="r_sd")
    P = pv[["baseline", "aug_frozen"]].dropna()
    w = stats.wilcoxon(P.baseline, P.aug_frozen)
    CP = d["crosspilot"]; X = CP[~CP.same_pilot]
    bars = [("vs the crystal's own\nprediction\n(consistency)",
             X.d_r_state_crystal.mean(), X.d_r_state_crystal.sem(), X.d_r_state_crystal.notna().sum()),
            ("vs the experimental\nmotif\n(accuracy)",
             X.d_r_state_exp.mean(), X.d_r_state_exp.sem(), X.d_r_state_exp.notna().sum())]
    fig = plt.figure(figsize=(7.2, 3.8))
    gs = fig.add_gridspec(1, 2, wspace=.44, width_ratios=[1.25, 1])
    a1 = fig.add_subplot(gs[0, 0])
    for tf, r in P.iterrows():
        up = r.aug_frozen > r.baseline
        a1.plot([0, 1], [r.baseline, r.aug_frozen], color=ALARM if up else MG,
                lw=1.1, marker="o", ms=3.2, zorder=3 if up else 2, alpha=.95 if up else .7)
        if up:
            a1.text(1.04, r.aug_frozen, tf, fontsize=5.8, va="center", color=ALARM)
    a1.plot([0, 1], [P.baseline.mean(), P.aug_frozen.mean()], color=HOLO,
            lw=2.6, marker="o", ms=5, zorder=5)
    a1.text(1.04, P.aug_frozen.mean(), "mean", fontsize=6.2, va="center",
            color=HOLO, fontweight="bold")
    a1.set_xticks([0, 1]); a1.set_xticklabels(["baseline", "augmented\n(frozen DNA)"], fontsize=6.4)
    a1.set_xlim(-.18, 1.42)
    a1.set_ylabel("spread of PWM agreement\nacross conformations (SD)")
    a1.text(.02, .04, "lower in %d of %d pilots\nWilcoxon p = %.3f"
            % ((P.aug_frozen < P.baseline).sum(), len(P), w.pvalue),
            transform=a1.transAxes, fontsize=6.2, color=HOLO)
    a1.text(1.015, .5, "lower = more consistent", transform=a1.transAxes,
            rotation=90, va="center", fontsize=6, color=MG)
    a1.set_title("Augmentation made predictions more consistent\n"
                 "across conformations (n = %d pilots)" % len(P))
    a2 = fig.add_subplot(gs[0, 1])
    xb = np.arange(2)
    a2.bar(xb, [b[1] for b in bars], yerr=[1.96 * b[2] for b in bars],
           color=[APO, HOLO], width=.6, capsize=3, error_kw=dict(lw=.9), zorder=3)
    a2.axhline(0, color="black", lw=.9)
    for i, b in enumerate(bars):
        a2.text(i, b[1] + 1.96 * b[2] + .0015, "%+.3f" % b[1],
                ha="center", va="bottom", fontsize=6.6)
    a2.set_xticks(xb); a2.set_xticklabels([b[0] for b in bars], fontsize=6.2)
    a2.set_ylabel("change in Pearson r,\naugmented − baseline")
    a2.set_ylim(0, max(b[1] + 1.96 * b[2] for b in bars) * 1.28); a2.set_xlim(-.6, 1.75)
    a2.text(.03, .97, "cross-pilot: one pilot's conformations scored under\n"
                      "another pilot's checkpoints (%d rows, 95%% CI).\n"
                      "The two intervals overlap — they are not\ndistinguishable from each other." % bars[0][3],
            transform=a2.transAxes, fontsize=5.8, va="top", color=MG)
    a2.set_title("On genuinely unseen complexes both gains are small\n"
                 "(≤ %.3f in r)" % max(b[1] for b in bars))
    for L, ax_ in zip("ab", [a1, a2]):
        letter(ax_, L)
    fig.savefig(os.path.join(out, "N8_augmentation.png"), bbox_inches="tight")
    plt.close(fig)
    pd.DataFrame([dict(metric=b[0].replace("\n", " "), mean=round(b[1], 4),
                       ci_lo=round(b[1] - 1.96 * b[2], 4),
                       ci_hi=round(b[1] + 1.96 * b[2], 4), n=b[3]) for b in bars]
                 ).to_csv(os.path.join(tabdir, "v5_S28_crosspilot.csv"), index=False)
    P.round(4).to_csv(os.path.join(tabdir, "v5_S28_pilot_sd.csv"))


def tables(d, AC, tabdir):
    """The two derived tables that are not a figure's by-product."""
    ok = set(d["three_way"].query("ok").pair)
    F = d["fnat_states"]
    rows = []
    for name, keep, per_pair in [("A. pooled over all states, 16 pairs", set(d["three_way"].pair), False),
                                 ("B. pooled over states, 13 ok-pairs", ok, False),
                                 ("C. mean of per-pair medians, 16 pairs", set(d["three_way"].pair), True),
                                 ("D. mean of per-pair medians, 13 ok-pairs", ok, True)]:
        for s in SRCS:
            g = F[(F.source == s) & (F.pair.isin(keep))]
            if per_pair:
                med = g.groupby("pair").fnat.median().mean()
                pas = g.groupby("pair").fnat.apply(lambda v: 100 * (v >= .5).mean()).mean()
                n = g.pair.nunique()
            else:
                med, pas, n = g.fnat.median(), 100 * (g.fnat >= .5).mean(), len(g)
            rows.append(dict(aggregation=name, source=NM[s].replace("\n", " "),
                             median=round(med, 3), pass_pct=round(pas, 1), n=n))
    pd.DataFrame(rows).to_csv(os.path.join(tabdir, "v5_S8_aggregations.csv"), index=False)
    (AC.merge(d["core"][["pair_id", "family", "apo_pdb", "holo_pdb", "n_intersect",
                         "n_trimmed", "n_core_positions"]], on="pair_id")
       .assign(protein=lambda x: x.pair_id.map(SHORT)).round(4)
       .to_csv(os.path.join(tabdir, "v5_pair_table.csv"), index=False))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(HERE, "plots", "talk"))
    ap.add_argument("--tables", default=os.path.join(HERE, "data", "talk"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True); os.makedirs(a.tables, exist_ok=True)
    style()
    d = load()
    AC = derive(d)
    assert set(SHORT) == set(d["core"].pair_id), "SHORT name map is out of sync with pairs"
    fig_N1(d, a.out); fig_N2(d, AC, a.out); fig_N3(d, a.out); fig_N4(AC, a.out)
    fig_N5(AC, a.out); fig_N6(AC, a.out); fig_N7(d, a.out); fig_N8(d, a.out, a.tables)
    tables(d, AC, a.tables)
    print("wrote 8 figures to %s and 4 tables to %s" % (a.out, a.tables))


if __name__ == "__main__":
    main()
