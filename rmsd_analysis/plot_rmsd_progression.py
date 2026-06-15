#!/usr/bin/env python
"""
Plot RMSD distributions from per_state_rmsds.csv (output of compute_rmsds.py).
Variant-aware: stages 0–2 are shared (variant="shared"); stage 3 has separate
"legacy" and "metal_cage" rows, so every stage-3 plot compares the two
augmentation variants.

Produces:
  - sampling_quality.png:    Stage 0 backbone RMSD vs reference, histogram per TF
                              (variant-agnostic — stage 0 is shared)
  - minimization_motion.png: Stage 3 delta_stage2_to_stage3, per TF, backbone vs
                              sidechain, two variant rows
  - bioemu_vs_min.png:        Scatter — Stage 0 RMSD vs Stage 3 min Δ (per state),
                              variants overlaid by colour
  - stagewise_progression.png: per-stage RMSD trajectory (vs reference, backbone)
                              per TF, with stage 3 split into legacy / metal_cage
  - sidechain_vs_backbone.png: scatter of sidechain vs backbone delta_stage2_to_stage3,
                              variants overlaid

Usage:
  python plot_rmsd_progression.py --csv rmsd_analysis/per_state_rmsds.csv \\
                                  --output-dir rmsd_analysis/plots
"""
import argparse
import csv
import math
from pathlib import Path


def import_matplotlib():
    global plt, np
    import matplotlib
    matplotlib.use("Agg")  # non-interactive backend
    import matplotlib.pyplot as plt_mod
    import numpy as np_mod
    plt = plt_mod
    np = np_mod


# TF display config: stable color and order
TF_ORDER = ["egr1", "tbp", "dux4"]
TF_COLORS = {"egr1": "#2c7fb8", "tbp": "#41ab5d", "dux4": "#e6550d"}
TF_LABELS = {"egr1": "EGR1 (1aay, zinc finger)",
             "tbp": "TBP (1tgh, β-saddle)",
             "dux4": "DUX4 (5z6z, homeodomain)"}

# Variant display config — same palette as the EGR1 augmentation plots.
VARIANT_ORDER = ["legacy", "metal_cage"]
VARIANT_LABELS = {
    "legacy":     "augmented (naive)",
    "metal_cage": "augmented (metal cage)",
}
VARIANT_COLORS = {
    "legacy":     "#F17AD5",
    "metal_cage": "#BD56EC",
}


FLOAT_COLS = ("rmsd_backbone_A", "rmsd_sidechain_A", "rmsd_total_A",
              "rmsd_interface_A")
INT_COLS   = ("state", "stage", "n_backbone_atoms", "n_sidechain_atoms",
              "n_total_atoms", "n_interface_atoms")


def load_rmsds(csv_path):
    """Load the per-state-RMSD CSV into a list of dicts with numeric values.
    Backward-compat: if the file lacks a 'variant' column, every row is
    treated as variant='shared'. Missing rmsd_interface_A / n_interface_atoms
    columns (older CSVs) are filled with NaN / 0."""
    rows = []
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        fields = reader.fieldnames or []
        has_variant = "variant" in fields
        for r in reader:
            for k in FLOAT_COLS:
                if k in r and r[k] not in ("", "nan", None):
                    r[k] = float(r[k])
                else:
                    r[k] = float("nan")
            for k in INT_COLS:
                if k in r and r[k] not in ("", None):
                    r[k] = int(r[k])
                else:
                    r[k] = 0
            if not has_variant:
                r["variant"] = "shared"
            rows.append(r)
    return rows


ATOM_KEY_CHOICES = ("backbone", "sidechain", "total", "interface")
ATOM_KEY_LABELS = {
    "backbone":  "Backbone",
    "sidechain": "Sidechain",
    "total":     "Total heavy-atom",
    "interface": "Interface heavy-atom",
}

# Set by main() from --atom-key. The single-metric plot functions read this
# (via _col()) rather than threading it through every signature; mixed plots
# (minimization_motion, sidechain_vs_backbone, stage_delta_violins) ignore it.
ATOM_KEY = "backbone"


def _col(atom_key=None):
    return f"rmsd_{atom_key or ATOM_KEY}_A"


def _atom_label():
    return ATOM_KEY_LABELS[ATOM_KEY]


def state_lookup(rows, *, tf, stage, variant, comparison, atom_key=None):
    """Return {state_id: value} for the filtered rows. If atom_key is None,
    defaults to the current global ATOM_KEY's column (use _col() for the
    full 'rmsd_<key>_A' name, or pass that string directly)."""
    if atom_key is None:
        atom_key = _col()
    out = {}
    for r in rows:
        if r["tf"] != tf: continue
        if r["stage"] != stage: continue
        if r["variant"] != variant: continue
        if r["comparison"] != comparison: continue
        v = r.get(atom_key, float("nan"))
        if not (isinstance(v, float) and math.isnan(v)):
            out[r["state"]] = v
    return out


def get(rows, *, tf=None, stage=None, comparison=None, variant=None,
        atom_key=None):
    """Filter rows and return a list of values for `atom_key` column (or full
    rows if atom_key=None). `variant` can be a single string or an iterable of
    accepted values."""
    if isinstance(variant, str):
        accepted = {variant}
    elif variant is None:
        accepted = None
    else:
        accepted = set(variant)

    out = []
    for r in rows:
        if tf is not None and r["tf"] != tf: continue
        if stage is not None and r["stage"] != stage: continue
        if comparison is not None and r["comparison"] != comparison: continue
        if accepted is not None and r["variant"] not in accepted: continue
        if atom_key is not None:
            v = r[atom_key]
            if not (isinstance(v, float) and math.isnan(v)):
                out.append(v)
        else:
            out.append(r)
    return out


def _annotate_stats(ax, values, *, ha="right", va="top", x=0.97, y=0.95,
                    prefix=""):
    if not values:
        return
    median = float(np.median(values))
    max_val = float(np.max(values))
    ax.text(x, y, f"{prefix}n = {len(values)}\nmedian = {median:.2f} Å\n"
                  f"max = {max_val:.2f} Å",
            transform=ax.transAxes, ha=ha, va=va, fontsize=9,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.85))


# ---------------------------------------------------------------------------
def plot_sampling_quality(rows, out_path):
    """Histogram of Stage 0 backbone RMSD vs reference, faceted per TF.
    Stage 0 is variant-shared, so this plot is unchanged from earlier."""
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharey=True)
    bins = np.linspace(0, 22, 45)
    for ax, tf in zip(axes, TF_ORDER):
        vals = get(rows, tf=tf, stage=0, comparison="vs_reference",
                   variant="shared", atom_key=_col())
        if not vals:
            ax.text(0.5, 0.5, "no data", transform=ax.transAxes, ha="center")
            ax.set_title(TF_LABELS[tf])
            continue
        ax.hist(vals, bins=bins, color=TF_COLORS[tf], edgecolor="white",
                linewidth=0.4, alpha=0.85)
        ax.axvline(float(np.median(vals)), color="black", linestyle="--",
                   linewidth=1)
        _annotate_stats(ax, vals)
        ax.set_title(TF_LABELS[tf], fontsize=11)
        ax.set_xlabel(f"{_atom_label()} RMSD vs reference crystal (Å)")
        ax.grid(axis="y", alpha=0.3)
    axes[0].set_ylabel("Number of BioEmu samples")
    fig.suptitle(f"BioEmu sampling quality: Stage 0 {_atom_label()} RMSD vs reference",
                 fontsize=12, y=1.01)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out_path}")


# ---------------------------------------------------------------------------
def plot_minimization_motion(rows, out_path):
    """Backbone vs sidechain RMSD during Stage 3 minimization. One row per
    variant, one column per TF — so each row is one augmentation strategy."""
    variants_present = [v for v in VARIANT_ORDER
                        if any(r["variant"] == v for r in rows)]
    if not variants_present:
        print(f"  no variant-tagged stage-3 rows; skipping {out_path}")
        return

    n_rows = len(variants_present)
    fig, axes = plt.subplots(n_rows, 3, figsize=(13, 4 * n_rows),
                             sharex=True, sharey=True, squeeze=False)
    bins = np.linspace(0, 2.0, 30)
    for row_i, vname in enumerate(variants_present):
        for ax, tf in zip(axes[row_i], TF_ORDER):
            bb = get(rows, tf=tf, stage=3,
                     comparison="delta_stage2_to_stage3", variant=vname,
                     atom_key="rmsd_backbone_A")
            sc = get(rows, tf=tf, stage=3,
                     comparison="delta_stage2_to_stage3", variant=vname,
                     atom_key="rmsd_sidechain_A")
            if not bb:
                ax.text(0.5, 0.5, "no data", transform=ax.transAxes,
                        ha="center")
                if row_i == 0:
                    ax.set_title(TF_LABELS[tf], fontsize=11)
                continue
            ax.hist(bb, bins=bins, color="#08519c", alpha=0.7,
                    label=f"backbone (med={np.median(bb):.2f} Å)",
                    edgecolor="white", linewidth=0.4)
            ax.hist(sc, bins=bins, color="#cb181d", alpha=0.6,
                    label=f"sidechain (med={np.median(sc):.2f} Å)",
                    edgecolor="white", linewidth=0.4)
            if row_i == 0:
                ax.set_title(TF_LABELS[tf], fontsize=11)
            if row_i == n_rows - 1:
                ax.set_xlabel("RMSD: Stage 3 vs Stage 2 (Å)")
            ax.legend(loc="upper right", fontsize=8)
            ax.grid(axis="y", alpha=0.3)
        # Row label on the left
        axes[row_i][0].set_ylabel(
            f"{VARIANT_LABELS[vname]}\n\nNumber of states", fontsize=10)
    fig.suptitle("Stage 3 minimization motion: backbone (restrained) vs "
                 "sidechain (free), by variant", fontsize=12, y=1.005)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out_path}")


# ---------------------------------------------------------------------------
def plot_bioemu_vs_min(rows, out_path):
    """Per-state scatter: x = Stage 0 backbone RMSD vs reference,
    y = Stage 3 backbone delta from Stage 2 — variants overlaid by color."""
    variants_present = [v for v in VARIANT_ORDER
                        if any(r["variant"] == v for r in rows)]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    for ax, tf in zip(axes, TF_ORDER):
        bioemu_by_state = {}
        for r in get(rows, tf=tf, stage=0, comparison="vs_reference",
                     variant="shared"):
            v = r.get(_col(), float("nan"))
            if not (isinstance(v, float) and math.isnan(v)):
                bioemu_by_state[r["state"]] = v

        plotted_any = False
        for vname in variants_present:
            min_by_state = {}
            for r in get(rows, tf=tf, stage=3,
                         comparison="delta_stage2_to_stage3", variant=vname):
                v = r.get(_col(), float("nan"))
                if not (isinstance(v, float) and math.isnan(v)):
                    min_by_state[r["state"]] = v
            common = sorted(set(bioemu_by_state) & set(min_by_state))
            if not common:
                continue
            x = [bioemu_by_state[s] for s in common]
            y = [min_by_state[s] for s in common]
            color = VARIANT_COLORS[vname]
            r_val = float(np.corrcoef(x, y)[0, 1]) if len(x) > 2 else float("nan")
            label = f"{VARIANT_LABELS[vname]}  (n={len(common)}, r={r_val:.2f})"
            ax.scatter(x, y, c=color, alpha=0.65, edgecolors="white",
                       linewidth=0.5, s=26, label=label)
            plotted_any = True

        if not plotted_any:
            ax.text(0.5, 0.5, "no data", transform=ax.transAxes, ha="center")
            ax.set_title(TF_LABELS[tf])
            continue
        ax.set_xlabel("Stage 0 vs reference (Å)\n[BioEmu sampling quality]")
        ax.set_ylabel("Stage 3 minimization Δ (Å)\n[backbone motion]")
        ax.set_title(TF_LABELS[tf], fontsize=11)
        ax.grid(alpha=0.3)
        ax.legend(loc="upper left", fontsize=7,
                  framealpha=0.85)
    fig.suptitle("Does BioEmu sampling quality predict minimization movement? "
                 "(by variant)", fontsize=12, y=1.01)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out_path}")


# ---------------------------------------------------------------------------
def plot_stagewise_progression(rows, out_path):
    """Per-TF box plots of vs_reference backbone RMSD at each stage.
    Stages 0/1/2 are shared; stage 3 is split into one box per variant.
    Two medians lines (one per variant) connect 0→1→2→3 so legacy and
    metal_cage progressions are directly comparable per TF."""
    variants_present = [v for v in VARIANT_ORDER
                        if any(r["variant"] == v and r["stage"] == 3
                               for r in rows)]
    if not variants_present:
        print(f"  no stage-3 variant rows; falling back to single-track plot")

    # X positions: stages 0,1,2 single box each; stage 3 has len(variants_present)
    # side-by-side boxes centered around x=3.
    stage_positions = {0: [0], 1: [1], 2: [2]}
    width = 0.6
    gap = 0.7
    n_v = max(len(variants_present), 1)
    s3_centers = [3 + (i - (n_v - 1) / 2) * gap for i in range(n_v)]
    stage_positions[3] = s3_centers

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), sharey=False)
    for ax, tf in zip(axes, TF_ORDER):
        all_positions = []
        all_labels = []
        # Shared stages 0,1,2
        median_by_x = {}
        for s in (0, 1, 2):
            vals = get(rows, tf=tf, stage=s, comparison="vs_reference",
                       variant="shared", atom_key=_col())
            data = vals if vals else [0]
            pos = stage_positions[s][0]
            bp = ax.boxplot([data], positions=[pos], widths=width,
                            patch_artist=True, showfliers=False,
                            medianprops=dict(color="black", linewidth=1.5))
            for patch in bp["boxes"]:
                patch.set_facecolor(TF_COLORS[tf])
                patch.set_alpha(0.55)
                patch.set_edgecolor("black")
                patch.set_linewidth(0.5)
            all_positions.append(pos)
            all_labels.append({0: "0\nBioEmu", 1: "1\nHPACKER",
                               2: "2\nDocked"}[s])
            if vals:
                median_by_x[pos] = float(np.median(vals))

        # Stage 3 per variant
        per_variant_medians = {}
        for v_i, vname in enumerate(variants_present):
            vals = get(rows, tf=tf, stage=3, comparison="vs_reference",
                       variant=vname, atom_key=_col())
            data = vals if vals else [0]
            pos = stage_positions[3][v_i]
            bp = ax.boxplot([data], positions=[pos], widths=width,
                            patch_artist=True, showfliers=False,
                            medianprops=dict(color="black", linewidth=1.5))
            for patch in bp["boxes"]:
                patch.set_facecolor(VARIANT_COLORS[vname])
                patch.set_alpha(0.85)
                patch.set_edgecolor("black")
                patch.set_linewidth(0.5)
            all_positions.append(pos)
            all_labels.append(f"3 — {VARIANT_LABELS[vname].split(' (')[1].rstrip(')')}")
            if vals:
                per_variant_medians[vname] = float(np.median(vals))

        # Median trend lines per variant: 0→1→2→3-variant
        shared_xs = [stage_positions[s][0] for s in (0, 1, 2)
                     if stage_positions[s][0] in median_by_x]
        shared_ys = [median_by_x[x] for x in shared_xs]
        for v_i, vname in enumerate(variants_present):
            if vname not in per_variant_medians:
                continue
            xs = shared_xs + [stage_positions[3][v_i]]
            ys = shared_ys + [per_variant_medians[vname]]
            ax.plot(xs, ys, color=VARIANT_COLORS[vname], linewidth=1.8,
                    marker="o", markersize=5, alpha=0.95, zorder=5,
                    label=VARIANT_LABELS[vname])

        ax.set_xticks(all_positions)
        ax.set_xticklabels(all_labels, fontsize=8)
        ax.set_xlabel("Pipeline stage")
        ax.set_ylabel(f"{_atom_label()} RMSD vs reference (Å)")
        ax.set_title(TF_LABELS[tf], fontsize=11)
        ax.grid(axis="y", alpha=0.3)
        ax.legend(loc="best", fontsize=7, framealpha=0.85)
    fig.suptitle(f"Per-stage {_atom_label()} RMSD vs reference crystal — "
                 "stage 3 split by augmentation variant",
                 fontsize=12, y=1.01)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out_path}")


# ---------------------------------------------------------------------------
def plot_sidechain_vs_backbone(rows, out_path):
    """Per-state scatter of sidechain vs backbone RMSD movement during
    Stage 3 minimization, variants overlaid by color."""
    variants_present = [v for v in VARIANT_ORDER
                        if any(r["variant"] == v for r in rows)]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    for ax, tf in zip(axes, TF_ORDER):
        plotted_any = False
        lim_max = 0.1
        for vname in variants_present:
            bb_by_state, sc_by_state = {}, {}
            for r in get(rows, tf=tf, stage=3,
                         comparison="delta_stage2_to_stage3", variant=vname):
                bv = r["rmsd_backbone_A"]
                sv = r["rmsd_sidechain_A"]
                if not (isinstance(bv, float) and math.isnan(bv)):
                    bb_by_state[r["state"]] = bv
                if not (isinstance(sv, float) and math.isnan(sv)):
                    sc_by_state[r["state"]] = sv
            common = sorted(set(bb_by_state) & set(sc_by_state))
            if not common:
                continue
            x = [bb_by_state[s] for s in common]
            y = [sc_by_state[s] for s in common]
            color = VARIANT_COLORS[vname]
            r_val = float(np.corrcoef(x, y)[0, 1]) if len(x) > 2 else float("nan")
            label = f"{VARIANT_LABELS[vname]}  (n={len(common)}, r={r_val:.2f})"
            ax.scatter(x, y, c=color, alpha=0.65, edgecolors="white",
                       linewidth=0.5, s=26, label=label)
            lim_max = max(lim_max, max(x), max(y))
            plotted_any = True
        if not plotted_any:
            ax.text(0.5, 0.5, "no data", transform=ax.transAxes, ha="center")
            ax.set_title(TF_LABELS[tf])
            continue
        ax.plot([0, lim_max * 1.1], [0, lim_max * 1.1],
                color="grey", linestyle=":", alpha=0.6, label="y = x")
        ax.set_xlabel("Backbone Δ Stage 2→3 (Å)")
        ax.set_ylabel("Sidechain Δ Stage 2→3 (Å)")
        ax.set_title(TF_LABELS[tf], fontsize=11)
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)
        ax.legend(loc="lower right", fontsize=7, framealpha=0.85)
        ax.grid(alpha=0.3)
    fig.suptitle("Stage 3 sidechain vs backbone movement, by variant",
                 fontsize=12, y=1.01)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out_path}")


# ---------------------------------------------------------------------------
# A. Per-state trajectory plot (parallel coordinates) ----------------------
def plot_state_trajectories(rows, out_path):
    """One thin line per state across stages 0,1,2 (shared) → stage 3, with
    stage 3 forking into legacy and metal_cage. Thicker medians overlaid."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), sharey=False)
    for ax, tf in zip(axes, TF_ORDER):
        # Shared per-state values at stage 0, 1, 2
        shared = {s: state_lookup(rows, tf=tf, stage=s, variant="shared",
                                   comparison="vs_reference")
                  for s in (0, 1, 2)}
        # Per-variant stage 3 values
        s3 = {v: state_lookup(rows, tf=tf, stage=3, variant=v,
                              comparison="vs_reference")
              for v in VARIANT_ORDER}

        states = sorted(set().union(*shared.values()) | set().union(*s3.values()))
        if not states:
            ax.text(0.5, 0.5, "no data", transform=ax.transAxes, ha="center")
            ax.set_title(TF_LABELS[tf]); continue

        # Faint per-state lines, one per variant fork
        x_shared = [0, 1, 2]
        for v in VARIANT_ORDER:
            color = VARIANT_COLORS[v]
            for s in states:
                ys = [shared[st].get(s) for st in (0, 1, 2)]
                y3 = s3[v].get(s)
                if y3 is None or any(y is None for y in ys):
                    continue
                ax.plot(x_shared + [3], ys + [y3], color=color, alpha=0.12,
                        linewidth=0.7)

        # Median trace per variant
        med_legend = []
        for v in VARIANT_ORDER:
            color = VARIANT_COLORS[v]
            med_y = []
            for s in (0, 1, 2):
                vals = list(shared[s].values())
                med_y.append(float(np.median(vals)) if vals else float("nan"))
            v_vals = list(s3[v].values())
            med_y.append(float(np.median(v_vals)) if v_vals else float("nan"))
            line, = ax.plot(x_shared + [3], med_y, color=color, linewidth=2.4,
                            marker="o", markersize=6, alpha=0.95, zorder=5,
                            label=VARIANT_LABELS[v])
            med_legend.append(line)

        ax.set_xticks([0, 1, 2, 3])
        ax.set_xticklabels(["0\nBioEmu", "1\nHPACKER", "2\nDocked", "3\nMinimized"],
                           fontsize=8)
        ax.set_xlabel("Pipeline stage")
        ax.set_ylabel(f"{_atom_label()} RMSD vs reference (Å)")
        ax.set_title(TF_LABELS[tf], fontsize=11)
        ax.grid(axis="y", alpha=0.3)
        ax.legend(handles=med_legend, loc="best", fontsize=8, framealpha=0.85)
    fig.suptitle(f"Per-state {_atom_label()} RMSD trajectories — every state, both variants",
                 fontsize=12, y=1.01)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out_path}")


# B. Variant agreement scatter --------------------------------------------
def plot_variant_agreement(rows, out_path):
    """Per state at stage 3, scatter metal_cage vs legacy for two quantities:
    vs_reference (top row) and delta_stage2_to_stage3 (bottom row)."""
    quantities = [
        ("vs_reference",            "Stage 3 vs reference (Å)"),
        ("delta_stage2_to_stage3",  f"Stage 2→3 {_atom_label()} Δ (Å)"),
    ]
    fig, axes = plt.subplots(len(quantities), 3, figsize=(14, 8.5))
    for row_i, (comp, ylabel) in enumerate(quantities):
        for ax, tf in zip(axes[row_i], TF_ORDER):
            legacy = state_lookup(rows, tf=tf, stage=3, variant="legacy",
                                  comparison=comp)
            metal  = state_lookup(rows, tf=tf, stage=3, variant="metal_cage",
                                  comparison=comp)
            common = sorted(set(legacy) & set(metal))
            if not common:
                ax.text(0.5, 0.5, "no data", transform=ax.transAxes, ha="center")
                if row_i == 0: ax.set_title(TF_LABELS[tf])
                continue
            x = np.array([metal[s] for s in common])
            y = np.array([legacy[s] for s in common])
            ax.scatter(x, y, c=TF_COLORS[tf], alpha=0.7, edgecolors="white",
                       linewidth=0.5, s=30)
            lim = max(x.max(), y.max()) * 1.05
            lim_min = min(x.min(), y.min()) * 0.95 if min(x.min(), y.min()) > 0 else 0
            ax.plot([lim_min, lim], [lim_min, lim], color="grey",
                    linestyle=":", alpha=0.7, label="y = x")
            r = float(np.corrcoef(x, y)[0, 1]) if len(x) > 2 else float("nan")
            ax.text(0.05, 0.95, f"n = {len(common)}\nPearson r = {r:.3f}",
                    transform=ax.transAxes, va="top",
                    bbox=dict(facecolor="white", edgecolor="none", alpha=0.85),
                    fontsize=9)
            # Annotate 3 most off-diagonal
            offs = np.abs(y - x)
            top = np.argsort(offs)[-3:]
            for idx in top:
                ax.annotate(f"{common[idx]}", (x[idx], y[idx]),
                            xytext=(4, 4), textcoords="offset points",
                            fontsize=7, color="dimgray")
            if row_i == 0:
                ax.set_title(TF_LABELS[tf], fontsize=11)
            ax.set_xlabel(f"metal_cage — {ylabel}")
            ax.set_ylabel(f"legacy — {ylabel}")
            ax.grid(alpha=0.3)
            ax.legend(loc="lower right", fontsize=8)
    fig.suptitle("Variant agreement at stage 3 — where do legacy and metal_cage diverge?",
                 fontsize=12, y=1.005)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out_path}")


# C. Signed Δ-between-variants histogram ---------------------------------
def plot_variant_delta_signed(rows, out_path):
    """Histogram of (legacy − metal_cage) per state, for vs_reference (row 1)
    and delta_stage2_to_stage3 (row 2)."""
    quantities = [
        ("vs_reference",           "Stage 3 vs reference"),
        ("delta_stage2_to_stage3", "Stage 2→3 Δ"),
    ]
    fig, axes = plt.subplots(len(quantities), 3, figsize=(14, 7.5))
    for row_i, (comp, title) in enumerate(quantities):
        for ax, tf in zip(axes[row_i], TF_ORDER):
            legacy = state_lookup(rows, tf=tf, stage=3, variant="legacy",
                                  comparison=comp)
            metal  = state_lookup(rows, tf=tf, stage=3, variant="metal_cage",
                                  comparison=comp)
            common = sorted(set(legacy) & set(metal))
            if not common:
                ax.text(0.5, 0.5, "no data", transform=ax.transAxes, ha="center")
                if row_i == 0: ax.set_title(TF_LABELS[tf])
                continue
            diffs = np.array([legacy[s] - metal[s] for s in common])
            lim = max(0.1, np.max(np.abs(diffs)) * 1.05)
            bins = np.linspace(-lim, lim, 30)
            ax.hist(diffs, bins=bins, color=TF_COLORS[tf], alpha=0.8,
                    edgecolor="white", linewidth=0.4)
            ax.axvline(0, color="black", linestyle="-", linewidth=1, alpha=0.6)
            med = float(np.median(diffs))
            ax.axvline(med, color="darkred", linestyle="--", linewidth=1.4,
                       label=f"median = {med:+.3f} Å")
            ax.text(0.97, 0.95,
                    f"n = {len(common)}\nmedian = {med:+.3f}\n"
                    f"mean = {diffs.mean():+.3f}\n"
                    f"frac > 0: {(diffs > 0).mean():.2f}",
                    transform=ax.transAxes, ha="right", va="top",
                    bbox=dict(facecolor="white", edgecolor="none", alpha=0.85),
                    fontsize=8)
            if row_i == 0:
                ax.set_title(TF_LABELS[tf], fontsize=11)
            ax.set_xlabel(f"legacy − metal_cage  ({title}, Å)")
            ax.set_ylabel("# states")
            ax.legend(loc="upper left", fontsize=8)
            ax.grid(axis="y", alpha=0.3)
    fig.suptitle("Signed Δ between variants — is one strategy systematically tighter?",
                 fontsize=12, y=1.005)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out_path}")


# D. State heatmap (sorted by stage 0) -----------------------------------
def plot_state_heatmap(rows, out_path):
    """3 rows (TFs) × 2 cols (variants). Each cell: a heatmap with rows=states
    (sorted by stage-0 RMSD), cols=stages 0,1,2,3; color=vs_reference RMSD."""
    fig, axes = plt.subplots(len(TF_ORDER), 2, figsize=(11, 14),
                             squeeze=False, sharex=False)
    for row_i, tf in enumerate(TF_ORDER):
        shared = {s: state_lookup(rows, tf=tf, stage=s, variant="shared",
                                   comparison="vs_reference")
                  for s in (0, 1, 2)}
        if not shared[0]:
            for c in (0, 1):
                axes[row_i][c].text(0.5, 0.5, "no data",
                                     transform=axes[row_i][c].transAxes,
                                     ha="center")
                axes[row_i][c].set_title(TF_LABELS[tf])
            continue
        # Sort states by stage-0 RMSD
        all_states_with_s0 = sorted(shared[0].keys(),
                                     key=lambda s: shared[0][s])
        # Compute shared color scale across both variants
        all_vals = []
        all_vals += list(shared[0].values()) + list(shared[1].values()) + list(shared[2].values())
        for v in VARIANT_ORDER:
            all_vals += list(state_lookup(rows, tf=tf, stage=3, variant=v,
                                           comparison="vs_reference").values())
        vmin, vmax = float(np.min(all_vals)), float(np.max(all_vals))

        for col_i, vname in enumerate(VARIANT_ORDER):
            ax = axes[row_i][col_i]
            s3 = state_lookup(rows, tf=tf, stage=3, variant=vname,
                              comparison="vs_reference")
            # Build matrix: rows = states (sorted), cols = stages
            states = [s for s in all_states_with_s0 if s in s3]
            mat = np.full((len(states), 4), np.nan)
            for i, s in enumerate(states):
                for c, stg in enumerate((0, 1, 2)):
                    if s in shared[stg]:
                        mat[i, c] = shared[stg][s]
                mat[i, 3] = s3[s]
            im = ax.imshow(mat, aspect="auto", cmap="viridis",
                           vmin=vmin, vmax=vmax)
            ax.set_xticks(range(4))
            ax.set_xticklabels(["0\nBioEmu", "1\nHPACKER", "2\nDocked",
                                "3\n" + VARIANT_LABELS[vname].split(" (")[1].rstrip(")")],
                               fontsize=8)
            ax.set_ylabel(f"States ({len(states)})  sorted by stage-0 RMSD" if col_i == 0
                          else "")
            ax.set_yticks([])
            if col_i == 1:
                cbar = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
                cbar.set_label("vs_reference RMSD (Å)", fontsize=8)
            if col_i == 0:
                ax.set_title(f"{TF_LABELS[tf]}  —  {VARIANT_LABELS[vname]}",
                             fontsize=10)
            else:
                ax.set_title(VARIANT_LABELS[vname], fontsize=10)
    fig.suptitle("Per-state RMSD heatmap — states sorted by BioEmu sampling quality",
                 fontsize=12, y=1.0)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out_path}")


# E. ECDF overlay per stage ----------------------------------------------
def plot_stage_ecdf(rows, out_path):
    """Empirical CDF of vs_reference backbone RMSD at each stage. Stages
    0/1/2 are shared (sequential colormap); stage 3 splits into legacy /
    metal_cage curves in the variant palette."""
    stage_colors = {0: "#cccccc", 1: "#888888", 2: "#444444"}
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), sharey=True)
    for ax, tf in zip(axes, TF_ORDER):
        for s in (0, 1, 2):
            vals = sorted(state_lookup(rows, tf=tf, stage=s, variant="shared",
                                        comparison="vs_reference").values())
            if not vals: continue
            y = np.arange(1, len(vals) + 1) / len(vals)
            ax.plot(vals, y, color=stage_colors[s], linewidth=1.8,
                    label=f"Stage {s}  (n={len(vals)})")
        for v in VARIANT_ORDER:
            vals = sorted(state_lookup(rows, tf=tf, stage=3, variant=v,
                                        comparison="vs_reference").values())
            if not vals: continue
            y = np.arange(1, len(vals) + 1) / len(vals)
            ax.plot(vals, y, color=VARIANT_COLORS[v], linewidth=2.2,
                    label=f"Stage 3 — {VARIANT_LABELS[v]}  (n={len(vals)})")
        ax.set_xlabel(f"{_atom_label()} RMSD vs reference (Å)")
        ax.set_ylabel("Empirical CDF")
        ax.set_title(TF_LABELS[tf], fontsize=11)
        ax.grid(alpha=0.3)
        ax.legend(loc="lower right", fontsize=7, framealpha=0.85)
    fig.suptitle("Per-stage ECDFs — tails visible (unlike box plots)",
                 fontsize=12, y=1.01)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out_path}")


# F. Δ-decomposition violins ---------------------------------------------
def plot_stage_delta_violins(rows, out_path):
    """Per-stage delta RMSD violins. Rows: backbone (top), sidechain (bottom).
    Cols: 3 TFs. X axis: 0→1, 1→2, 2→3 legacy, 2→3 metal_cage."""
    subset_keys = [("rmsd_backbone_A", "backbone Δ (Å)"),
                   ("rmsd_sidechain_A", "sidechain Δ (Å)")]
    transitions = [
        ("shared",     "delta_stage0_to_stage1", "0→1\nHPACKER", "#888888"),
        ("shared",     "delta_stage1_to_stage2", "1→2\nDocked",  "#888888"),
        ("legacy",     "delta_stage2_to_stage3", "2→3\nlegacy",  VARIANT_COLORS["legacy"]),
        ("metal_cage", "delta_stage2_to_stage3", "2→3\nmetal_cage", VARIANT_COLORS["metal_cage"]),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(15, 8), sharex="col")
    for row_i, (key, ylabel) in enumerate(subset_keys):
        for ax, tf in zip(axes[row_i], TF_ORDER):
            data, labels, colors = [], [], []
            for variant, comp, label, color in transitions:
                stage = int(comp.split("_to_stage")[1])
                vals = get(rows, tf=tf, stage=stage, variant=variant,
                           comparison=comp, atom_key=key)
                if vals:
                    data.append(vals)
                    labels.append(label)
                    colors.append(color)
            if not data:
                ax.text(0.5, 0.5, "no data", transform=ax.transAxes, ha="center")
                if row_i == 0: ax.set_title(TF_LABELS[tf])
                continue
            parts = ax.violinplot(data, showmeans=False, showmedians=True,
                                   widths=0.8)
            for body, color in zip(parts["bodies"], colors):
                body.set_facecolor(color)
                body.set_alpha(0.75)
                body.set_edgecolor("black")
                body.set_linewidth(0.5)
            for k in ("cmedians", "cbars", "cmins", "cmaxes"):
                if k in parts:
                    parts[k].set_color("black")
                    parts[k].set_linewidth(1.0)
            ax.set_xticks(range(1, len(labels) + 1))
            ax.set_xticklabels(labels, fontsize=8)
            ax.set_ylabel(ylabel)
            if row_i == 0:
                ax.set_title(TF_LABELS[tf], fontsize=11)
            ax.grid(axis="y", alpha=0.3)
    fig.suptitle("Per-stage motion decomposition — where do the structures move?",
                 fontsize=12, y=1.005)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out_path}")


# G. Rank stability scatter ----------------------------------------------
def _spearman(x, y):
    """Spearman ρ via NumPy: rank-correlate."""
    if len(x) < 3:
        return float("nan")
    rx = np.argsort(np.argsort(x))
    ry = np.argsort(np.argsort(y))
    return float(np.corrcoef(rx, ry)[0, 1])


def plot_rank_stability(rows, out_path):
    """For each (TF, variant): scatter state rank in stage 0 RMSD vs state
    rank in stage-3 vs_reference RMSD. Annotate Spearman ρ."""
    fig, axes = plt.subplots(len(VARIANT_ORDER), 3, figsize=(14, 8.5),
                             squeeze=False)
    for row_i, vname in enumerate(VARIANT_ORDER):
        for ax, tf in zip(axes[row_i], TF_ORDER):
            stage0 = state_lookup(rows, tf=tf, stage=0, variant="shared",
                                   comparison="vs_reference")
            s3 = state_lookup(rows, tf=tf, stage=3, variant=vname,
                              comparison="vs_reference")
            common = sorted(set(stage0) & set(s3))
            if not common:
                ax.text(0.5, 0.5, "no data", transform=ax.transAxes, ha="center")
                if row_i == 0: ax.set_title(TF_LABELS[tf])
                continue
            x = np.array([stage0[s] for s in common])
            y = np.array([s3[s] for s in common])
            rx = np.argsort(np.argsort(x))
            ry = np.argsort(np.argsort(y))
            rho = _spearman(x, y)
            ax.scatter(rx, ry, c=VARIANT_COLORS[vname], alpha=0.7,
                       edgecolors="white", linewidth=0.5, s=30)
            lim = len(common) - 1
            ax.plot([0, lim], [0, lim], color="grey", linestyle=":",
                    alpha=0.7, label="ρ = 1 (perfect rank)")
            ax.text(0.05, 0.95, f"Spearman ρ = {rho:.3f}\nn = {len(common)}",
                    transform=ax.transAxes, va="top",
                    bbox=dict(facecolor="white", edgecolor="none", alpha=0.85),
                    fontsize=9)
            if row_i == 0:
                ax.set_title(TF_LABELS[tf], fontsize=11)
            ax.set_xlabel("State rank by stage-0 RMSD\n(0 = best BioEmu sample)")
            ax.set_ylabel(f"State rank by stage-3 RMSD\n(0 = best, {vname})")
            ax.grid(alpha=0.3)
            ax.legend(loc="lower right", fontsize=8)
    fig.suptitle("Rank stability — does BioEmu quality predict final quality?",
                 fontsize=12, y=1.005)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out_path}")


# H. State survival bar --------------------------------------------------
def plot_state_survival(rows, out_path):
    """Bar chart: number of states present at each stage, per variant."""
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharey=True)
    for ax, tf in zip(axes, TF_ORDER):
        # Stages 0/1/2 are shared; stage 3 split per variant
        counts_shared = []
        for s in (0, 1, 2):
            counts_shared.append(len(state_lookup(
                rows, tf=tf, stage=s, variant="shared",
                comparison="vs_reference")))
        s3_counts = {v: len(state_lookup(
            rows, tf=tf, stage=3, variant=v, comparison="vs_reference"))
            for v in VARIANT_ORDER}

        positions, heights, colors, labels = [], [], [], []
        for s, c in zip((0, 1, 2), counts_shared):
            positions.append(s); heights.append(c); colors.append(TF_COLORS[tf])
            labels.append({0: "0\nBioEmu", 1: "1\nHPACKER",
                           2: "2\nDocked"}[s])
        gap = 0.45
        for i, v in enumerate(VARIANT_ORDER):
            offset = (i - (len(VARIANT_ORDER) - 1) / 2) * gap
            positions.append(3 + offset)
            heights.append(s3_counts[v])
            colors.append(VARIANT_COLORS[v])
            labels.append(f"3 — {VARIANT_LABELS[v].split(' (')[1].rstrip(')')}")
        bars = ax.bar(positions, heights, width=0.42, color=colors,
                      edgecolor="black", linewidth=0.5)
        for bar, h in zip(bars, heights):
            ax.text(bar.get_x() + bar.get_width()/2, h + 0.3, str(h),
                    ha="center", va="bottom", fontsize=8)
        ax.set_xticks(positions)
        ax.set_xticklabels(labels, fontsize=7)
        ax.set_ylabel("# states present")
        ax.set_title(TF_LABELS[tf], fontsize=11)
        ax.grid(axis="y", alpha=0.3)
    fig.suptitle("State survival through the pipeline — how many states make it to each stage?",
                 fontsize=12, y=1.01)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out_path}")


# ---------------------------------------------------------------------------
def main():
    global ATOM_KEY
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--csv", required=True,
                        help="Path to per_state_rmsds.csv (output of compute_rmsds.py)")
    parser.add_argument("--output-dir", required=True,
                        help="Where to write PNG plots")
    parser.add_argument("--atom-key", choices=ATOM_KEY_CHOICES,
                        default="backbone",
                        help="Which RMSD column drives the single-metric "
                             "plots (default: backbone). 'interface' requires "
                             "the CSV to include rmsd_interface_A (run "
                             "compute_rmsds.py with --interface-dir). The "
                             "backbone-vs-sidechain mixed plots ignore this.")
    args = parser.parse_args()
    ATOM_KEY = args.atom_key

    import_matplotlib()
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading {args.csv}")
    rows = load_rmsds(args.csv)
    print(f"  loaded {len(rows)} rows")
    variants = sorted({r["variant"] for r in rows})
    print(f"  variants present: {variants}")
    print(f"  atom_key:         {ATOM_KEY} (column {_col()})")

    # Sanity check for interface mode
    if ATOM_KEY == "interface":
        has_if = any(not (isinstance(r.get("rmsd_interface_A"), float)
                          and math.isnan(r["rmsd_interface_A"]))
                     for r in rows)
        if not has_if:
            raise SystemExit(
                "atom_key=interface requested but rmsd_interface_A is empty/"
                "missing in the CSV. Re-run compute_rmsds.py with "
                "--interface-dir pointing at the JSONs from "
                "identify_interface_residues.py.")

    # Non-default atom_key gets suffixed filenames so multiple runs don't overwrite.
    sfx = "" if ATOM_KEY == "backbone" else f"_{ATOM_KEY}"

    print("\nGenerating plots:")
    # Original five
    plot_sampling_quality(rows, out_dir / f"sampling_quality{sfx}.png")
    plot_minimization_motion(rows, out_dir / "minimization_motion.png")   # bb-vs-sc, ignores atom_key
    plot_bioemu_vs_min(rows, out_dir / f"bioemu_vs_min{sfx}.png")
    plot_stagewise_progression(rows, out_dir / f"stagewise_progression{sfx}.png")
    plot_sidechain_vs_backbone(rows, out_dir / "sidechain_vs_backbone.png")   # bb-vs-sc
    # Additions (A–H)
    plot_state_trajectories(rows, out_dir / f"state_trajectories{sfx}.png")
    plot_variant_agreement(rows, out_dir / f"variant_agreement{sfx}.png")
    plot_variant_delta_signed(rows, out_dir / f"variant_delta_signed{sfx}.png")
    plot_state_heatmap(rows, out_dir / f"state_heatmap{sfx}.png")
    plot_stage_ecdf(rows, out_dir / f"stage_ecdf{sfx}.png")
    plot_stage_delta_violins(rows, out_dir / "stage_delta_violins.png")   # bb-vs-sc
    plot_rank_stability(rows, out_dir / f"rank_stability{sfx}.png")
    plot_state_survival(rows, out_dir / "state_survival.png")   # counts only
    print("\nDone.")


if __name__ == "__main__":
    main()
