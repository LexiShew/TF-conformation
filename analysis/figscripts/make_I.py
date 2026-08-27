"""
make_I.py — I-series figures: interface RMSD structure.

Uses perstate_metrics.csv (stage3 rows). Pilot-agnostic; ordering by fnat
pass-rate via fig_common.

Figures:
  I1_iRMSD_distributions.png : per-pilot iRMSD_global distribution (violin + box),
      pass-rate order.
  I2_iRMSD_seg.png           : scatter iRMSD_seg_max (y) vs iRMSD_seg_mean (x),
      all states, colored by pilot; y=x line — all points sit ABOVE y=x
      (distortion is localized to one segment).
  I4_interface_geometry.png  : per-pilot crystal interface geometry — n_iface_res
      and n_segments (one constant value per pilot; modal/first per pilot).

  I3 is intentionally omitted: stage2 per-state metrics are not present in the
  current perstate_metrics.csv (only stage3 rows exist), so the stage2-vs-stage3
  comparison that I3 would show has no data. Noted here and in the run report.
"""
from fig_common import *
import matplotlib.pyplot as plt
import numpy as np

apply_style()
FNAT_GATE = 0.5

ps = pd.read_csv(os.path.join(DATA_DIR, "perstate_metrics.csv"))
ps = ps[ps.stage == "stage3"].copy()

pilots   = sorted(ps.pilot.unique())
passrate = {tf: float((ps[ps.pilot == tf].fnat >= FNAT_GATE).mean()) for tf in pilots}
nstates  = {tf: int((ps.pilot == tf).sum()) for tf in pilots}
order    = order_by_passrate(pilots, passrate)
pcol     = pilot_color_ordered(pilots)


def _xlabels(ax, order):
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([short_of(t) for t in order], rotation=45, ha="right")


# ---------------------------------------------------------------------------
# I1 — per-pilot iRMSD_global distribution: violin + box, pass-rate order
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8.2, 3.6))
data = [ps[ps.pilot == tf].iRMSD_global.values for tf in order]
vp = ax.violinplot(data, positions=range(len(order)), widths=0.8, showextrema=False)
for b in vp["bodies"]:
    b.set_facecolor(TEAL_R[0]); b.set_edgecolor("none"); b.set_alpha(0.55)
bp = ax.boxplot(data, positions=range(len(order)), widths=0.22, patch_artist=True,
                showfliers=False, medianprops=dict(color="white", lw=1.2),
                whiskerprops=dict(color=TEAL_R[2], lw=0.8),
                capprops=dict(color=TEAL_R[2], lw=0.8),
                boxprops=dict(facecolor=TEAL_R[2], edgecolor=TEAL_R[2], lw=0.5))
ax.set_ylabel("interface RMSD to crystal  (Å, global)")
_xlabels(ax, order)
ax.set_title("Interface RMSD distribution by pilot (stage3 ensemble states)")
cap = "Pilots ordered by descending fnat pass-rate.  n states: " + \
      ", ".join(f"{short_of(t)}={nstates[t]}" for t in order)
fig.text(0.01, -0.06, cap, fontsize=5.6, ha="left", va="top", wrap=True)
savefig(fig, "I1_iRMSD_distributions.png")
plt.close(fig)

# ---------------------------------------------------------------------------
# I2 — iRMSD_seg_max (y) vs iRMSD_seg_mean (x), colored by pilot; y=x line
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(5.6, 4.8))
for tf in order:
    g = ps[ps.pilot == tf]
    ax.scatter(g.iRMSD_seg_mean, g.iRMSD_seg_max, s=9, color=pcol[tf],
               alpha=0.7, linewidths=0, label=short_of(tf))
lim = max(ps.iRMSD_seg_max.max(), ps.iRMSD_seg_mean.max()) * 1.05
ax.plot([0, lim], [0, lim], color="0.4", lw=0.9, ls="--", zorder=1)
ax.text(lim * 0.62, lim * 0.55, "y = x", rotation=45, fontsize=6, color="0.4",
        ha="center", va="center")
n_above = int((ps.iRMSD_seg_max >= ps.iRMSD_seg_mean).mean() * 100)
ax.text(0.03, 0.97, f"{n_above}% of states above y=x\n(max ≥ mean by construction)",
        transform=ax.transAxes, ha="left", va="top", fontsize=6.5,
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7", lw=0.5))
ax.set_xlim(0, lim); ax.set_ylim(0, lim)
ax.set_xlabel("per-segment iRMSD, mean over segments  (Å)")
ax.set_ylabel("per-segment iRMSD, worst segment  (Å)")
ax.set_title("Distortion is localized: worst segment ≫ mean segment")
ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False,
          fontsize=6, handletextpad=0.2, labelspacing=0.25)
savefig(fig, "I2_iRMSD_seg.png")
plt.close(fig)

# ---------------------------------------------------------------------------
# I4 — crystal interface geometry: n_iface_res and n_segments per pilot
# ---------------------------------------------------------------------------
niface = {tf: int(ps[ps.pilot == tf].n_iface_res.mode().iloc[0]) for tf in pilots}
nseg   = {tf: int(ps[ps.pilot == tf].n_segments.mode().iloc[0]) for tf in pilots}
fig, ax1 = plt.subplots(figsize=(8.0, 3.6))
x = np.arange(len(order)); w = 0.38
ax1.bar(x - w/2, [niface[t] for t in order], width=w, color=TEAL,
        label="interface residues", edgecolor="white", linewidth=0.4)
ax1.set_ylabel("n interface residues", color=TEAL)
ax1.tick_params(axis="y", labelcolor=TEAL)
ax2 = ax1.twinx()
ax2.spines["top"].set_visible(False)
ax2.bar(x + w/2, [nseg[t] for t in order], width=w, color=GREY,
        label="interface segments", edgecolor="white", linewidth=0.4)
ax2.set_ylabel("n interface segments", color=GREY)
ax2.tick_params(axis="y", labelcolor=GREY)
ax2.grid(False)
ax1.set_xticks(x); ax1.set_xticklabels([short_of(t) for t in order], rotation=45, ha="right")
ax1.set_title("Crystal interface geometry by pilot (constant per pilot)")
h1, l1 = ax1.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
ax1.legend(h1 + h2, l1 + l2, loc="upper right", frameon=False, fontsize=6.5)
savefig(fig, "I4_interface_geometry.png")
plt.close(fig)

print("interface geometry (crystal, per pilot):")
for t in order:
    print(f"  {t:10s} n_iface={niface[t]:3d}  n_seg={nseg[t]:2d}")
print("I3 skipped: no stage2 rows in perstate_metrics.csv (stage3 only).")
print("rendered: I1_iRMSD_distributions.png I2_iRMSD_seg.png I4_interface_geometry.png")
