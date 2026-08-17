"""
make_F.py — F-series figures: the fnat gate behaviour.

Uses perstate_metrics.csv (stage3 rows only). Pilot-agnostic: the pilot list,
per-pilot fnat pass-rate ordering, and colors all come from fig_common, which
discovers pilots from disk. Nothing here hardcodes the pilot set.

FNAT PASS RATE per pilot = fraction of stage3 states with fnat >= 0.5.

Figures (all ordered by descending fnat pass-rate unless noted):
  F1_fnat_distributions.png : per-pilot per-state fnat strip over a light violin,
      0.5 gate floor line, points colored pass (TEAL) / fail (ALARM). n in caption.
  F2_passrate_bars.png      : bar per pilot = pass rate, frame counts labeled, TEAL bars.
  F3_fnat_vs_iRMSD.png      : scatter fnat vs iRMSD_global (all states, all pilots),
      colored by pilot; Spearman rho annotated.
  F4_interface_size.png     : n_iface_res (x) vs pass-rate (y), one point per pilot;
      Spearman rho annotated; TBP and LEF1 labeled.
"""
from fig_common import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr

apply_style()
FNAT_GATE = 0.5

# ---------------------------------------------------------------------------
# Load + pass rate
# ---------------------------------------------------------------------------
ps = pd.read_csv(os.path.join(DATA_DIR, "perstate_metrics.csv"))
ps = ps[ps.stage == "stage3"].copy()

pilots = sorted(ps.pilot.unique())
passrate = {tf: float((ps[ps.pilot == tf].fnat >= FNAT_GATE).mean()) for tf in pilots}
nstates  = {tf: int((ps.pilot == tf).sum()) for tf in pilots}
order    = order_by_passrate(pilots, passrate)
pcol     = pilot_color_ordered(pilots)


def _xlabels(ax, order):
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([short_of(t) for t in order], rotation=45, ha="right")


# ---------------------------------------------------------------------------
# F1 — per-pilot per-state fnat: strip over light violin, colored pass/fail
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8.2, 3.6))
rng = np.random.default_rng(0)
for i, tf in enumerate(order):
    v = ps[ps.pilot == tf].fnat.values
    if len(v) >= 3 and np.ptp(v) > 0:
        vp = ax.violinplot([v], positions=[i], widths=0.78, showextrema=False)
        for b in vp["bodies"]:
            b.set_facecolor(GREY_R[0]); b.set_edgecolor("none")
            b.set_alpha(0.55)
    xj = i + rng.uniform(-0.16, 0.16, size=len(v))
    passed = v >= FNAT_GATE
    ax.scatter(xj[passed],  v[passed],  s=7, c=TEAL,  alpha=0.75, linewidths=0, zorder=3)
    ax.scatter(xj[~passed], v[~passed], s=7, c=ALARM, alpha=0.75, linewidths=0, zorder=3)
ax.axhline(FNAT_GATE, color="0.35", lw=1.0, ls="--", zorder=2)
ax.text(len(order) - 0.4, FNAT_GATE + 0.015, "fnat gate = 0.5", ha="right", va="bottom",
        fontsize=6, color="0.35")
ax.set_ylabel("fnat  (fraction native interface contacts)")
ax.set_ylim(-0.02, 1.02)
_xlabels(ax, order)
ax.set_title("Per-state interface fidelity by pilot — fnat 0.5 gate")
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([0],[0],marker="o",ls="",mfc=TEAL,mec="none",ms=5,label="pass (fnat ≥ 0.5)"),
                   Line2D([0],[0],marker="o",ls="",mfc=ALARM,mec="none",ms=5,label="fail (fnat < 0.5)")],
          loc="lower left", frameon=False, ncol=1)
cap = "Pilots ordered by descending pass-rate.  n states per pilot: " + \
      ", ".join(f"{short_of(t)}={nstates[t]}" for t in order)
fig.text(0.01, -0.06, cap, fontsize=5.6, ha="left", va="top", wrap=True)
savefig(fig, "F1_fnat_distributions.png")
plt.close(fig)

# ---------------------------------------------------------------------------
# F2 — pass-rate bars, TEAL, frame counts labeled
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.4, 3.4))
vals = [passrate[t] for t in order]
bars = ax.bar(range(len(order)), vals, color=TEAL, width=0.72, edgecolor="white", linewidth=0.5)
for i, t in enumerate(order):
    ax.text(i, passrate[t] + 0.015, f"{passrate[t]*100:.0f}%\nn={nstates[t]}",
            ha="center", va="bottom", fontsize=5.4, color="0.2", linespacing=0.95)
ax.set_ylim(0, 1.14)
ax.set_ylabel("fnat pass-rate  (fraction of states with fnat ≥ 0.5)")
_xlabels(ax, order)
ax.set_title("Fraction of ensemble states passing the fnat gate, by pilot")
savefig(fig, "F2_passrate_bars.png")
plt.close(fig)

# ---------------------------------------------------------------------------
# F3 — fnat vs iRMSD_global, all states, colored by pilot; Spearman rho
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(5.6, 4.4))
for tf in order:
    g = ps[ps.pilot == tf]
    ax.scatter(g.iRMSD_global, g.fnat, s=9, color=pcol[tf], alpha=0.7,
               linewidths=0, label=short_of(tf))
ax.axhline(FNAT_GATE, color="0.4", lw=0.9, ls="--", zorder=1)
rho, p = spearmanr(ps.iRMSD_global, ps.fnat)
ax.text(0.97, 0.95, f"Spearman ρ = {rho:.2f}\n(n = {len(ps)} states)",
        transform=ax.transAxes, ha="right", va="top", fontsize=7,
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7", lw=0.5))
ax.set_xlabel("interface RMSD to crystal  (Å, global)")
ax.set_ylabel("fnat")
ax.set_ylim(-0.02, 1.02)
ax.set_title("Interface fidelity collapses with interface RMSD")
ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False,
          fontsize=6, handletextpad=0.2, labelspacing=0.25)
savefig(fig, "F3_fnat_vs_iRMSD.png")
plt.close(fig)

# ---------------------------------------------------------------------------
# F4 — interface size vs pass-rate, one point per pilot; TBP/LEF1 labeled
# ---------------------------------------------------------------------------
# crystal interface size = modal n_iface_res per pilot (constant per pilot)
niface = {tf: int(ps[ps.pilot == tf].n_iface_res.mode().iloc[0]) for tf in pilots}
fig, ax = plt.subplots(figsize=(5.6, 4.4))
xs = [niface[t] for t in order]; ys = [passrate[t] for t in order]
for tf in order:
    ax.scatter(niface[tf], passrate[tf], s=42, color=pcol[tf],
               edgecolor="white", linewidth=0.6, zorder=3)
rho, p = spearmanr(xs, ys)
sig = "n.s." if p >= 0.05 else f"p={p:.3g}"
ax.text(0.03, 0.05, f"Spearman ρ = {rho:.2f}  ({sig})",
        transform=ax.transAxes, ha="left", va="bottom", fontsize=7,
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7", lw=0.5))
for tf in ("tbp", "lef1"):
    if tf in niface:
        ax.annotate(tf.upper(), (niface[tf], passrate[tf]),
                    textcoords="offset points", xytext=(7, 4), fontsize=7,
                    color="0.15", fontweight="bold")
ax.set_xlabel("interface size  (n interface residues, crystal)")
ax.set_ylabel("fnat pass-rate")
ax.set_ylim(-0.03, 1.05)
ax.set_title("Interface size does not predict fidelity")
savefig(fig, "F4_interface_size.png")
plt.close(fig)

print("PASS-RATE TABLE (descending)")
for t in order:
    print(f"  {t:10s} {passrate[t]*100:6.1f}%   n={nstates[t]:3d}  n_iface={niface[t]}")
print("rendered: F1_fnat_distributions.png F2_passrate_bars.png F3_fnat_vs_iRMSD.png F4_interface_size.png")
