"""
make_P.py — P-series figures: family structure & augmentation accuracy.

Uses perentry_accuracy.csv (12 benchmarked pilots x 130 entries, seed-averaged)
and family_annotation.csv (fixed 130-entry benchmark, motif-level family per
entry — pilot-INDEPENDENT). Each benchmark entry maps to exactly one family at
the motif level (verified: 0 entries with >1 family), so a single entry->family
map is applied to every pilot's rows.

Figures:
  P1_family_table.png       : pilot -> family map (fig_common family_of + config PDB).
  P2_baseline_by_family.png : mean base_pearsonr by motif family (mean over pilot
      models), sorted descending.
  P3_augeffect_by_family.png: THE KEY FIGURE. Per motif family, distribution of
      per-entry ΔPearson (aug - base) pooled across pilots, median + %positive.
      Bars colored by sign (TEAL net-positive, GREY net-negative) with a 0 line.
      n per family in caption. Uses MOTIF-level family (family_annotation), not
      per-pilot family.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "..", "common"))
from fig_common import *
import matplotlib.pyplot as plt
import numpy as np

apply_style()

pe = pd.read_csv(os.path.join(DATA_DIR, "perentry_accuracy.csv"))
fa = pd.read_csv(os.path.join(DATA_DIR, "family_annotation.csv"))

# entry -> motif family (pilot-independent; one family per entry)
emap = fa.drop_duplicates("entry").set_index("entry")["family"]
assert fa.groupby("entry")["family"].nunique().max() == 1, "entry has >1 family"
pe = pe.copy()
pe["family"] = pe.entry.map(emap)
pe["dP"] = pe.aug_pearsonr - pe.base_pearsonr
assert pe.family.isna().sum() == 0, "unmapped entries"

# ---------------------------------------------------------------------------
# P1 — pilot -> family table (pilot's OWN crystal family + PDB), pass-rate order
# ---------------------------------------------------------------------------
# order pilots by fnat pass-rate (from perstate) so P1 shares the book ordering
ps = pd.read_csv(os.path.join(DATA_DIR, "perstate_metrics.csv"))
ps = ps[ps.stage == "stage3"]
allp = sorted(ps.pilot.unique())
passrate = {tf: float((ps[ps.pilot == tf].fnat >= 0.5).mean()) for tf in allp}
order = order_by_passrate(allp, passrate)

rows = []
for tf in order:
    cfg = parse_pilot_config(tf)
    rows.append([tf.upper(), family_of(tf), cfg.get("PDB_ID", "—").upper(),
                 f"{passrate[tf]*100:.0f}%"])
col_labels = ["pilot", "DBD family", "PDB", "fnat pass"]

fig, ax = plt.subplots(figsize=(6.4, 0.34 * len(rows) + 0.8))
ax.axis("off")
tbl = ax.table(cellText=rows, colLabels=col_labels, cellLoc="left", loc="center")
tbl.auto_set_font_size(False); tbl.set_fontsize(7.5); tbl.scale(1, 1.35)
for (r, c), cell in tbl.get_cells().items() if hasattr(tbl, "get_cells") else tbl._cells.items():
    cell.set_edgecolor("0.85"); cell.set_linewidth(0.5)
    if r == 0:
        cell.set_facecolor(GREY); cell.set_text_props(color="white", fontweight="bold")
    elif r % 2 == 0:
        cell.set_facecolor("#F2F4F7")
ax.set_title("Pilot transcription factors and their DBD families", pad=10, fontsize=8)
savefig(fig, "P1_family_table.png")
plt.close(fig)

# ---------------------------------------------------------------------------
# P2 — mean base_pearsonr by motif family (mean over pilot models), sorted desc
# ---------------------------------------------------------------------------
g = pe.groupby("family")
base_mean = g.base_pearsonr.mean().sort_values(ascending=False)
base_sem  = g.base_pearsonr.sem().reindex(base_mean.index)
nent      = g.entry.nunique().reindex(base_mean.index)

fig, ax = plt.subplots(figsize=(7.6, 4.0))
y = np.arange(len(base_mean))[::-1]
ax.barh(y, base_mean.values, xerr=base_sem.values, color=GREY, height=0.7,
        edgecolor="white", linewidth=0.4,
        error_kw=dict(ecolor="0.35", lw=0.8, capsize=2))
ax.set_yticks(y); ax.set_yticklabels(base_mean.index, fontsize=6.5)
for yi, fam in zip(y, base_mean.index):
    ax.text(base_mean[fam] + base_sem[fam] + 0.01, yi, f"n={nent[fam]}",
            va="center", ha="left", fontsize=5.6, color="0.35")
ax.set_xlabel("baseline PWM accuracy  (mean Pearson r, pooled over pilots)")
ax.set_xlim(0, max(0.85, (base_mean + base_sem).max() + 0.08))
ax.set_title("Baseline DeepPBS accuracy by motif family")
fig.text(0.01, -0.04,
         "Family assigned at the MOTIF level (family_annotation), pooled across all pilot models. "
         "n = benchmark entries per family. Error bars = SEM.",
         fontsize=5.6, ha="left", va="top", wrap=True)
savefig(fig, "P2_baseline_by_family.png")
plt.close(fig)

# ---------------------------------------------------------------------------
# P3 — ΔPearson by motif family, pooled across pilots. THE KEY FIGURE.
# ---------------------------------------------------------------------------
fam_stats = []
for fam, gg in pe.groupby("family"):
    fam_stats.append((fam, gg.dP.median(), (gg.dP > 0).mean(), len(gg)))
fam_stats.sort(key=lambda t: t[1], reverse=True)   # by median ΔPearson desc
fams   = [f[0] for f in fam_stats]
meds   = [f[1] for f in fam_stats]
posfr  = [f[2] for f in fam_stats]
ns     = [f[3] for f in fam_stats]

fig, ax = plt.subplots(figsize=(8.0, 4.6))
y = np.arange(len(fams))[::-1]
# violin of per-entry dP per family
data = [pe[pe.family == f].dP.values for f in fams]
vp = ax.violinplot(data, positions=y, vert=False, widths=0.85, showextrema=False)
for b, med in zip(vp["bodies"], meds):
    b.set_facecolor(TEAL_R[0] if med > 0 else GREY_R[0]); b.set_edgecolor("none")
    b.set_alpha(0.5)
# median markers, colored by sign
for yi, med in zip(y, meds):
    ax.scatter(med, yi, s=34, color=(TEAL if med > 0 else GREY), zorder=4,
               edgecolor="white", linewidth=0.6)
ax.axvline(0, color=ALARM, lw=1.0, ls="--", zorder=2)
ax.set_yticks(y); ax.set_yticklabels(fams, fontsize=6.5)
for yi, fam, pf, n in zip(y, fams, posfr, ns):
    xr = max(pe[pe.family == fam].dP.max(), 0) + 0.005
    ax.text(xr + 0.004, yi, f"{pf*100:.0f}%+  n={n}", va="center", ha="left",
            fontsize=5.4, color="0.3")
ax.set_xlabel("Δ Pearson r  (augmented − baseline), per benchmark entry")
ax.set_title("Augmentation effect on PWM accuracy by motif family")
# headroom for right-side annotation
xmax = max(d.max() for d in data); xmin = min(d.min() for d in data)
ax.set_xlim(xmin - 0.02, xmax + 0.09)
fig.text(0.01, -0.03,
         "Δ per benchmark entry (aug − base Pearson r), pooled across all pilot models by MOTIF-level "
         "family (family_annotation), NOT per-pilot family.  Dashed rose line = 0 (no effect); markers = "
         "family median; ‘%+’ = share of entries with Δ>0.  n = entries per family.",
         fontsize=5.6, ha="left", va="top", wrap=True)
savefig(fig, "P3_augeffect_by_family.png")
plt.close(fig)

print("BASELINE ACCURACY BY FAMILY (mean base Pearson r):")
for fam in base_mean.index:
    print(f"  {fam:24s} {base_mean[fam]:.3f}  n={nent[fam]}")
print("\nAUG EFFECT BY FAMILY (median ΔPearson, %positive, n):")
for fam, med, pf, n in fam_stats:
    print(f"  {fam:24s} med={med:+.4f}  {pf*100:5.1f}%+  n={n}")
print("\nrendered: P1_family_table.png P2_baseline_by_family.png P3_augeffect_by_family.png")
