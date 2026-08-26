"""fig2_augmentation_delta — regenerates fig2_augmentation_delta.png
Run from this directory after placing perentry_condition.csv and perseed_perentry.csv here.
"""
from _common import *

apply_figure_style()

ps = pseed  # per-seed summary (tf,arm,dna,seed,mean_*)
pilot_family = {
 'ets1':'ETS','tbp':'TBP / β-saddle','egr1':'C2H2 zinc finger','engrailed':'Homeodomain',
 'foxa':'Forkhead','lef1':'HMG-box','csl':'CSL / Rel-like','err':'Nuclear receptor',
 'nfat':'Rel / NFAT','runx':'Runt'}

ps['family'] = ps.tf.map(pilot_family)
ps['pdb'] = ps.tf.map(pilot_pdb)


arms = [('baseline','frozen','Baseline (no augmentation)',C_BASE),
        ('augmented','frozen','Augmented · frozen DNA',C_FROZ),
        ('augmented','relax','Augmented · relaxed DNA',C_RELAX)]
dnas = [('frozen',C_FROZ,'Frozen DNA',+0.16),('relax',C_RELAX,'Relaxed DNA',-0.16)]

piv2 = ps.pivot_table(index=['tf','dna','seed'], columns='arm', values='mean_pearsonr')
piv2['delta'] = piv2['augmented'] - piv2['baseline']

def seed_delta(tf, dna):
    try:
        return piv2.xs((tf, dna), level=('tf','dna'))['delta'].dropna().values
    except KeyError:
        return np.array([])

stat = {}
for dna, _, _, _ in dnas:
    for tf in order:
        d = seed_delta(tf, dna)
        if len(d) == 0:
            stat[(tf, dna)] = (np.nan, np.nan)
            continue
        m = d.mean()
        h = (sps.t.ppf(0.975, len(d)-1) * d.std(ddof=1) / np.sqrt(len(d))) if len(d) > 1 else 0
        stat[(tf, dna)] = (m, h)

ys = np.arange(len(order))[::-1]

# FIG2 final: paired ΔPearson, no divider
fig, ax = plt.subplots(figsize=(7.2, 5.2))
for dna, col, lab, dy in dnas:
    for yi, tf in zip(ys, order):
        m, h = stat[(tf, dna)]
        if np.isnan(m): continue
        lo_, hi_ = m-h, m+h
        XL = (-0.075, 0.075)
        cl_lo = max(lo_, XL[0]); cl_hi = min(hi_, XL[1])
        ax.plot([cl_lo, cl_hi], [yi+dy]*2, color=col, lw=1.8, zorder=3, solid_capstyle='butt')
        if lo_ >= XL[0]:
            ax.plot([lo_, lo_], [yi+dy-0.06, yi+dy+0.06], color=col, lw=1.2, zorder=3)
        else:
            ax.annotate('', xy=(XL[0], yi+dy), xytext=(XL[0]+0.006, yi+dy),
                        arrowprops=dict(arrowstyle='-|>', color=col, lw=1.2))
        if hi_ <= XL[1]:
            ax.plot([hi_, hi_], [yi+dy-0.06, yi+dy+0.06], color=col, lw=1.2, zorder=3)
        else:
            ax.annotate('', xy=(XL[1], yi+dy), xytext=(XL[1]-0.006, yi+dy),
                        arrowprops=dict(arrowstyle='-|>', color=col, lw=1.2))
        if XL[0] <= m <= XL[1]:
            ax.plot(m, yi+dy, 'o', color=col, ms=5.5, zorder=4, mec='white', mew=0.5)
ax.axvline(0, color='0.4', lw=1.0, zorder=1)
if len(pe[(pe.tf=='csl')&(pe.dna=='relax')])==0:
    ax.text(0.003, ys[order.index('csl')]-0.16, 'n.d.', color=C_RELAX, fontsize=6, va='center', ha='left', style='italic')
ax.set_yticks(ys); ax.set_yticklabels([ylab[t] for t in order])
ax.set_xlabel('Augmentation effect  ΔPearson = augmented − baseline   (seed-paired mean ± 95% CI; arrow = CI off-scale)')
ax.set_title('Augmentation is small and mostly negative; relaxing the DNA shifts it\nfurther negative for the rigid families (ETS1, TBP, FOXA), not positive', loc='left')
ax.set_xlim(-0.075, 0.075); ax.margins(y=0.03); set_frame(ax)
ax.annotate('← hurts        helps →', xy=(1.0, -0.12), xycoords='axes fraction', ha='right', va='top', fontsize=6, color='0.5')
handles = [Line2D([0],[0], marker='o', color=c, markersize=6, label=l, lw=0, mec='white') for _, c, l, _ in dnas]
ax.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2, frameon=False, fontsize=6.5,
          handletextpad=0.3, columnspacing=1.8, title='DNA treatment (both augmented)', title_fontsize=6.5)
fig.tight_layout()
fig.savefig('fig2_augmentation_delta.png', dpi=300, bbox_inches='tight')
