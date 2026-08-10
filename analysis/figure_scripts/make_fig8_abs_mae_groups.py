"""fig8_abs_mae_groups — regenerates fig8_abs_mae_groups.png
Run from this directory after placing perentry_condition.csv and perseed_perentry.csv here.
"""
from _common import *

apply_figure_style()

ORANGE = ['#F3C6A2','#E8935A','#B85E22']


def six_group_vfig(metric, title, fname, ylim, ylabel):
    fig, axes = plt.subplots(2, 1, figsize=(12.5, 9.2), sharex=True)
    xs = np.arange(len(order9))
    for axi, (dna, augcols) in enumerate([('frozen', TEAL), ('relax', GREEN)]):
        ax = axes[axi]
        for xi, tf in zip(xs, order9):
            if dna == 'relax' and tf == 'csl' and len(grp6(tf,dna,'baseline','all',metric))==0:
                ax.text(xi, np.mean(ylim), 'n.d.\n(no relaxed\nDNA run)', color=GREEN[2], fontsize=6,
                        va='center', ha='center', style='italic')
                continue
            for (arm, sub), dx in zip(specs, offs):
                v = grp6(tf, dna, arm, sub, metric)
                if len(v) == 0:
                    continue
                si = subsets.index(sub)
                fc = GREY[si] if arm == 'baseline' else augcols[si]
                bp = ax.boxplot([v], positions=[xi+dx], widths=BW, patch_artist=True,
                                showfliers=False, manage_ticks=False, zorder=3 if sub=='same' else 2)
                for p in bp['boxes']:
                    p.set_facecolor(fc); p.set_edgecolor('0.35'); p.set_linewidth(0.5)
                for el in ['whiskers','caps']:
                    for ln in bp[el]: ln.set_color('0.5'); ln.set_linewidth(0.6)
                for md in bp['medians']: md.set_color('white'); md.set_linewidth(0.9)
            ax.axvline(xi, color='0.9', lw=0.5, zorder=0)
        ax.set_title(f"{'Frozen' if dna=='frozen' else 'Relaxed'} DNA", loc='left', fontsize=9)
        ax.set_ylim(*ylim); ax.set_xlim(-0.6, len(order9)-0.4); set_frame(ax)
        ax.set_ylabel(ylabel)
    axes[1].set_xticks(xs)
    axes[1].set_xticklabels([lab2[t] for t in order9], rotation=0, ha='center', fontsize=7.5)
    axes[0].tick_params(labelbottom=False)
    fig.suptitle(title, fontsize=10, y=0.99, x=0.02, ha='left')
    handles = [Patch(facecolor=GREY[0], edgecolor='0.35', label='baseline · all 130'),
               Patch(facecolor=GREY[1], edgecolor='0.35', label='baseline · other-family'),
               Patch(facecolor=GREY[2], edgecolor='0.35', label='baseline · same-family'),
               Patch(facecolor=TEAL[0], edgecolor='0.35', label='augmented · all 130   (green in Relaxed panel)'),
               Patch(facecolor=TEAL[1], edgecolor='0.35', label='augmented · other-family'),
               Patch(facecolor=TEAL[2], edgecolor='0.35', label='augmented · same-family')]
    fig.tight_layout(rect=[0, 0.08, 1, 1])
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 0.06), ncol=3, frameon=False, fontsize=7.5)
    fig.savefig(fname, dpi=300, bbox_inches='tight')
    return fig

for fg in []:
    pass

figM = six_group_vfig('mean_mae',
    "Per-entry error (MAE): baseline vs augmented, split by entry subset (all / other-family / same-family)",
    'fig8_abs_mae_groups.png', (0.15, 1.25), 'MAE per entry  (lower = better ↓)')

for ax in figM.axes[:2]:
    ax.yaxis.set_label_coords(-0.055, 0.5)
figM.subplots_adjust(left=0.10)
figM.subplots_adjust(hspace=0.22)

for ax in figM.axes[:2]:
    ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0, 1.2])

figM.savefig('fig8_abs_mae_groups.png', dpi=300, bbox_inches='tight')
print("saved fig8_abs_mae_groups.png")
