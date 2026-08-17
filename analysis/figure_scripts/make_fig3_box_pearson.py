"""fig3_box_pearson — regenerates fig3_box_pearson.png
Run from this directory after placing perentry_condition.csv and perseed_perentry.csv here.
"""
from _common import *

apply_figure_style()

arms = [('baseline', 'frozen', 'Baseline (no augmentation)', C_BASE),
        ('augmented', 'frozen', 'Augmented · frozen DNA', C_FROZ),
        ('augmented', 'relax', 'Augmented · relaxed DNA', C_RELAX)]


def entry_vals(tf, arm, dna, metric='mean_pearsonr'):
    s = pe[(pe.tf == tf) & (pe.arm == arm) & (pe.dna == dna)][metric].values
    return s


def box_figure(metric, better_txt, title, fname, xlabel, xlim=None):
    fig, ax = plt.subplots(figsize=(7.6, 6.0))
    ys = np.arange(len(order))[::-1]; off = {0: +0.26, 1: 0.0, 2: -0.26}; bw = 0.22
    for ai, (arm, dna, lab, col) in enumerate(arms):
        boxpos = []; data = []
        for yi, tf in zip(ys, order):
            v = entry_vals(tf, arm, dna, metric)
            if len(v) == 0: continue
            data.append(v); boxpos.append(yi + off[ai])
        bp = hbox(ax, data, positions=boxpos, widths=bw, patch_artist=True,
                        showfliers=False, manage_ticks=False, zorder=2)
        for p in bp['boxes']: p.set_facecolor(col); p.set_alpha(0.55); p.set_edgecolor(col); p.set_linewidth(0.8)
        for el in ['whiskers', 'caps']:
            for ln in bp[el]: ln.set_color(col); ln.set_linewidth(0.9)
        for md in bp['medians']: md.set_color('white'); md.set_linewidth(1.3)
    if len(entry_vals('csl', 'augmented', 'relax', metric))==0:
        ax.text(entry_vals('csl', 'baseline', 'frozen', metric).mean(), ys[order.index('csl')] - 0.26,
                'n.d.', color=C_RELAX, fontsize=6, va='center', ha='center', style='italic')
    ax.set_yticks(ys); ax.set_yticklabels([ylab[t] for t in order])
    ax.set_xlabel(xlabel); ax.set_title(title, loc='left'); ax.margins(y=0.02)
    if xlim: ax.set_xlim(*xlim)
    set_frame(ax)
    ax.annotate(better_txt, xy=(1.0, -0.11), xycoords='axes fraction', ha='right', va='top',
                fontsize=6, color='0.5')
    from matplotlib.lines import Line2D
    handles = [Line2D([0], [0], marker='s', color='none', markerfacecolor=c, markersize=7, label=l, alpha=0.7) for _, _, l, c in arms]
    ax.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, -0.14), ncol=3,
              frameon=False, fontsize=6.5, handletextpad=0.3, columnspacing=1.5)
    fig.tight_layout(); fig.savefig(fname, dpi=300, bbox_inches='tight')
    r = fig.canvas.get_renderer()
    tt = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text) if t.get_text().strip() and t.get_visible()]
    return sum(1 for i, (a, ba) in enumerate(tt) for b, bb in tt[i + 1:] if ba.overlaps(bb))


ov3 = box_figure('mean_pearsonr', 'higher r = better →',
    'Per-entry accuracy distribution (130 benchmark motifs): the three arms\nbarely differ against the entry-to-entry spread',
    'fig3_box_pearson.png', 'Pearson r per benchmark entry  (box = IQR, whisker = 1.5·IQR, n=130 entries)',
    xlim=(0.15, 0.98))
print("fig3 overlaps:", ov3)
