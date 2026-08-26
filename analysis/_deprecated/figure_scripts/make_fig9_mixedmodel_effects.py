"""fig9_mixedmodel_effects — regenerates fig9_mixedmodel_effects.png

Forest plot of the per-pilot augmentation effect on the pilot's OWN family,
estimated with a linear mixed-effects model that models the per-seed x per-entry
observations directly (crossed random intercepts on entry and seed). This is the
analysis that recovers power at n=5 (CSL, RUNX reach p<0.05) where the
seed-averaged Wilcoxon test is floor-locked.

Per-pilot model (same-family entries only):
    metric ~ C(arm) * C(dna)        [dna present only if the pilot has a relaxed run]
    crossed random intercepts:  ~0+C(entry_key)  and  ~0+C(seed)   (variance components)
Effects reported:
    frozen effect  = C(arm)[T.augmented]                       (arm main effect)
    relaxed effect = arm main + arm:dna interaction            (linear combination)
p-values are BH-FDR corrected within each metric across all pilot x dna tests.

Run from this directory after placing perseed_perentry.csv (+ perentry_condition.csv,
perseed_summary.csv) here.
"""
from _common import *
try:
    import statsmodels.formula.api as smf
    from statsmodels.stats.multitest import multipletests
except ModuleNotFoundError:
    import sys
    sys.exit("fig9 requires statsmodels (mixed-effects models). Install into your "
             "env, e.g.  pip install statsmodels   (or conda install -c conda-forge statsmodels). "
             "The other 8 figures do not need it.")

apply_figure_style()

DNA = [('froz', C_FROZ, 'Frozen DNA', +0.16), ('relax', GREEN[2], 'Relaxed DNA', -0.16)]
ys = list(range(len(order9)))[::-1]

# ---- fit per-pilot mixed models on same-family entries ----
def fit_pilot(tf, metric):
    """Return dict of effects/SEs/p for frozen and relaxed arm effects."""
    d = md[md.selffam & (md.tf == tf)].copy()
    if len(d) == 0:
        return None
    has_relax = d.dna.nunique() > 1
    d['grp'] = 1
    vc = {'entry': '0+C(entry_key)', 'seed': '0+C(seed)'}
    out = {}
    if has_relax:
        # baseline as reference so C(arm)[T.augmented] = augmented - baseline (the augmentation effect);
        # frozen as reference dna so the arm main effect is the frozen-DNA augmentation effect.
        f = f"{metric} ~ C(arm, Treatment('baseline'))*C(dna, Treatment('frozen'))"
        r = smf.mixedlm(f, d, groups='grp', vc_formula=vc).fit(reml=True, method='lbfgs')
        armname = [n for n in r.params.index if n.startswith('C(arm') and ':' not in n and 'augmented' in n][0]
        intname = [n for n in r.params.index if ':' in n and 'augmented' in n]
        b_arm, se_arm = r.params[armname], r.bse[armname]
        # frozen = arm main effect (reference dna = frozen, alphabetically before relax)
        out['froz'] = (b_arm, se_arm, r.pvalues[armname])
        if intname:
            iname = intname[0]
            b_int, se_int = r.params[iname], r.bse[iname]
            # relaxed = arm main + interaction; SE via covariance of the linear combo
            cov = r.cov_params()
            var = cov.loc[armname, armname] + cov.loc[iname, iname] + 2*cov.loc[armname, iname]
            se_rel = np.sqrt(var)
            eff_rel = b_arm + b_int
            from scipy.stats import norm
            p_rel = 2*(1 - norm.cdf(abs(eff_rel/se_rel)))
            out['relax'] = (eff_rel, se_rel, p_rel)
        else:
            out['relax'] = (np.nan, np.nan, np.nan)
    else:
        f = f"{metric} ~ C(arm, Treatment('baseline'))"
        r = smf.mixedlm(f, d, groups='grp', vc_formula=vc).fit(reml=True, method='lbfgs')
        armname = [n for n in r.params.index if 'augmented' in n][0]
        out['froz'] = (r.params[armname], r.bse[armname], r.pvalues[armname])
        out['relax'] = (np.nan, np.nan, np.nan)
    return out

rows = {}
row2se = {}
for metric, pc in [('m_pearsonr', 'P'), ('m_mae', 'MAE')]:
    for tf in order9:
        res = fit_pilot(tf, metric)
        if res is None:
            continue
        for dna in ['froz', 'relax']:
            eff, se, p = res[dna]
            rows.setdefault(tf, {})[f'{pc}_{dna}_effect'] = eff
            rows[tf][f'{pc}_{dna}_p'] = p
            row2se[(tf, dna, pc)] = se

mmtab = pd.DataFrame([{'tf': tf, **v} for tf, v in rows.items()])
mmtab['family'] = mmtab.tf.map(pilot_entryfam)

# BH-FDR within each metric across all (pilot x dna) tests
for pc in ['P', 'MAE']:
    pcols = [(tf, dna) for tf in mmtab.tf for dna in ['froz', 'relax']]
    pvals = [mmtab[mmtab.tf == tf].iloc[0][f'{pc}_{dna}_p'] for tf, dna in pcols]
    mask = [pd.notna(x) for x in pvals]
    padj = np.full(len(pvals), np.nan)
    if any(mask):
        padj_valid = multipletests(np.array(pvals)[mask], method='fdr_bh')[1]
        j = 0
        for i, m in enumerate(mask):
            if m:
                padj[i] = padj_valid[j]; j += 1
    bh = {pcols[i]: padj[i] for i in range(len(pcols))}
    mmtab[f'{pc}_froz_p_BH'] = mmtab.tf.map(lambda t: bh[(t, 'froz')])
    mmtab[f'{pc}_relax_p_BH'] = mmtab.tf.map(lambda t: bh[(t, 'relax')])

# ---- forest plot ----
fig, axes = plt.subplots(1, 2, figsize=(11, 6.2), sharey=True)

def plot_metric(ax, tab, pcol, title, xlab):
    for dna, col, lab, dy in DNA:
        for yi, tf in zip(ys, order9):
            row = tab[tab.tf == tf].iloc[0]
            eff = row[f'{pcol}_{dna}_effect']
            p = row[f'{pcol}_{dna}_p_BH']
            if pd.isna(eff):
                ax.text(0, yi + dy, 'n.d.', color=col, fontsize=6, va='center', ha='center', style='italic')
                continue
            se = row2se[(tf, dna, pcol)]
            sig = (p < 0.05)
            ax.errorbar(eff, yi + dy, xerr=1.96 * se, fmt='o', color=col, ms=6 if sig else 4,
                        mfc=col if sig else 'white', mec=col, capsize=2, lw=1.2, zorder=3)
    ax.axvline(0, color='0.4', lw=1, zorder=1)
    ax.set_yticks(ys)
    ax.set_yticklabels([ylab9(t) for t in order9])
    ax.set_title(title, loc='left', fontsize=9)
    ax.set_xlabel(xlab)
    set_frame(ax)

plot_metric(axes[0], mmtab, 'P', 'Accuracy (Pearson r)', 'augmentation effect on own family  (\u2192 better)')
plot_metric(axes[1], mmtab, 'MAE', 'Error (MAE)', 'augmentation effect on own family  (\u2190 better)')

mx = 0
for tf in order9:
    for dna in ['froz', 'relax']:
        for pc in ['P', 'MAE']:
            se = row2se.get((tf, dna, pc))
            e = mmtab[mmtab.tf == tf].iloc[0][f'{pc}_{dna}_effect']
            if pd.notna(se) and pd.notna(e):
                mx = max(mx, abs(e) + 1.96 * se)
lim = round(mx + 0.01, 2)
axes[0].set_xlim(-lim, lim)
axes[1].set_xlim(-lim, lim)

handles = [Line2D([0], [0], marker='o', color=C_FROZ, lw=0, label='Frozen DNA'),
           Line2D([0], [0], marker='o', color=GREEN[2], lw=0, label='Relaxed DNA'),
           Line2D([0], [0], marker='o', color='0.4', mfc='0.4', lw=0, label='filled = BH-significant (p<0.05)'),
           Line2D([0], [0], marker='o', color='0.4', mfc='white', lw=0, label='open = n.s.')]
fig.suptitle("Mixed-effects augmentation effect on the pilot's own family (crossed random effects: entry + seed)\n"
             "CSL and RUNX now significant despite n=5 — power recovered vs the seed-averaged Wilcoxon test",
             fontsize=9.5, y=1.03, x=0.02, ha='left')
fig.tight_layout(rect=[0, 0.09, 1, 1])
fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 0.06), ncol=4, frameon=False, fontsize=7.5)
fig.savefig('fig9_mixedmodel_effects.png', dpi=300, bbox_inches='tight')
