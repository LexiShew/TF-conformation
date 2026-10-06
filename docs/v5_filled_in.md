# Presentation v5 — missing data and figures, filled in

**2026-10-06.** Every **?** and **○** in v5 re-derived from the tables on endeavour, plus all eight `NEW`
figures built. Two of the re-derivations change what a slide says; the rest confirm it.

---

## Part A — The three **?** numbers

### A1 · S8 — the aggregation question answers itself

All four candidate aggregations, from the 2,221 per-state fnat rows:

| aggregation | source | median fnat | pass ≥ 0.5 | n |
|---|---|---|---|---|
| A. pooled over all states, 16 pairs | holo re-docked | 0.844 | 90.0% | 309 |
| A. pooled over all states, 16 pairs | deposited apo | 0.636 | 71.1% | 384 |
| A. pooled over all states, 16 pairs | BioEmu | 0.393 | 23.3% | 1528 |
| B. pooled over states, 13 ok-pairs | holo re-docked | 0.859 | 97.7% | 257 |
| B. pooled over states, 13 ok-pairs | deposited apo | 0.553 | 61.6% | 289 |
| B. pooled over states, 13 ok-pairs | BioEmu | 0.381 | 17.8% | 1242 |
| C. mean of per-pair medians, 16 pairs | holo re-docked | 0.795 | 91.9% | 16 |
| C. mean of per-pair medians, 16 pairs | deposited apo | 0.581 | 68.2% | 16 |
| C. mean of per-pair medians, 16 pairs | BioEmu | 0.399 | 23.5% | 16 |
| D. mean of per-pair medians, 13 ok-pairs | holo re-docked | 0.851 | 97.7% | 13 |
| D. mean of per-pair medians, 13 ok-pairs | deposited apo | 0.537 | 60.8% | 13 |
| D. mean of per-pair medians, 13 ok-pairs | BioEmu | 0.377 | 17.9% | 13 |

**Pick "pooled over states, 13 holo-control pairs."** The reason is in the table: changing the
*aggregation* barely moves anything (pooled 0.859/0.553/0.381 vs per-pair medians 0.851/0.537/0.377),
but changing the *pair set* moves apo by 0.08 and its pass rate by 9 points. The choice that matters is
whether the 3 pairs whose holo control fails are in or out — and S7 already says they are "not
interpretable, not negative," so they should be out.

Say on the slide: **"13 pairs whose holo control clears the gate; pooled over states."** Then the
aggregation is not a question anyone can ask you.

**v4's 0.839 / 0.629 / 0.396 and 100 / 92 / 11.5% reproduce under none of the four.** Drop them.

### A2 · S27 — conformation vs seed spread, verified

| quantity | value |
|---|---|
| conformation SD (within TF and seed, across states) | **0.234** |
| seed SD (within TF and state, across seeds) | **0.149** |
| ratio | **1.57×** |
| n conformations | **701** |

v4's "0.235 vs 0.16" essentially reproduces; the seed figure was 0.149, not 0.16. **The ✓ can go on
this slide now.** One methods point worth keeping in your head: pooling TFs inflates the conformation
SD to 0.319 because it absorbs between-TF differences. The within-TF number above is the honest one.

### A3 · S28 — v4 mixed two units in one sentence

| unit | baseline → augmented | lower in | Wilcoxon p |
|---|---|---|---|
| **pilot** (n = 12) | **0.247 → 0.202** | **10/12** | **0.034** |
| pilot × seed (n = 58) | 0.232 → 0.177 | 49/58 | 6.7 × 10⁻⁶ |

v4 quoted "0.235 → 0.180, 11/12, p = 0.027" — the *means* match the pilot×seed unit while the *count
and p-value* are pilot-shaped. It is two units in one row.

**Use the pilot unit.** Seeds within a pilot share data and conformations, so they are not independent
replicates of the effect; the pilot×seed p-value is pseudoreplicated. This is also consistent with the
pseudoreplication aside you already have in backup. v5's table is already correct — keep it.

---

## Part B — Two corrections that change a slide

### B1 · S28's "+0.025" is not an accuracy number ⚠

The cross-pilot delta is `d_r_state_crystal`. In `per_state_all_conditions.csv` the crystal row itself
has `r_vs_crystal = 1.0` — the column measures agreement with **the crystal's own prediction**, not with
the experimental motif. The accuracy column is `r_vs_exp`.

| cross-pilot delta (augmented − baseline) | mean | 95% CI |
|---|---|---|
| consistency with the crystal's prediction (`d_r_state_crystal`) | **+0.025** | [+0.014, +0.036] |
| accuracy against the experimental motif (`d_r_state_exp`) | **+0.014** | [+0.002, +0.025] |

So "**accuracy** on a genuinely unseen complex: +0.025" conflates the two. Both are positive and both
are small; **their intervals overlap heavily, so they are not distinguishable from each other.** N8
panel b plots both, labelled.

Also: 660 donor×target×seed rows exist, **638** are non-null. Quote 638.

### B2 · The within-pilot and cross-pilot numbers use different estimators

v5 says "within-pilot the same quantity is +0.283." Matched estimators:

| | mean | median |
|---|---|---|
| within-pilot | +0.246 | **+0.283** |
| cross-pilot | **+0.025** | +0.020 |

The +0.283 is a **median** and the +0.025 is a **mean**. Presented side by side as "the same quantity"
they are not comparable. Use mean/mean (**+0.246 vs +0.025**) or median/median (+0.283 vs +0.020).

---

## Part C — The ○ tables, now ✓

### S12 — the gate (all 16 pairs)

| protein | within | between | separation [95% CI] | ratio | verdict |
|---|---|---|---|---|---|
| THAP1 | 1.02 | 4.30 | +3.28 [+3.02, +3.45] | 4.22 | informative |
| TrpR | 1.70 | 4.41 | +2.71 [+2.02, +3.40] | 2.60 | informative |
| Dead ringer | 0.77 | 3.43 | +2.66 [+2.49, +2.77] | 4.45 | informative |
| MeCP2 MBD | 1.56 | 4.08 | +2.52 [+2.26, +2.76] | 2.62 | informative |
| Repressor C | 0.93 | 2.66 | +1.73 [+1.57, +1.85] | 2.85 | informative |
| Lac repressor | 0.60 | 2.10 | +1.49 [+1.37, +1.59] | 3.46 | informative |
| WRKY4 | 1.18 | 2.46 | +1.28 [+1.11, +1.46] | 2.08 | informative |
| FOXD3 | 1.16 | 2.26 | +1.10 [+0.93, +1.26] | 1.95 | informative |
| NHP6A | 0.51 | 1.55 | +1.05 [+0.96, +1.12] | 3.08 | informative |
| Rok | 0.75 | 1.75 | +1.00 [+0.89, +1.11] | 2.33 | informative |
| VND/NK-2 | 0.41 | 1.15 | +0.75 [+0.69, +0.81] | 2.83 | informative |
| YB-1 | 0.62 | 1.23 | +0.61 [+0.54, +0.68] | 2.00 | informative |
| c-Myb | 1.25 | 1.70 | +0.44 [+0.36, +0.52] | 1.35 | informative |
| SarA | 0.80 | 1.15 | +0.35 [+0.26, +0.45] | 1.44 | informative |
| MazE | 0.90 | 1.07 | +0.17 [-0.03, +0.30] | 1.19 | **gated out** |
| TRF1 | 1.11 | 1.28 | +0.17 [+0.09, +0.23] | 1.15 | **gated out** |

The gate is reproduced exactly by **ratio ≥ 1.25** (equivalently separation ≥ 0.25 Å) — not by the CI
excluding zero and not by the bootstrap p. Worth knowing: c-Myb (1.36) and SarA (1.44) are informative
at ratios *below* 2, so "ratio ≥ 2" would be wrong on the slide.

### S14 / S16 — position, residual, and subspace quality

| protein | family | % across | resid ratio | PC1 share | k for 90% |
|---|---|---|---|---|---|
| Rok | Rok | 2% | 3.7× | 0.70 | 7 |
| FOXD3 | Forkhead (FOXD3/Genesis) | 18% | 4.3× | 0.60 | 11 |
| SarA | MarR/winged helix (SarA) | 25% | 10.2× | 0.42 | 8 |
| NHP6A | HMG-box | 42% | 3.0× | 0.78 | 4 |
| YB-1 | Cold shock (YB-1) | 46% | 6.9× | 0.59 | 10 |
| c-Myb | MYB/SANT (c-Myb) | 47% | 1.9× | 0.30 | 17 |
| THAP1 | THAP zinc finger | 51% | 9.1× | 0.88 | 2 |
| Lac repressor | HTH (Lac repressor) | 58% | 3.1× | 0.84 | 4 |
| MeCP2 MBD | MBD (methyl-CpG) | 61% | 2.5× | 0.76 | 4 |
| VND/NK-2 | Homeodomain | 64% | 3.5× | 0.75 | 4 |
| TrpR | HTH repressor (TrpR) | 65% | 2.9× | 0.63 | 8 |
| Repressor C | HTH / phage repressor | 84% | 4.4× | 0.64 | 6 |
| Dead ringer | ARID (Dead ringer) | 85% | 5.5× | 0.88 | 2 |
| WRKY4 | WRKY zinc finger | 94% | 3.0× | 0.66 | 8 |

Median **55%**, range **2–94%**, all 14 residual ratios **> 1** (1.9–10.2×). Near-subspace subset
(ratio < 4, n = 8): median **60%**, same 2–94% range. Every v5 figure in this block is confirmed.

Gated-out pairs have PC1 shares of 0.32 (MazE) and 0.25 (TRF1) against 0.30–0.88 for the informative
set — so the "no direction dominates" claim holds for MazE and TRF1, but **c-Myb's PC1 share is 0.30**,
below MazE's. The clean version of S14's claim is the *k* column: k = 2–4 for the sharp pairs against
k = 10–15 for the gated-out ones.

### S18 — the identical-input noise control

| pair | PC1 span | sampling noise | as % of span | gap to holo |
|---|---|---|---|---|
| NHP6A | 10.86 Å | 0.173 Å (KS p = 0.64) | 1.6% | **36× noise** |
| VND/NK-2 | 6.51 Å | 0.188 Å (KS p = 0.33) | 2.9% | **13× noise** |

v5's "2–3% of the apo→holo axis" is right. **"15–25× the sampling noise" should be 13–36×.** And Rok
sits 0.23 Å from its apo bundle against ~0.18 Å of noise, so "within ~1 noise unit" is right.

### S19 — the trim plateau

| cutoff | NHP6A core / separation | VND core / separation |
|---|---|---|
| 2.0 Å | 7 pos / 0.16 Å | 4 pos / 0.12 Å |
| 2.5 Å | 50 / 1.04 | 35 / 0.74 |
| **3.0 Å** | **54 / 1.05** | **36 / 0.75** |
| 3.5 Å | 56 / 1.11 | 41 / 0.74 |
| 4.0 Å | 60 / 1.15 | 48 / 0.70 |
| 5.0 Å | 64 / 1.37 | 52 / 0.49 |
| none | 93 / **9.87** | 77 / **1.18** |

2.0 Å degeneracy (7 and 4 positions) ✓. Plateau 2.5–4.0 Å ✓. No-trim blow-up ✓ — NHP6A's untrimmed
separation is **9.87 Å** against 1.05 Å on its core. (S10's "9.47 Å / 0.81 Å" is a different quantity —
per-entry RMSD, not separation — and I did not re-derive it.)

### S21 — the mechanism split

| protein | family | apo fnat |
|---|---|---|
| c-Myb | MYB/SANT (c-Myb) | 0.765 |
| SarA | MarR/winged helix (SarA) | 0.765 |
| WRKY4 | WRKY zinc finger | 0.727 |
| Dead ringer | ARID (Dead ringer) | 0.676 |
| FOXD3 | Forkhead (FOXD3/Genesis) | 0.667 |
| YB-1 | Cold shock (YB-1) | 0.640 |
| Repressor C | HTH / phage repressor | 0.619 |
| MeCP2 MBD | MBD (methyl-CpG) | 0.577 |
| Lac repressor | HTH (Lac repressor) | 0.500 |
| Rok | Rok | 0.469 |
| TrpR | HTH repressor (TrpR) | 0.406 |
| THAP1 | THAP zinc finger | 0.393 |
| VND/NK-2 | Homeodomain | 0.349 |
| NHP6A | HMG-box | 0.194 |

Every value in v5's S21 reproduces.

---

## Part D — The eight figures

All built from the tables, no new compute. **Protein names replace the `cNN` codes throughout** (the
codes stay in N2 for cross-reference to `pairs.csv`).

| id | slide | what it shows |
|---|---|---|
| N1 | S8, S21 | three-way fnat violins, 13 holo-control pairs, 0.5 gate, n per source |
| N2 | S9 | pair roster: protein, family, PDBs, state counts, shared core, % across |
| N3 | S11 | three-step core construction + per-pair intersect → core reduction |
| N4 | S16 | % across lollipop, shaded by residual ratio |
| N5 | S12 | separation with bootstrap CIs, informative vs gated out |
| N6 | S17 | residual ratio bars with the 1× parity line |
| N7 | S20 | per-pilot BioEmu distribution, best frame, and each deposited apo source |
| N8 | S28 | paired SD slopes (pilot unit) + both cross-pilot deltas |

Name map used (from the `family` column where the deposited title is generic):

Dead ringer (ARID) · FOXD3 (Forkhead) · Lac repressor · MeCP2 MBD · NHP6A (HMG-box) · Repressor C
(HTH/phage) · Rok · SarA (MarR) · THAP1 · TrpR · VND/NK-2 (homeodomain) · WRKY4 · YB-1 (cold shock) ·
c-Myb (MYB/SANT) · MazE · TRF1 (Telobox)

Two names in v5 I could not verify from the deposited records: **"Mu repressor C"** (the entry title is
"REPRESSOR PROTEIN C", family "HTH / phage repressor" — I used "Repressor C") and **"MeCP2"** for c04
(title "Methyl-CpG Binding Protein", family "MBD (methyl-CpG)" — I kept MeCP2 MBD since 1ig4 is the
MeCP2 MBD, but it is your curation not the file's).

### Defects found while building

- **N2** — `% across` for the two gated-out pairs is computed but meaningless: MazE's is **−307%**,
  because its apo–holo separation (0.171 Å, CI crossing zero) is the denominator. The figure prints
  "— gated out" for both. Do not let that number onto a slide.
- **N5** — "informative (ratio ≥ 2)" is false; the lowest informative ratio is 1.36.
- **N8** — the two cross-pilot deltas have overlapping CIs, so "smaller still" is not supportable.

### Still to do by hand

Three items are not table-driven: the **S7** three-way design schematic, the **S30** in-flight table,
and the **S25** PyMOL render of 2LE4/9QPF superposed on 1GT0. For S25, the numbering bridge is
UniProt = auth **+0** (9QPF), **+37** (2LE4), **+38** (1GT0 chain D) — superpose on UniProt 47–103.

## Reproduce

```bash
cd /project2/rohs_102/shewchuk/TF-conformation
# S8 aggregations
python -c "import pandas as pd,glob,os; \
rows=[(os.path.basename(f).split('_fnat_')[0], f.split('_fnat_')[1][:-4], v) \
 for f in glob.glob('analysis/analyses/nmr_rmsd/fnat/*_fnat_*.csv') \
 for v in pd.read_csv(f).fnat.dropna()]; d=pd.DataFrame(rows,columns=['pair','src','fnat']); \
print(d.groupby('src').fnat.agg(['median','size']))"
# S28 metric identity
python -c "import pandas as pd; d=pd.read_csv('analysis/analyses/conf_vs_crystal_pwm/per_state_all_conditions.csv'); \
print(d[d.kind=='crystal'].r_vs_crystal.describe())"   # -> 1.0, so r_vs_crystal is not accuracy
```
