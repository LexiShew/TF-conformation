# Chapter 19½ — Is BioEmu actually the apo state? (the assumption, tested)

> **What this chapter is.** Chapter 2 introduces BioEmu's ensemble as a **computational apo (free) state**
> and the crystal as **holo**, and calls that contrast "the conceptual lever the rest of the book pulls on."
> It is stated as framing and never measured. This chapter measures it.
>
> **Why it can be measured at all.** The pipeline's pilots have only a bound crystal — there is no
> experimental free state to compare against, which is why the assumption went untested. The
> `structures/nmr/` set breaks that: it holds **apo/holo NMR pairs**, two deposited bundles of the *same*
> protein, one free and one DNA-bound. That is the only place in the project where the apo/holo axis is
> defined by experiment rather than by assumption.
>
> **Sourcing.** All values are marked **[nmr]** and come from `analysis/analyses/nmr_rmsd/`:
> `pair_bundle_summary_*.csv` (the informativeness gate), `pair_pca_summary_*.csv` (position on the
> apo→holo axis), `fnat/fnat_three_way.csv` (the interface-fidelity three-way), and
> `fnat_pilots/pilot_mechanism.csv` (the pilots). Scripts are co-located; `README.md` there is the
> method of record.
>
> **The headline, up front.** Two findings, one negative for the framing and one positive and new:
> (1) BioEmu is **not** reliably a free-state ensemble — where it lands between apo and holo varies by
> system across the entire range; and (2) docking a *deposited apo structure* onto the bound DNA and
> scoring `fnat` gives a **mechanism read-out that needs no BioEmu at all**, which is the independent
> measurement Chapter 17 currently lacks.

---

## 19½.1 — The pair set and the informativeness gate

Sixteen apo/holo pairs have BioEmu ensembles for **both** members: five original (trf1, nhp6a, vnd, maze,
rok) and eleven found by PDB search (solution-NMR protein+DNA entries matched to protein-only entries in
the same 95% sequence cluster). `pairs.csv` carries all of them plus 21 further candidates at
`status: needs_ensembles`.

Every RMSD in a pair is computed on **one shared atom set** — the 4-way sequence intersection of the two
bundles and the two ensembles, minus expression tags, minus disordered termini trimmed inward while the
intra-bundle spread exceeds 3.0 Å. This is not cosmetic. Untrimmed, the numbers measure floppy tails:
`5zuz` reads **9.47 Å** median RMSD-to-closest-state untrimmed and **0.81 Å** on its core. **[nmr]**

A pair earns its place only if the apo/holo separation clears the bundles' own heterogeneity — NMR spread
is precision as much as dynamics. Gate: bootstrap CI on `median(between) − median(within)` excluding zero,
and a ratio ≥ 1.25. **14 of 16 pass**; trf1 (ratio 1.15) and maze (1.19) do not, and no conclusion should
rest on them. **[nmr]**

Two independent checks agree with that gate: a Cartesian PCA fitted on the deposited states alone puts
70–88% of the variance on PC1 for the informative pairs and only 25–32% for the two that fail; and the
`t` projection places their apo and holo states at 0 and 1 with non-overlapping ranges. **[nmr]**

## 19½.2 — BioEmu is not reliably the apo state

Fit a PCA on the **deposited states only** (apo + holo NMR), project the ensembles in as test data, and
read where they land on PC1 as `%across` — 0% at the apo bundle, 100% at holo.

| | %across | | | %across |
|---|---|---|---|---|
| rok | 2% | | c04 MeCP2 | 61% |
| c09 FOXD3 | 18% | | vnd | 64% |
| c10 SarA | 25% | | c22 TrpR | 65% |
| nhp6a | 42% | | c01 Mu repressor | 84% |
| c17 YB-1 | 46% | | c11 Dead ringer | 85% |
| c03 c-Myb | 47% | | c05 WRKY4 | 94% |
| c20 THAP1 | 51% | | | |
| c13 Lac repressor | 58% | | | |

**[nmr]**

`%across` spans **2% to 94%, median 55%**, and restricting to the pairs whose ensembles stay nearest the
deposited subspace does not change it (median 60%, same range). Two pairs sit at the apo end, three at the
holo end, nine in between. **BioEmu is not systematically apo-biased; where its ensemble lands varies by
system across the whole range.** c05 (WRKY4) at 94% is essentially *at* the bound state.

The noise floor makes that a real spread rather than scatter. For nhp6a and vnd the two construct sequences
are **byte-identical**, so their two ensembles are independent draws from one distribution: PC1 medians
differ by 0.17 Å and 0.19 Å and KS cannot separate them (p = 0.64, 0.33) — about **2–3% of the apo→holo
axis**. The displacements above are 15–25× that. **[nmr]**

What *is* consistent: the ensembles are wider than the deposited bundles and sit outside their subspace, by
**2–10× the bundles' own leave-one-out residual**. That displacement, not a directional apo bias, is the
reproducible BioEmu signature. **[nmr]**

> **A superseded reading the successor should not resurrect.** On the first three pairs alone (nhp6a, rok,
> vnd) this looked like "BioEmu reaches neither basin and never gets to holo." All three happened to land
> apo-of-centre. At n = 14 it does not survive. A related scalar — projecting onto the apo-medoid→holo-medoid
> *line* and reading the leftover distance `p` — also overstated the off-axis story, because the deposited
> bundles themselves score p/L 0.15–0.40 on it. The PCA residual against the whole deposited subspace is the
> fair version.

## 19½.3 — Interface fidelity: BioEmu vs the real free state

The holo member of every pair carries its own DNA, so the pipeline's Stage-2 dock and `fnat_gate` can be
pointed at structures where the apo/holo answer is known independently. Three sources through the *same*
unmodified dock, all trimmed to one residue set:

| source | median fnat | median pass @ `FNAT_FLOOR` 0.5 |
|---|---|---|
| deposited **holo** states, DNA stripped, re-docked | 0.84 | **100%** |
| deposited **apo** states | 0.58 | **81%** |
| **BioEmu** ensemble | 0.38 | **6%** |

**[nmr]**

Each row answers a different objection. Holo-re-docked shows the **dock is lossless and the 0.5 floor is
reachable** — a structure that *is* the bound state survives the rigid-body round trip. Real apo mostly
**clears the gate**, so the gate is not rejecting free-state conformations as such. What remains is BioEmu:
below the genuine free state of the same protein in **10 of 13 pairs (Wilcoxon p = 0.002)**.

One fairness correction. BioEmu contributes ~80–100 broadly sampled frames against 20 models of one refined
structure. Comparing BioEmu's **best k frames** against each pair's k apo models narrows the gap to 0.46 vs
0.58, no longer significant (p = 0.15), and BioEmu's single best frame reaches 0.5 in **10 of 13** pairs.
**[nmr]** So BioEmu does produce interface-competent conformations; they are a small minority of what it
generates. That is a sampling-breadth statement, not an accuracy failure, and it is the defensible framing.

**`%across` does not predict fnat** (Spearman +0.36, p = 0.23). The tempting bridge — "fnat pass rate as a
cheap pipeline-side proxy for where the ensemble sits on the apo/holo axis" — is not supported. **[nmr]**

> **Method notes that bit, and will bite again.** (i) Scoring deposited holo states *with* their own DNA is
> near-circular — `fnat`'s native contact set is defined from model 1 of that same bundle — so it measures
> deposition self-consistency, not what a docked structure can reach. Strip the DNA and re-dock. (ii) The
> gate writes `NA` when the *reference* has no native contacts (`fnat = rec/len(nat) if nat else nan`), so
> an all-`NA` file means the reference chain selection is wrong, not that the ensemble scored zero. Both
> scripts now raise `BrokenReference` on that. (iii) `c03` (c-Myb) has a broken holo reference — 0.00
> ceiling against an apo of 0.77 — and is excluded throughout.

## 19½.4 — The apo-fnat mechanism read-out

Strip BioEmu out of §19½.3 and one column stands on its own: **how well a deposited apo structure of the
protein reproduces the bound interface**. It splits the pair set cleanly.

| induced-fit-like (apo fails the gate) | apo pass | | selection-like (apo passes) | apo pass |
|---|---|---|---|---|
| c22 Trp repressor | 0% | | c13 Lac repressor | 56% |
| nhp6a NHP6A | 0% | | c01 Mu repressor C | 81% |
| c20 THAP1 | 5% | | c04 MeCP2 MBD | 84% |
| vnd VND/NK-2 | 30% | | c05 WRKY4, c09 FOXD3, c10 SarA, c11 Dead ringer, c17 YB-1 | 100% |
| rok Rok DBD | 35% | | | |

**[nmr]**

**This is the measurement Chapter 17 does not currently have.** Ch. 17 infers mechanism from the *sign of
the augmentation effect*; this measures the protein-side half of the same question directly, from deposited
structures, with no model in the loop. It is not a replacement — Ch. 17's revised discriminator is a DNA
property (ensemble bend fluctuation) and this is a protein property — but it is an independent axis that
the sign-based inference can be tested against rather than assumed.

**It agrees with the literature where the literature is clear.** Trp repressor: without tryptophan the HTH
collapses into the empty binding pocket, and the largest apo-vs-holo differences are in helix D, the
N-terminus of helix E and the intervening turn — the DNA-reading heads themselves. NHP6A: the free protein
carries an unstructured 17-residue basic N-terminal tail that wraps the major groove on binding, a standard
induced-fit example. Both land at 0%.

## 19½.5 — The pilots

Five of the thirteen pilots have a deposited apo structure of the same protein (95% cluster, no nucleic
acid): **engrailed, ets1, hsf, nfat, runx**. The other eight — csl, dux4, egr1, err, foxa, irf, lef1, tbp —
have none, so the read-out cannot be computed for them.

| pilot | #apo | apo fnat | apo pass | BioEmu | pass | best frame |
|---|---|---|---|---|---|---|
| nfat | 1 | 0.51 | 60% | 0.40 | 1% | 0.50 |
| engrailed | 2 | 0.81 | 100% | 0.31 | 0% | 0.46 |
| ets1 | 7 | 0.74 | 100% | 0.37 | 1% | 0.52 |
| hsf | 2 | 0.79 | 100% | 0.45 | 18% | 0.55 |
| runx | 8 | 0.73 | 100% | 0.29 | 5% | 0.54 |

**[nmr]**

**BioEmu falls below its own apo reference in 5 of 5 pilots.** Wherever the free state is demonstrably
binding-competent — real apo clearing the gate at 60–100% — the BioEmu ensemble passes at 0–18%. The
deficit found on the NMR pairs reproduces on the systems the pipeline actually trains on.

**runx must be reclassified, and how the data flagged it is the point.** Its eight apo entries disagree:
the six crystals pass at 100%, both NMR bundles (`1CMO` 43 models, `1CO1` 10) at **0%** (fnat 0.25–0.38).
The literature resolves it against the crystals. The Runt domain has a documented **S-switch** — loop L11
closed in the free protein, open on CBFβ/DNA binding — and ~25% of amide resonances are unobservable in
free Runt from conformational exchange, appearing only on DNA binding. The apo crystals used here are
precisely the *"structural switch and specifically bound chloride ions modulate DNA binding"* structures,
with chloride occupying the DNA-binding site: crystallised in a binding-competent state. **runx is
induced-fit-like**; `1CMO`/`1CO1` are the meaningful measurement. **Prefer NMR apo over crystal apo, and
treat a crystal/NMR disagreement as the signal that the free state is dynamic — not as noise to median
away.**

## 19½.6 — Honest limits

- **The read-out measures the DBD fold, not the protein.** The shared-core intersection trims exactly the
  flanking and terminal regions where induced fit most often lives — autoinhibitory modules, N-terminal
  arms, hinge helices. **ets1 is the worked example:** Ets-1 genuinely changes on DNA binding (helix HI-1
  of the autoinhibitory module unfolds; free Ets-1 sits in a rigid-inactive/flexible-active equilibrium),
  but that module is flanking and was trimmed away, so the 103 residues scored are the preformed ETS
  domain. "Selection-like" is right *for the fold* and the literature agrees — but "is the DBD preformed"
  and "does this protein undergo induced fit" are different claims, and the mechanism hypothesis is about
  the second.
- **Do not read it per family.** engrailed and vnd are both homeodomains and land at 100% and 30%. The same
  split appears in HMG-box, where NHP6A is induced-fit while the SOX2 DBD has been reported to bind by
  conformational selection.
- **No discrimination among the pilots.** All five available pilots read selection-like (runx after
  reclassification excepted), because the induced-fit end is not represented among the pilots that happen
  to have apo structures. The NMR pairs span both ends, so this is a **coverage** problem, not a method
  problem. lef1 — HMG-box, the family where NHP6A scores 0% — is among the eight with no apo structure.
- **Crystal apo ≠ solution apo** (runx), and **n is small**: 14 informative pairs, 5 pilots.
- NMR references have their own coordinate precision, so absolute fnat here is not numerically comparable
  to the pilots' crystal-referenced fnat — only internally comparable within each table.

## Figures — Chapter 19½

- *(create)* **The `%across` strip** — the 14 informative pairs ordered on the apo→holo axis, with the
  2–3% sampling-noise floor drawn as a band. The one figure that retires "BioEmu = apo."
- *(create)* **The three-way fnat panel** — holo-redocked / real apo / BioEmu, box per source, `FNAT_FLOOR`
  marked. Carries §19½.3 on its own.
- *(create)* **The mechanism split** — apo pass rate per pair, sorted, with the induced-fit/selection
  divide marked and the literature-confirmed cases (TrpR, NHP6A, RUNX1) flagged.
- **Existing:** `analysis/analyses/nmr_rmsd/plots/pair_pca_*_ca.png` (deposited-subspace PCA, per pair —
  PC1 position + residual), `pair_side_by_side_all_ca_core.png` (apo vs holo per pair, shared scales),
  `rmsd_overview_ca_core.png` (all 35 ensembles ranked). Caption the PCA panels with their variance
  fraction; PC2 holds only 5–10% and carries little.
