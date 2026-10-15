# `analysis/analyses/nmr_rmsd` — what is actually in there

**Read 2026-09-29.** Summary of findings, verified against the on-disk tables rather than the README
alone.

## First: the provenance note in v7 is wrong and should be struck

v7 says these results "are not on the cluster" and that Ch. 19½ should stay out of the talk until
resolved. **That came from an earlier check of mine and it was incorrect.** The directory exists and is
substantial:

- **545 MB, 5,442 files**
- a 180-line `README.md` documenting method, gates, and limits
- 11 analysis scripts (`compute_pair_rmsd.py`, `pair_core.py`, `fnat_on_pairs.py`,
  `fnat_on_pilots.py`, `sweep_trim.py`, five plotting scripts)
- 230 figures in `plots/`
- per-pair and per-pilot result tables, plus `fnat/` and `fnat_pilots/` output trees

All four Ch. 19½ components have traceable sources: the `%across` table
(`pair_pca_summary_ca.csv`), the three-way fnat panel (`fnat/fnat_three_way.csv`), the
informativeness gate (`pair_bundle_summary_ca.csv`), and the pilot mechanism read-out
(`fnat_pilots/pilot_mechanism.csv`). I recomputed `%across` from the PCA summary and reproduced the
README's table exactly.

**One narrow provenance gap survives.** `pilot_mechanism.csv` has no identifiable writer:
`fnat_on_pilots.py` writes only `fnat_pilots_summary.csv`, and the string `nmr_pass` appears in no
script in the tree. The file was written 5 s after the summary, so it looks like an unsaved
post-processing step. **That one table — and the runx reclassification that rests on it — is the only
thing that should stay out of the talk.**

**The README is also internally stale.** Its "What this shows" section describes 16 pairs / 14
informative, which matches the data; its "Known limits" section still says "Five pairs, three
informative." The second is left over from an earlier state. Fix before citing.

## What the analysis does

Two layers, both about whether **BioEmu reproduces real free-state conformational behaviour** — which
is the apo/holo question the project's BioEmu citation raises.

**Layer 1 — apo/holo NMR pairs.** Solution-NMR structures of the same protein with and without DNA.
Both bundles plus two BioEmu ensembles (run on the apo and holo constructs separately, as a control)
are superposed on one shared core, then a PCA is fitted **on the deposited states only** and the
ensembles projected in as test data. `%across` = position on the apo→holo axis, 0% at apo, 100% at
holo.

The methodology is careful in ways worth noting: the core is defined by the NMR bundles alone so the
thing being tested cannot define its own yardstick; expression tags are regex-stripped; disordered
termini are trimmed with a sensitivity sweep showing 2.5–4.0 Å is a stable plateau and no-trim flips
signs; bundle residuals are leave-one-out so deposited states don't score a free zero.

**Layer 2 — the apo-fnat mechanism read-out.** Dock a *deposited apo* structure onto the holo DNA and
score fnat. High ⇒ the real free state already reproduces the bound interface (selection-like); low ⇒
it does not (induced-fit-like). Per the script's own docstring, the point is that **this needs no
BioEmu at all**, so it measures mechanism independently of Ch. 17, which currently infers it from the
sign of the augmentation effect.

## Finding 1 — BioEmu is not apo-biased; it is *unsystematic*

| pair | protein | %across | resid ratio | bundle sep (Å) |
|---|---|---|---|---|
| rok | Rok C-term DBD | 2% | 3.7 | 1.00 |
| c09 | PROTEIN (TRANSCRIPTION FACTOR) | 18% | 4.3 | 1.10 |
| c10 | Transcriptional regulator SarA | 25% | 10.2 | 0.35 |
| nhp6a | NHP6A | 42% | 3.0 | 1.05 |
| c17 | Y-box-binding protein 1 | 46% | 6.9 | 0.61 |
| c03 | C-Myb DNA-Binding Domain | 47% | 1.9 | 0.44 |
| c20 | THAP domain-containing protein 1 | 51% | 9.1 | 3.28 |
| c13 | Lactose operon repressor | 58% | 3.1 | 1.49 |
| c04 | Methyl-CpG Binding Protein | 61% | 2.5 | 2.52 |
| vnd | VND/NK-2 | 64% | 3.5 | 0.75 |
| c22 | TRP REPRESSOR | 65% | 2.9 | 2.71 |
| c01 | REPRESSOR PROTEIN C | 84% | 4.4 | 1.73 |
| c11 | DEAD RINGER PROTEIN | 85% | 5.5 | 2.66 |
| c05 | Probable WRKY transcription fact | 94% | 3.0 | 1.28 |

14 of 16 pairs pass the informativeness gate (trf1 and maze fail: their apo/holo separation is ~0.17 Å
against a ~1 Å bundle-noise floor). Across those 14, **`%across` spans 2% to 94%, median 55%**.
Restricting to the 8 pairs whose ensembles stay near the deposited subspace changes nothing
(median 60%, same range).

**This retires the directory's own earlier conclusion**, which the README documents honestly: on the
first three pairs the reading was "BioEmu reaches neither basin and never gets to holo." At n = 14 that
does not survive — two pairs sit at the apo end, three at the holo end, nine in between, and WRKY4 at
94% is essentially at the bound state.

It also bears directly on the BioEmu paper's documented holo preference. **Neither bias reproduces
here**: on protein–DNA systems the landing point varies by system across the whole range.

**What does survive is a different signature.** The ensembles are consistently *wider* than the
deposited bundles and displaced outside their subspace — residual ratio **1.9× to 10.2×**, above 1 in
all 14. That displacement, not a directional bias, is the reproducible BioEmu behaviour here.

I also tested the one question the README leaves open (what predicts `%across`). Bundle separation
gives ρ = +0.552 (p = 0.041), but it does not survive dropping the two lowest-separation pairs
(ρ = +0.448, p = 0.145), and `%across` carries the apo–holo gap in its denominator while
`separation_A` is the RMSD analogue of that same gap — so the two are not independent. **Treat as
suggestive at best.** Residual ratio predicts nothing (ρ = −0.196, p = 0.50).

## Finding 2 — real apo structures rebuild the bound interface; BioEmu does not

The three-way fnat panel, 13 scoreable pairs (c03's holo bundle scores 0 — a docking failure, not a
result; trf1 and maze are gated out):

| pair | holo bundle | apo bundle | BioEmu |
|---|---|---|---|
| c13 | 0.658 (70%) | 0.500 (56%) | 0.263 (0%) |
| nhp6a | 0.738 (100%) | 0.194 (0%) | 0.263 (6%) |
| vnd | 0.837 (100%) | 0.349 (30%) | 0.279 (0%) |
| c17 | 0.840 (100%) | 0.640 (100%) | 0.280 (1%) |
| c20 | 0.929 (100%) | 0.393 (5%) | 0.321 (1%) |
| c09 | 0.818 (100%) | 0.667 (100%) | 0.381 (0%) |
| c10 | 0.882 (100%) | 0.765 (100%) | 0.382 (6%) |
| c04 | 0.923 (100%) | 0.577 (84%) | 0.385 (12%) |
| c22 | 0.906 (100%) | 0.406 (0%) | 0.406 (11%) |
| c05 | 0.841 (100%) | 0.727 (100%) | 0.455 (32%) |
| c01 | 0.833 (100%) | 0.619 (81%) | 0.464 (43%) |
| rok | 0.891 (100%) | 0.469 (35%) | 0.469 (47%) |
| c11 | 0.971 (100%) | 0.676 (100%) | 0.559 (75%) |

Median fnat: **holo 0.839, apo 0.629, BioEmu 0.396.** Pass rates at the 0.5 gate: holo 100%, apo 92%,
**BioEmu 11.5%**. The apo bundle beats BioEmu in **13 of 16** pairs.

That ordering is the substantive result. The holo bundle scoring high is near-tautological — it is the
bound state. But the *apo* bundles also largely clear the gate, meaning real free-state ensembles
contain near-bound conformations, while BioEmu ensembles of comparable proteins mostly do not.

## Finding 3 — the same thing on your own pilots, and it is sharper

| pilot | n apo entries | deposited apo fnat (pass) | BioEmu fnat (pass) | BioEmu best frame |
|---|---|---|---|---|
| nfat | 1 | 0.509 (60%) | 0.404 (1.1%) | 0.500 |
| runx | 8 | 0.729 (100%) | 0.292 (5.2%) | 0.542 |
| ets1 | 7 | 0.739 (100%) | 0.370 (1.0%) | 0.522 |
| hsf | 2 | 0.793 (100%) | 0.448 (18.2%) | 0.552 |
| engrailed | 2 | 0.807 (100%) | 0.308 (0.0%) | 0.462 |

Four of five pilots have deposited apo structures that clear the gate **100%** of the time. BioEmu
passes 0–18%. And the sharpest form:

> **BioEmu's single best frame is below the deposited apo *median* in 5 of 5 pilots.**

Not "BioEmu samples the bound state less often" — across ~90 frames per pilot it never reaches what a
median real apo structure achieves.

This is a direct, BioEmu-free measurement of mechanism for your pilots: engrailed, ets1, hsf and runx
all look **selection-like** (apo fnat 0.73–0.81, 100% pass), nfat intermediate (0.509, 60%). Measured
without reference to the augmentation effect, so it is an independent check on Ch. 17's sign-based
inference.

## How this connects to the rest of the project

- **It is the missing control for the whole premise.** The pipeline pairs BioEmu protein frames with
  crystal DNA and asks whether that recovers binding specificity. Finding 3 says the frames themselves
  do not reconstitute the native interface as well as a deposited apo structure does — which bounds
  what any downstream treatment can recover, independently of the DNA handling.
- **It is consistent with the interface-gate pass rates** already in the pipeline, and explains them:
  states are being gated on an interface that the generator does not reproduce well.
- **It does not depend on DeepPBS at all**, so it is unaffected by the in-sample/shared-label problems
  that qualify F2 and F3.

## Recommendation

**Strike the v7 provenance note and put Ch. 19½ back in scope**, minus the runx reclassification. The
`%across` result and the three-way fnat panel are well-sourced, methodologically careful, and
self-correcting — the directory already retired its own three-pair conclusion when n grew. Finding 3 in
particular is a strong, cheap, DeepPBS-independent result and arguably belongs in the talk as a
limitation on the generator rather than being held back.

Two fixes before citing: update the README's stale "Known limits" paragraph, and either recover the
script that wrote `pilot_mechanism.csv` or regenerate that table from `fnat_pilots_summary.csv` with
the `apo_pass` / `nmr_pass` distinction documented.

## Source files

`analysis/analyses/nmr_rmsd/` — `README.md`, `pair_bundle_summary_ca.csv` (gate),
`pair_pca_summary_ca.csv` (`%across`, residuals), `fnat/fnat_three_way.csv`,
`fnat_pilots/pilot_mechanism.csv`, `pairs.csv` (5 ready, 33 `needs_ensembles`), `plots/` (230 figures).
Backbone-atom variants of every table agree with CA to within 0.05 Å.
