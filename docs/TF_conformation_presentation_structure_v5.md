# Presentation structure v5 — ~60 minutes

**Supersedes v4.** Different purpose, so a different shape. v4 was built to defend a position: it led
with a verdict table, spent Act 4 on a retraction, and treated the negatives as the spine. This one is
built to **bring the room up to date and get direction out of them**. The spine is now the apo/holo NMR
experiment — the newest and most generative thing in the project — and everything else is sized relative
to it.

**Framing rule:** this is an exploration, not a paper draft. Say so in the first minute and keep saying
it. Numbers are presented as "here is what we see on 16 systems", never as "here is the result".

**Tone rules**

- Negatives get one slide each, stated plainly, then we move on. No act is built around one.
- Every section ends on a question or an option, not a conclusion.
- `n` on every slide with a statistic. Say "16 pairs", "5 pilots", "12 pilots" out loud each time.
- When a number is soft, say it is soft. The room is going to propose directions; they need to know
  which numbers can bear weight.

**Title:** _"What apo ensembles do and don't know about the bound state — an exploration on 16 NMR pairs
and 13 TFs."_

---

## Part 1 — The arc

| act | question                                              | slides        | min     |
| --- | ----------------------------------------------------- | ------------- | ------- |
| 0   | Where we left off in August, in three slides          | S1–S4         | 4       |
| 1   | Do apo protein frames even rebuild the DNA interface? | S5–S8         | 7       |
| 2   | **The NMR pair experiment — how it works**            | S9–S14        | 12      |
| 3   | **Where BioEmu ensembles actually sit**               | S15–S20       | 12      |
| 4   | What this is telling us about TFs and DNA             | S21–S25       | 9       |
| 5   | What the augmentation experiment taught (briskly)     | S26–S29       | 6       |
| 6   | Where I want your input                               | S30–S32       | 10–12   |
|     |                                                       | **32 slides** | **~60** |

Act 2 + Act 3 = 24 minutes on the NMR work, which is 40% of the talk and the right weight for the thing
that is newest and least settled.

---

## Part 2 — Slide by slide

Format: **Say** = the one sentence that has to land · **Figure** = exact repo path, or `NEW` with the table
it comes from · **Notes** = what to have in your head.

Provenance marks on numbers: **✓** recomputed 2026-10-06 from the named table · **○** quoted from v4 or
the August deck, not re-verified · **?** needs re-deriving before the talk.

---

### Act 0 — Where we left off (4 min)

**S1 — Title + the frame** [1]
**Say:** "Two months ago I showed you a pipeline and a hypothesis. Today I have a new measurement, it is
pointed at a question none of us had asked, and I want your help deciding where it goes next."
**Figure:** none.
**Notes:**

- Set the exploration framing explicitly: nothing here is a paper, several numbers will move.
- Promise the directions slide at the end so the room saves its ideas for it (or interrupts — say which
  you prefer; interruptions are better in Acts 2–3).

**S2 — The project in one slide** [1]
**Say:** A co-crystal is one conformation; specificity is a distribution over sequences. We asked whether
feeding DeepPBS an ensemble instead of one structure changes what it predicts.
**Figure:** reuse **old deck S5** (goal schematic) or **old S29** (5-stage pipeline strip).
**Notes:** thirty seconds on DeepPBS, thirty on BioEmu. Do not re-teach the pipeline; S3 handles it.

**S3 — The pipeline, compressed** [1]
**Say:** BioEmu → dock to crystal DNA → minimize → fnat gate → featurize → DeepPBS, with surviving-state
counts per pilot.
**Figure:** reuse **old deck S29/S41** (stage strip) + `analysis/figures/F2_passrate_bars.png`.
**Notes:**

- "If you were here in August this is the 55 slides you already sat through. Ask and I'll open any."
- 13 pilots, ~90–100 frames each, 12 with full frozen×relaxed grids.

**S4 — What changed since August** [1]
**Say:** Three things: a new measurement on deposited NMR pairs, the augmentation experiment finished and
interpreted, and one label bug that cost us TBP and bought us NDT80.
**Figure:** none — a three-row list, deliberately bland. The content is in the acts.
**Notes:** do not preview conclusions here. This is a table of contents, not a verdict table.

---

### Act 1 — Do apo frames rebuild the interface? (7 min)

**S5 — The question nobody had asked** [1.5]
**Say:** The whole pipeline rests on an assumption: that a free-state protein conformation, placed on the
bound DNA, reconstitutes something like the native interface. We never tested it.
**Figure:** reuse **old deck S7/S8** (apo vs holo schematic).
**Notes:**

- BioEmu is trained on and samples apo. DeepPBS reads interface geometry. The gap between those two
  facts is the whole talk.
- This is the slide that reframes an engineering detail as a biological question.

**S6 — fnat, in thirty seconds** [1]
**Say:** Fraction of native protein–DNA contacts recovered; 0.5 is the gate the pipeline already uses.
**Figure:** reuse **old deck S72** (fnat definition, DockQ reference) + `analysis/figures/F1_fnat_distributions.png`.
**Notes:** iRMSD and fnat co-vary at ρ = −0.90 ○ (old S75), so one number suffices for this talk.

**S7 — The test, and three things to compare** [1.5]
**Say:** Take solution-NMR structures of the _same protein_ with and without DNA, dock each state onto the
bound DNA, and score fnat. No BioEmu, no DeepPBS, no training anywhere in this measurement.
**Figure:** `NEW` — schematic, three rows: holo bundle with DNA stripped and re-docked (positive control) ·
deposited apo bundle · BioEmu ensemble. Hand-drawn is fine.
**Notes:**

- The holo re-dock control is what makes the other two readable: it measures what the docking procedure
  costs you when the conformation is already right.
- Flag honestly: the control fails for 3 of 16 pairs (c03 is the worst — median fnat 0.000 ✓), and those
  pairs are therefore not interpretable, not negative.

**S8 — The three-way comparison** ⭐ [3]
**Say:** Re-docked holo states clear the gate almost always, deposited apo states clear it about two-thirds
of the time, BioEmu frames clear it about a fifth of the time.
**Figure:** `NEW` **N1** — violin or strip per source with the 0.5 floor drawn. Source:
`analysis/analyses/nmr_rmsd/fnat/*_fnat_{holo_states,apo_states,bioemu}.csv`.

| source                               | pooled median fnat | pass ≥ 0.5 | n states |
| ------------------------------------ | ------------------ | ---------- | -------- |
| holo bundle, DNA stripped, re-docked | 0.844 ✓            | 90% ✓      | 309      |
| deposited apo bundle                 | 0.636 ✓            | 71% ✓      | 384      |
| BioEmu ensemble                      | 0.393 ✓            | 23% ✓      | 1528     |

Restricted to the 13 pairs whose holo control itself passes: 0.859 / 98%, 0.553 / 62%, 0.381 / 18% ✓.

**Notes:**

- **?** v4's S10 quoted 0.839 / 0.629 / 0.396 and 100% / 92% / 11.5%. Those do not reproduce from the
  saved per-state tables under either pooling or per-pair medians. Re-derive and pick one aggregation
  before this slide goes in — it is the number the room will remember.
- Say what each row buys: row 1 = the dock is not the problem; row 2 = the gate is not rejecting
  free-state conformations as a class; row 3 = so the gap is specific to the generated ensembles.
- End forward, not down: "that is a measurement about ensembles, and it is also a measurement about
  proteins — which is what Act 4 is about."

---

### Act 2 — How the NMR pair experiment works (12 min)

This act is method, and it is worth the time: it is the part the room can poke holes in, and the part that
makes Act 3 believable. Invite interruptions here explicitly.

**S9 — The pair set** [2]
**Say:** 16 apo/holo pairs mined from the PDB — solution-NMR protein+DNA entries matched to protein-only
solution-NMR entries in the same 95% sequence cluster — plus 33 more already identified and queued.
**Figure:** `NEW` **N2** — pair roster table from `analysis/analyses/nmr_rmsd/pairs.csv` +
`fnat/fnat_three_way.csv`: pair id, protein, family, apo/holo PDB, #states each side.
**Notes:**

- Filters: ≥10 models both sides, ≥84% coverage of the shorter sequence, DNA-binding regulators.
- Families span HMG-box, homeodomain, forkhead, HTH, MYB/SANT, ARID, MBD, WRKY, THAP, cold-shock,
  MarR/winged-helix, Telobox — four of them map onto existing pilots (c09↔foxa, c12↔lef1, c32↔nfat,
  c08↔trf1), which is what ties this back to the benchmark.
- "apo" means no nucleic acid in the entry — it does **not** mean vetted: some carry peptides or ligands,
  differ in oligomeric state, or cover different domain boundaries. Say this out loud once.

**S10 — Why you cannot just use the deposited RMSDs** [2]
**Say:** Each bundle is a different construct with different disordered ends, so per-entry numbers land on
different atom sets and a histogram of them measures construct length as much as conformation.
**Figure:** reuse `analysis/analyses/nmr_rmsd/plots/rmsd_apo_holo_ca.png` (per-entry, ragged) beside
`plots/rmsd_apo_holo_ca_core.png` (shared core) — the before/after of the fix.
**Notes:**

- One pair reads 9.47 Å untrimmed and 0.81 Å on its core ○. That single comparison is the argument.
- This is the methods slide that earns the rest; do not cut it.

**S11 — Building one honest atom set per pair** ⭐ [3]
**Say:** One core per pair, shared by all four members — apo bundle, holo bundle, and both BioEmu
ensembles — defined by the deposited structures only.
**Figure:** `NEW` **N3** — three-step schematic: (1) four-way sequence intersection via difflib,
(2) minus expression tags matched by regex, (3) minus disordered termini trimmed inward while either
bundle's intra-bundle spread exceeds 3.0 Å. Annotate with core sizes from
`pair_core_summary_ca.csv`.
**Notes:**

- Step 2 is not cosmetic: 5zux/5zuz share a C-terminal `LEHHHHHH` that the intersection keeps because
  _both_ members have it; it would contribute pure noise to every RMSD.
- Only ends are trimmed — interior mobile segments stay in (NHP6A keeps an internal loop reaching ~7 Å).
- **The core is defined by the two NMR bundles only.** The ensembles are what is being tested and must
  not define the yardstick. If one objection lands tonight, let it be this one.

**S12 — Does the pair separate its own states? (the gate)** [2]
**Say:** A pair earns its place only if the apo–holo separation clears the bundles' own internal
heterogeneity — NMR spread is precision as much as it is dynamics.
**Figure:** `NEW` **N5** — within vs between bundle RMSD per pair with bootstrap CIs and the
ratio threshold. Source: `pair_bundle_summary_ca.csv`.

| pair  | within | between | separation [95% CI]  | ratio | verdict       |
| ----- | ------ | ------- | -------------------- | ----- | ------------- |
| nhp6a | 0.50   | 1.55    | +1.05 [+0.96, +1.12] | 3.08  | informative ○ |
| vnd   | 0.41   | 1.15    | +0.75 [+0.68, +0.81] | 2.83  | informative ○ |
| rok   | 0.75   | 1.75    | +1.00 [+0.89, +1.11] | 2.33  | informative ○ |
| trf1  | 1.11   | 1.28    | +0.17 [+0.09, +0.23] | 1.15  | gated out ○   |
| maze  | 0.90   | 1.07    | +0.17 [−0.03, +0.31] | 1.19  | gated out ○   |

**Notes:**

- 14 of 16 pass. trf1's CI excludes zero only because n is large — 0.17 Å on a 1.1 Å noise floor.
- Bootstrap unit is states, not state-pairs; self-comparisons from resampling are dropped, not scored 0.
- The new pairs widen the dynamic range a lot: THAP1 and TrpR separate their bundles by 4.3–4.4 Å against
  1.1–1.8 Å for the original three ○. That range is what was missing for asking whether BioEmu tracks a
  transition at all.

**S13 — The coordinate system: a deposited-state subspace** ⭐ [2]
**Say:** Superpose all four members on the shared core, fit a PCA on the deposited states only, then
project the ensembles in as test data — so the conformations that were actually observed define the axes
and the thing being tested never influences them.
**Figure:** `analysis/analyses/nmr_rmsd/plots/pair_pca_nhp6a_ca.png` (clean single pair).
**Notes:**

- Two readouts per structure, both in Å and directly comparable to an RMSD: position on PC1/PC2, and
  **residual** — the part of the structure the subspace cannot reproduce.
- Bundle residuals are leave-one-out, so states that defined the space don't score a free zero.
- This replaced a classical-MDS embedding that retained only 35–67% of the positive eigenvalue mass ○
  (RMSD after superposition is not a Euclidean metric). Mention in one sentence; it is a credibility beat.

**S14 — The subspace validates itself** [1]
**Say:** For the three original informative pairs a single PC dominates and it _is_ the apo→holo
direction; for the two gated-out pairs no direction dominates. The gate and the geometry agree without
being told to.
**Figure:** `plots/pair_pca_all_ca.png`.

| pair  | PC1 share of deposited variance ✓ | k for 90% ✓ |
| ----- | --------------------------------- | ----------- |
| nhp6a | 78%                               | 4           |
| vnd   | 75%                               | 4           |
| rok   | 70%                               | 7           |
| trf1  | 25%                               | 15          |
| maze  | 32%                               | 10          |

**Notes:** c11 (ARID) and c20 (THAP1) are even sharper — 88% on PC1, k = 2 ✓. Those are the pairs where
"% of the way from apo to holo" means the most.

---

### Act 3 — Where BioEmu ensembles actually sit (12 min)

**S15 — The question, now answerable** [1]
**Say:** With a coordinate system built from observed conformations, we can ask where a generated ensemble
lands: at the free state, at the bound state, somewhere between, or off the map entirely.
**Figure:** `plots/pair_pca_rok_ca.png` and `plots/pair_pca_c05_ca.png` side by side — the two extremes
(rok parks at apo, WRKY4 sits essentially at holo).
**Notes:** define `%across` once: 0% at the apo bundle, 100% at holo. And `resid ratio`: ensemble residual
over the bundles' own LOO residual — how far outside the observed subspace the ensemble sits.

**S16 — Position along the apo→holo axis, all 14 informative pairs** ⭐ [3]
**Say:** It spans 2% to 94%, median 55%. BioEmu is not systematically apo-biased — where its ensemble
lands varies by system across the entire range.
**Figure:** `NEW` **N4** — lollipop/dot plot, pairs sorted by `%across`, marker shaded by residual ratio,
0% and 100% anchored. Source: `pair_pca_summary_ca.csv` + `pair_projection_ca.csv`.

| pair  | protein            | %across ○ | resid ratio ○ |
| ----- | ------------------ | --------- | ------------- |
| rok   | Rok C-term DBD     | 2%        | 4             |
| c09   | FOXD3 / Genesis    | 18%       | 4             |
| c10   | SarA               | 25%       | 10            |
| nhp6a | NHP6A              | 42%       | 3             |
| c17   | YB-1               | 46%       | 7             |
| c03   | c-Myb              | 47%       | 2             |
| c20   | THAP1              | 51%       | 9             |
| c13   | Lac repressor      | 58%       | 3             |
| c04   | MeCP2 MBD          | 61%       | 2             |
| vnd   | VND/NK-2           | 64%       | 4             |
| c22   | Trp repressor      | 65%       | 3             |
| c01   | Mu repressor C     | 84%       | 4             |
| c11   | Dead ringer (ARID) | 85%       | 5             |
| c05   | WRKY4              | 94%       | 3             |

**Notes:**

- **This is the slide that supersedes what I would have told you at n=3.** On nhp6a, rok and vnd alone the
  reading was "BioEmu reaches neither basin and never gets to holo". At n=14 that does not survive. Say it
  as a lesson about small n, cheerfully — it is the best argument for the 33 queued pairs.
- Restricting to the 8 pairs whose ensembles stay near the subspace (ratio < 4) changes nothing: median
  60%, same 2–94% range ○.
- Two pairs at the apo end, three at the holo end, nine in between.
- c05/WRKY4 at 94% is the case the three-pair set said did not occur.

**S17 — What _does_ reproduce: the ensembles are wider** ⭐ [2.5]
**Say:** In every informative pair, the ensemble sits further outside the deposited subspace than the
bundles do — 1.9× to 10.2× the bundles' own leave-one-out residual.
**Figure:** `NEW` **N6** — per-pair bar of ensemble residual / bundle LOO residual, ordered, with 1× line.
Source: `pair_pca_summary_ca.csv` (recomputed ratios ✓: c03 1.9 · c04 2.5 · c22 2.9 · nhp6a 3.0 ·
c05 3.0 · c13 3.1 · vnd 3.5 · rok 3.7 · c09 4.3 · c01 4.4 · c11 5.5 · c17 6.9 · c20 9.1 · c10 10.2).
**Notes:**

- Direction-free and system-independent: this, not an apo bias, is the reproducible signature.
- Two readings, and I don't know which is right — say that: (a) the generator over-disperses; (b) the
  deposited bundles under-represent real breadth because NMR refinement regularizes toward a mean.
  Distinguishing them is on the directions slide.
- BioEmu's own paper reports a _holo_ preference on cryptic pockets (85% vs 49%) ○; neither that bias nor
  an apo bias reproduces here. Worth one sentence — the room will ask.

**S18 — The control that makes this trustworthy** [2]
**Say:** Two pairs give BioEmu byte-identical input sequences on the apo and holo sides, so their
disagreement is pure sampling noise — and it is about 2–3% of the apo→holo axis.
**Figure:** `NEW` — small two-panel PC1 density overlay for nhp6a and vnd (apo-seq vs holo-seq). Source:
`pair_projection_ca.csv`.
**Notes:**

- nhp6a and vnd: identical 93- and 77-residue inputs; PC1 medians differ by 0.17 and 0.19 Å; KS cannot
  separate them (p = 0.64, 0.33); IQRs match ○.
- So nhp6a's and vnd's failure to reach holo is 15–25× the sampling noise ○, and rok parking at apo is
  within ~1 noise unit of the bundle.
- Every other construct difference in the set is flanking-only; no pair has an interior substitution ○.
- One real caution: a tag can be dropped from the measurement core but not from the _generation_. maze's
  two ensembles differ 1.7× in residual and the worse one is the **untagged** 50-mer — truncation, not the
  His-tag ○. Strip tags from the input sequence, not just the analysis.

**S19 — How much of this is the trimming?** [1.5]
**Say:** Swept the spread cutoff from 2 to 5 Å plus a no-trim control; 2.5–4.0 Å is a stable plateau for
every pair and the default 3.0 Å sits in it.
**Figure:** `analysis/analyses/nmr_rmsd/plots/trim_sweep_ca.png`.
**Notes:**

- vnd's negative delta holds at −0.07 to −0.18 Å across the plateau ○ — the sign is not a trimming artifact.
- 2.0 Å is degenerate: the core collapses to 7 positions (nhp6a) and 4 (vnd) ○.
- No-trim flips signs (vnd −0.10 → +1.59; nhp6a +0.34 → −0.41) ○ — with 12–20 Å of disordered tail left
  in, you are measuring tails. That is the control that justifies trimming at all.

**S20 — And on our own pilots** [2]
**Say:** Same measurement, our 13 pilots, using whatever apo structures exist for them: BioEmu's single
best frame out of ~90 is below the deposited apo median in 5 of 5 pilots that have apo entries.
**Figure:** `NEW` **N7** — per-pilot: BioEmu fnat distribution with its max marked, versus each deposited
apo source (NMR bundles and apo crystals shown differently). Source:
`analysis/analyses/nmr_rmsd/fnat_pilots/fnat_pilots_summary.csv`.

| pilot     | deposited apo median ✓  | BioEmu median ✓ | BioEmu best ✓ |
| --------- | ----------------------- | --------------- | ------------- |
| engrailed | 0.769 (2JWT, 25 states) | 0.308           | 0.462         |
| ets1      | 0.696 (1R36, 25 states) | 0.370           | 0.522         |
| hsf       | 0.862 (2LDU, 20 states) | 0.448           | 0.552         |
| nfat      | 0.510 (1NFA, 10 states) | 0.404           | 0.500         |
| runx      | 0.250 (1CMO, 43 states) | 0.292           | 0.542         |

**Notes:**

- Not a sampling-efficiency story: the best of ~90 frames does not reach the median of 10–43 deposited
  states. The conformations are not in the distribution.
- **runx is the interesting exception and needs care:** its NMR apo bundles score _low_ (0.250/0.375,
  0% pass) while its seven apo _crystal_ structures score 0.667–0.875 ✓. That is either real
  solution-vs-lattice disagreement or an alignment problem on 1CMO. Keep it in backup; `pilot_mechanism.csv`
  is mid-regeneration.
- ets1 also has 7 apo entries — the best-covered pilot, and the obvious next pair to generate (c02,
  1r36/2stt, needs only the apo side).

---

### Act 4 — What this is telling us about TFs and DNA (9 min)

This act is the biology, and it is where the room's expertise actually is. Keep it open-ended.

**S21 — A model-free read-out of binding mechanism** ⭐ [2]
**Say:** Strip the generator out entirely and this becomes a structural measurement on deposited data: how
much of the DNA-binding interface is already formed in the free state?
**Figure:** reuse `analysis/figures/M1_apo_holo_mechanism.png` + the apo column of **N1**.
**Notes:**

- High apo fnat ⇒ the interface is preformed ⇒ conformational selection. Low apo fnat ⇒ it is built on
  binding ⇒ induced fit. One number, no model, no training.
- It splits the set: c03/c-Myb 0.765, c10/SarA 0.765, c05/WRKY4 0.727, c11/ARID 0.676, c09/FOXD3 0.667 at
  the preformed end; nhp6a 0.194, vnd 0.349, c20/THAP1 0.393, c22/TrpR 0.406 at the other ✓.
- It agrees with the literature where the literature is clear: TrpR and NHP6A are both documented
  induced-fit; MeCP2, WRKY4 and FOXD3 come out preformed ○.

**S22 — The two mechanisms, and why DNA shape is the other half** [2]
**Say:** Conformational selection and induced fit are statements about _both_ partners — and we have the
DNA side measured too.
**Figure:** reuse **old deck S94** (selection vs induced fit schematic) + **old S93** (the lab's
intrinsic-vs-induced MGW framing, Jiang et al. Biophys J 2026).
**Notes:**

- Frame the pairing explicitly: apo fnat reports whether the _protein_ is preformed; MGW fluctuation
  reports whether the _DNA_ is reshaped. Nobody has crossed those two axes on the same systems yet.
- That crossing is a directions-slide item, and it is the one I would most like to be talked into.

**S23 — What the DNA does in our ensembles** [2]
**Say:** Minimized ensembles keep the crystal's intrinsic minor-groove-width profile, and relaxing the DNA
adds fluctuation without moving the mean.
**Figure:** `analysis/analyses/dna_relax/figures/crystal_vs_ensemble_mgw.png` +
`engrailed_mgwfl_exemplar.png` (or reuse **old S97/S98/S103**).
**Notes:**

- Engrailed is the clean exemplar: MGW essentially identical frozen vs relaxed ○ — intrinsic shape held.
- pyCurves throughout (Jinsen's), two axis conventions, 13 pilots.
- Keep this slide descriptive. The correlation between DNA fluctuation and augmentation benefit is the
  thing we stopped claiming — that is S29, one line, later.

**S24 — AF3 vs an ensemble, on DNA flexibility** [1.5]
**Say:** AF3 gives one structure and it reads as nearly rigid DNA; the ensembles show 4–5× more
minor-groove fluctuation across every pilot.
**Figure:** `analysis/analyses/dna_relax/figures/af3_vs_ensemble_mgwfl.png` (+ `D1_diversity.png` or old
S23 for the protein-side diversity point).
**Notes:**

- AF3 MGW-FL ~0.09–0.24 Å vs ensemble medians 4.6×/4.9× higher, paired p = 2.9e-3 / 2.2e-4 ○.
- Reproduces the lab's own Biophys J Fig 3 on our systems — good place to say the pipeline agrees with
  something external.
- This answers "why BioEmu and not AF3?" before it is asked.

**S25 — A published test case I can actually be wrong about** ⭐ [1.5]
**Say:** Sox2 has a known right answer and a known wrong one — 2LE4 (2011) supported induced fit, 9QPF
(Orsetti et al., NAR 2025) finds a three-helical bundle identical to the bound state at 0.8 Å backbone
RMSD, and concludes conformational selection.
**Figure:** `NEW` — two apo structures superposed on the bound Sox2 structure, side by side. PyMOL,
`scripts/structure_viz/pymol_lib.py`.
**Notes:**

- 2LE4's RDC Q-factor is 0.76 against 0.08 for the refined structure ○.
- Can the apo-fnat read-out tell them apart from interface geometry alone? Half a day of work, queued.
- **Pre-register the likely failure, out loud:** the disorder trim removes about half the native contact
  set, and the helix that distinguishes the two structures makes only 2 of 72 DNA contacts ○. If both
  pass, that is a quantified blind spot — the read-out scores the ordered half of the interface — and
  that is worth knowing about a method. Saying this before the result is the most useful thing on this
  slide.

---

### Act 5 — What the augmentation experiment taught (6 min)

Deliberately brisk. Four slides, no act structure, and it ends on the fix.

**S26 — The experiment, and the one design fact that explains everything** [1.5]
**Say:** Baseline vs augmented, paired by seed, differing only in the data directory — and **every
conformation of a protein carries that protein's single experimental PWM.**
**Figure:** reuse **old deck S81** (paired-training design).
**Notes:** flag the shared label now, in one sentence. It is the whole interpretation two slides later.

**S27 — Conformation matters more than training noise** [1.5]
**Say:** Across 701 conformations, the spread in predicted PWM agreement across _conformations_ exceeds
the spread across independently trained _models_.
**Figure:** reuse **old S88** (how predictions differ) or `analysis/analyses/struct_pwm/figures/struct_pwm_context.png`.
**Notes:**

- Conformation SD 0.235 vs seed SD 0.16, n = 701 ○ **?** — re-derive from
  `analysis/analyses/conf_vs_crystal_pwm/per_state_all_conditions.csv` before use.
- Exemplar logos if time: EGR1's worst state loses the poly-G signature; TBP's best is indistinguishable
  from the crystal ○ (`struct_pwm/figures/struct_pwm_*.png`).

**S28 — What augmentation actually bought** ⭐ [2]
**Say:** The model became markedly more consistent across conformations and only slightly more accurate on
complexes it had never seen — which is exactly what we asked for when we gave every conformation the same
label.
**Figure:** `NEW` **N8** — paired baseline→augmented across-conformer SD per pilot (12 slopes), plus a
single panel for the cross-pilot delta. Source: `conf_vs_crystal_pwm/per_tf_all_conditions.csv` and
`crosspilot_deltas.csv`.
**Notes:**

- Across-conformer SD 0.247 → 0.202, lower in 10 of 12 pilots, Wilcoxon p = 0.034 ✓. HSF1 is the lone
  exception (0.392 → 0.518 ✓).
- Accuracy on a genuinely unseen complex: **+0.025** ✓ — cross-pilot design, 660 donor×target×seed rows,
  where pilot X's conformations are scored under donor pilot Y's checkpoints.
- Volunteer the in-sample structure rather than waiting to be asked: within-pilot the same quantity is
  +0.283 ✓, and all states are in-sample for their own augmented checkpoint, verified per entry. But the
  inputs are **not** near-copies — median distance to own crystal is 0.83 of the distance unrelated
  training entries sit at, and for err (1.07) and nfat (1.34) the states are _further_ than unrelated
  entries ✓. So this is label-side memorization with structurally distinct inputs.
- **?** v4's S18 quoted 0.235 → 0.180, 11/12, p = 0.027 — a different aggregation unit than the per-pilot
  table gives. Pick one before the talk.
- Close forward: the objective maps every conformation to one motif, so it cannot learn conformational
  sensitivity by construction. The fix is per-frame labels, which is a directions item.

**S29 — One thing we stopped claiming** [1]
**Say:** In August I showed a correlation between DNA bend fluctuation and augmentation benefit. When I
rebuilt bend from local base-pair-step geometry, every reconstruction came out flat, so I can't say what
the correlating quantity is and I've stopped using it.
**Figure:** reuse **old deck S99** with a "withdrawn" overlay — or no figure at all, which is faster.
**Notes:**

- One slide, 60 seconds, no defensiveness and no detail. The five-beat version and the full axis
  enumeration go to backup; pull them out only if asked.
- The useful residue, stated positively: at n = 12 with seed noise at two-thirds of the between-pilot
  spread, this design cannot resolve effects of this size — which is why the NMR pairs (n = 16 → 47,
  no training required) are where the effort went.

---

### Act 6 — Where I want your input (10–12 min)

**S30 — What's running right now** [1]
**Say:** Three things in flight, so you know what not to suggest.
**Figure:** none — a three-row table.
**Notes:**

- **NDT80 pilot**, 5 seeds training as of this morning, eval queued. New Pfam family (PF05224), 5
  within-family test anchors, and it exists because a label audit found TBP's PWM label pointed at an
  NDT80 matrix — so label, structure and anchors finally agree. One line on the bug, as infrastructure
  working rather than as a confession.
- **33 more apo/holo pairs** identified and filtered; each costs two BioEmu runs.
- **Sox2 blind test**, half a day.

**S31 — Six directions, and I want you to argue about them** ⭐ [8–10]
**Say:** Put it up and stop talking.
**Figure:** none. Six numbered items, one line each.

1. **Is the width signature a generator artifact or a deposition artifact?** Ensembles sit 1.9–10.2×
   outside the deposited subspace. NMR refinement regularizes toward a mean; BioEmu may over-disperse.
   What would separate those — RDC back-calculation on ensemble members? Comparison against MD on the
   same construct?
2. **Folded vs unfolded.** The free state really does contain partially unfolded members, and NMR
   determination selects for a determinable fold. Some of the fnat gap may be BioEmu correctly sampling a
   disordered tail. How would you split "generator is wrong" from "metric measures foldedness"?
3. **Which families should the next 33 pairs prioritize?** The set is choosable and four pairs already
   overlap our pilots. Breadth across folds, or depth on families where mechanism is contested?
4. **Cross the protein and DNA axes.** We have apo fnat (is the protein preformed?) and MGW fluctuation
   (is the DNA reshaped?) on overlapping systems and have never plotted one against the other. Is that
   the interesting 2×2, and what would each quadrant predict?
5. **Per-frame labels for augmentation — whose labels?** Self-labelling with DeepPBS is cheap but
   circular; physics-based scoring and SELEX-derived labels are both expensive. A λ-blend sweep between
   the experimental PWM and each conformation's own prediction is the cheap test. Better label source?
6. **Is "is the DBD preformed?" a paper on its own?** It needs no model and no training, it reproduces
   known induced-fit cases, and at n = 47 it would be a survey. Or is it a methods section inside the
   augmentation story?

**S32 — Backup index** [0]
**Say:** nothing; leave S31 up. This slide exists so you can jump.

---

## Part 3 — Backup slides

| slide                                                                                 | pulled out when                           |
| ------------------------------------------------------------------------------------- | ----------------------------------------- |
| Old deck 25–80 (full methods)                                                         | any pipeline question                     |
| Old deck 9–23 (AF3 vs BioEmu per TF)                                                  | "why BioEmu and not AF3?"                 |
| The five-beat bend-IQR withdrawal + full axis enumeration                             | "what happened to the August mechanism?"  |
| Power argument: seed noise = ⅔ of between-pilot spread, ~60 pilots for 80% power ○    | "why not more pilots?"                    |
| Pseudoreplication aside (old S86)                                                     | "why n = 5 and not 650?"                  |
| TBP label audit detail (`output/stage7_eval/tbp_label_audit/`)                        | "what happened to TBP?"                   |
| runx NMR-vs-crystal apo disagreement                                                  | "is the apo read-out robust?"             |
| Sox2 interface breakdown, 72 contacts by segment                                      | "why would both pass?"                    |
| Per-pair PCA panels, all 16 (`plots/pair_pca_*_ca.png`)                               | "show me pair X"                          |
| Trim-survival / core sizes (`pair_core_summary_ca.csv`)                               | "what does trimming cost?"                |
| Near-duplicate and leakage audits (`crystal_leakage.json`, `nearduplicate_ratio.csv`) | "isn't +0.025 leakage?"                   |
| Stage-3 minimization set, 20 panels (`analyses/stage3/figures/`)                      | "do the minimized structures make sense?" |

---

## Part 4 — Figure manifest

**Exists, use as-is** (paths relative to repo root)

| id                           | path                                                                                                                           |
| ---------------------------- | ------------------------------------------------------------------------------------------------------------------------------ |
| pipeline/gate                | `analysis/figures/F1_fnat_distributions.png`, `F2_passrate_bars.png`                                                           |
| per-entry vs core RMSD       | `analysis/analyses/nmr_rmsd/plots/rmsd_apo_holo_ca.png`, `rmsd_apo_holo_ca_core.png`                                           |
| single-pair PCA              | `analysis/analyses/nmr_rmsd/plots/pair_pca_{nhp6a,rok,c05,vnd,c22,c04}_ca.png`                                                 |
| all-pair PCA                 | `analysis/analyses/nmr_rmsd/plots/pair_pca_all_ca.png`                                                                         |
| trim sweep                   | `analysis/analyses/nmr_rmsd/plots/trim_sweep_ca.png`                                                                           |
| bundle/ensemble side-by-side | `analysis/analyses/nmr_rmsd/plots/pair_side_by_side_all_ca_core.png`                                                           |
| mechanism                    | `analysis/figures/M1_apo_holo_mechanism.png`                                                                                   |
| DNA shape                    | `analysis/analyses/dna_relax/figures/crystal_vs_ensemble_mgw.png`, `af3_vs_ensemble_mgwfl.png`, `engrailed_mgwfl_exemplar.png` |
| diversity                    | `analysis/figures/D1_diversity.png`                                                                                            |
| PWM exemplars                | `analysis/analyses/struct_pwm/figures/struct_pwm_{ets1,tbp,csl}.png`, `struct_pwm_context.png`                                 |

**To make** — all eight are table-driven, no new compute. Suggested home:
`analysis/analyses/nmr_rmsd/make_talk_figs.py` (N1–N7) and
`analysis/analyses/conf_vs_crystal_pwm/make_talk_figs.py` (N8); import `analysis/common/fig_common.py`
and the repo-root `palette.py`, per `analysis/docs/LAYOUT.md`.

| id  | slide   | figure                                                 | source table                                                             |
| --- | ------- | ------------------------------------------------------ | ------------------------------------------------------------------------ |
| N1  | S8, S21 | three-way fnat violins with 0.5 floor                  | `nmr_rmsd/fnat/*_fnat_{holo_states,apo_states,bioemu}.csv`               |
| N2  | S9      | pair roster table                                      | `nmr_rmsd/pairs.csv`, `fnat/fnat_three_way.csv`                          |
| N3  | S11     | core-construction schematic, annotated with core sizes | `nmr_rmsd/pair_core_summary_ca.csv`                                      |
| N4  | S16     | `%across` lollipop, shaded by residual ratio           | `pair_pca_summary_ca.csv`, `pair_projection_ca.csv`                      |
| N5  | S12     | within vs between separation with bootstrap CIs        | `pair_bundle_summary_ca.csv`                                             |
| N6  | S17     | ensemble/bundle residual ratio bars                    | `pair_pca_summary_ca.csv`                                                |
| N7  | S20     | per-pilot BioEmu vs deposited apo fnat                 | `nmr_rmsd/fnat_pilots/fnat_pilots_summary.csv`                           |
| N8  | S28     | paired SD slopes + cross-pilot delta                   | `conf_vs_crystal_pwm/per_tf_all_conditions.csv`, `crosspilot_deltas.csv` |

Plus two hand-drawn schematics (S7 three-way design, and the S30 in-flight table) and one PyMOL render
(S25 Sox2 2LE4/9QPF superposition).

---

## Part 5 — Numbers to re-derive before this is presentable

Three slides currently rest on numbers whose aggregation I could not reproduce. None of them changes a
direction; all of them are the kind of thing someone writes down.

| slide | quoted in v4                                                 | recomputed 2026-10-06                                           | action                                                                                                 |
| ----- | ------------------------------------------------------------ | --------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------ |
| S8    | holo 0.839 / apo 0.629 / BioEmu 0.396; pass 100 / 92 / 11.5% | pooled over 16 pairs: 0.844 / 0.636 / 0.393; pass 90 / 71 / 23% | pick one aggregation (pooled states vs per-pair medians vs informative-only) and state it on the slide |
| S28   | SD 0.235 → 0.180, 11/12, p = 0.027                           | per-pilot table: 0.247 → 0.202, 10/12, p = 0.034                | likely per-seed vs per-pilot unit; decide which                                                        |
| S27   | conformation SD 0.235 vs seed SD 0.16                        | not re-derived                                                  | recompute from `per_state_all_conditions.csv`                                                          |

Also: `%across`, residual ratios and the gate CIs in the tables above are quoted from the theme README,
which states that two of its tables predate the expansion from 5 to 16 pairs. The 14-pair `%across` table
is current; the 2×2 / delta table is not. Regenerate before putting either on a slide.

---

## Part 6 — Timing and cuts

| act                                    | slides  | min     |
| -------------------------------------- | ------- | ------- |
| 0 Where we left off                    | S1–S4   | 4       |
| 1 Do apo frames rebuild the interface? | S5–S8   | 7       |
| 2 How the pair experiment works        | S9–S14  | 12      |
| 3 Where the ensembles sit              | S15–S20 | 12      |
| 4 TFs and DNA                          | S21–S25 | 9       |
| 5 Augmentation, briskly                | S26–S29 | 6       |
| 6 Directions                           | S30–S32 | 10–12   |
| **total**                              | **32**  | **~60** |

**If running long, cut in this order:** S19 (trim sweep) → S24 (AF3 vs ensemble) → S27 (conformation vs
seed noise) → S14 (subspace self-validation, fold one line into S13) → S2 (project in one slide, fold
into S1).

**Never cut:** S8 (three-way), S11 (the shared core — it is what makes the rest honest), S16 (2–94%),
S17 (width signature), S25 (Sox2), S31 (directions).

**If the room gets interactive during Act 2 or 3, let it.** That is the outcome this talk is for. Protect
only S16, S17 and S31; everything else can be told in two sentences.

---

## Part 7 — Delivery

- Say "exploration" in the first minute and again before S31. The room should feel invited, not briefed.
- Rehearse S8, S11, S16 and S25 aloud. S11 is the one where reconstructing the logic live would cost you.
- Effect sizes before p-values in speech. `n` on every statistic slide.
- On every negative: state it once, say what it rules out, move on. Do not apologize and do not elaborate
  unless asked — the backup slides exist for that.
- When a number is soft, say "this one will move." You have three of them (Part 5) and naming them buys
  more credibility than hiding them costs.
- S31 stays up through the discussion. Take notes on it directly.

---

## Part 8 — What not to include

- The full five-beat bend-IQR withdrawal in the main line. One slide (S29), backup for the rest.
- The runx apo reclassification as a result — `pilot_mechanism.csv` is mid-regeneration.
- TBP as a decisive example of anything (label quarantined).
- The 2×2 / delta table from the NMR README — it predates the 16-pair expansion.
- The within-pilot +0.283 as an accuracy claim; it is in-sample and only interesting next to +0.025.
- A pilot-candidate shopping list. It competes with S31 for the same discussion time.
