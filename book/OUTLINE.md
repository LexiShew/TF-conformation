# The Ensemble and the Groove
## A book-length account of the TF-conformation project — DETAILED OUTLINE

**Repo:** `/project2/rohs_102/shewchuk/TF-conformation` (endeavour HPC)
**Outline status:** v0.1 draft. Built from a full read of the repo: 20 markdown docs,
26 pilot configs, the F/I/R/S/P/M/D figure series, the fig1–9 benchmark series, the
`analysis/dna_relax/` structure suite, `rmsd_analysis/`, `af3/`, `deck/deck_spec.json`,
`docs/concerns.txt`, and the git history (2026-04-24 → 2026-07-30).
**Working title alternatives:** *The Ensemble and the Groove* · *A Negative Result, Carefully*
· *Conformational Augmentation: What the Protein Cannot Tell You About the DNA*

---

## 0 · The through-line (read this first)

Every chapter should serve one narrative spine. The spine is:

> A structure→specificity model (DeepPBS) sees one crystal pose. Proteins are ensembles.
> The obvious fix — train on an ensemble — **does not work**, and the *pattern of its
> failure* is more informative than a success would have been: the model is blind to DNA
> deformation, so free-protein conformers can only help in the narrow regime where
> recognition needs no DNA deformation and the free protein already reaches its bound pose.
> Fixing the DNA side (co-relaxation) does not rescue it either — and the attribution
> analysis says why: relaxed-DNA training moves the model *away* from real contact atoms.

Three acts:
1. **Act I — The bet.** Why an ensemble ought to help. Building the machine to test it.
2. **Act II — The verdict.** It doesn't help. Thirteen pilots, three arms, five seeds.
3. **Act III — The autopsy.** Where the information actually goes, and what would work.

**Narrative devices to thread throughout**
- **The two protagonists:** ETS1 (rigid, reaches its bound pose, the only consistently
  net-positive pilot) and TBP (the β-saddle that *also* reaches its pose perfectly yet is
  hurt worst, because it bends DNA ~80°). ETS1 and TBP are the book's control pair — they
  falsify every protein-side explanation and force the DNA-side one.
- **The recurring villain:** the shared PWM label. Introduced in Ch. 6, indicted in Ch. 11,
  put on trial in Ch. 17.
- **The recurring instrument:** the fnat gate. It is a quality filter in Act I, a *measuring
  device* in Act II (pass-rate = rigidity readout), and a suspect in Act III (it selects
  crystal-like frames, so what diversity survives?).
- **Honest-negative discipline.** The `REVIEW_figure_scripts.md` reviewer pass — which found
  a pseudoreplication error in the project's own fig9 — should be dramatized, not buried.
  It is the book's credibility.

---

# PART I — WHY CONFORMATION SHOULD MATTER

## Chapter 1 · The Reading Problem
*Goal: make a non-specialist feel why a protein finding one site in 3×10⁹ bp is hard.*

- 1.1 Transcription factors and the specificity problem; the combinatorics of the genome.
- 1.2 **Base readout vs shape readout.** Direct H-bond donors/acceptors in the major groove
  vs sequence-dependent minor-groove width, electrostatic potential, Arg-anchor recognition.
  Position DNA shape theory (Rohs lab lineage) as the intellectual home of the project.
- 1.3 The PWM as the field's unit of specificity — what it does and does not encode
  (independence assumption, no shape, no conformational state).
- 1.4 **What a co-crystal is and is not.** A crystal is one frame, selected by
  crystallizability, sampled at 100 K, in a lattice. The book's central tension in one image.
- 1.5 Chapter close: the question. *If specificity is read from a 3D interface, and the 3D
  interface fluctuates, what does a model trained on one frame per complex actually learn?*

**Existing figures:** none — this chapter is prose + textbook schematics.
**Figures to create:**
- **N1.1** Schematic: base readout vs shape readout on one duplex (major-groove H-bonds;
  minor-groove width + Arg insertion). Hand-drawn/BioRender style, cool palette.
- **N1.2** "One frame vs the ensemble" conceptual figure: crystal pose (grey) with a faded
  BioEmu ensemble behind it — a real render from `output/stage1_bioemu/`, used as the book's
  visual motif. *Source: `analysis/scripts/render_ensembles.py`.*

---

## Chapter 2 · DeepPBS — A Machine That Reads Structures
*Goal: the reader must understand the model well enough to know what it is blind to.*

- 2.1 The task: protein–DNA co-crystal → predicted PWM. Why this framing (structure-to-
  specificity) is different from sequence-to-specificity models.
- 2.2 Architecture in plain language: the bipartite protein→DNA contact graph, the
  symmetric-aware readout, the shape channels (`prot_shape` full model).
- 2.3 The 130-entry general benchmark; `id.txt` vs `valid0.txt` vs the general set; folds.
- 2.4 **The blindness that defines the book.** DeepPBS holds the DNA geometry fixed as an
  input feature. It cannot represent protein-induced DNA deformation as a *variable*.
  State this explicitly and early — it is the load-bearing premise of the mechanism
  hypothesis in Ch. 14.
- 2.5 Feature ablations that exist but were never run (shape-only / groove / dnaseqInfo
  pred_configs) — flagged here, revisited in the "unrun experiments" appendix.

**Existing figures:** deck slide "DeepPBS in one slide" (`deck/deck_spec.json`) — redraw.
**Figures to create:**
- **N2.1** DeepPBS data flow: structure → contact graph → PWM, with the DNA-geometry input
  boxed in a different colour and labelled "held fixed."
- **N2.2** Baseline accuracy landscape: per-entry Pearson over the 130-entry benchmark,
  sorted, showing the ~0.63 mean and the long low tail. *Data: `analysis/data/perentry_accuracy.csv`.*

---

## Chapter 3 · The Ensemble Hypothesis
*Goal: state the bet precisely enough to be falsifiable.*

- 3.1 Conformational selection vs induced fit — the two textbook mechanisms, and what each
  predicts for a free-state ensemble.
- 3.2 Generative structure ensembles arrive: BioEmu (and where AlphaFold3 sits — Ch. 12).
  What BioEmu samples: an **apo-like, free-state** backbone ensemble.
- 3.3 **The bet, stated formally.** If free-state conformers are binding-relevant, then
  training DeepPBS on {crystal + N docked conformers} should beat training on {crystal}.
- 3.4 The three ways the bet could fail, written down *before* the results:
  (a) conformers are not binding-competent; (b) conformers are binding-competent but
  carry no new information (label degeneracy); (c) the missing variable is not the protein.
  Chapter 14 will land on (c) — but the reader should have the menu in hand from the start.
- 3.5 Prior art: data augmentation in structural ML, MD-ensemble docking, ensemble docking's
  mixed record. Honest framing — this is not a novel idea, it is a carefully tested one.

**Figures to create:**
- **N3.1** Decision-tree schematic of the three failure modes with the eventual verdict
  greyed out (revealed in Ch. 14). A "roadmap" figure the reader can return to.

---

# PART II — BUILDING THE MACHINE

## Chapter 4 · Seven Stages
*Goal: the pipeline as an engineered artifact. This is the "nitty gritty" chapter.*

- 4.1 Design principle: **one pilot = one config file** (`config/pilots/<tf>.sh`), one SLURM
  DAG (`scripts/pipeline/run_pilot.sh`), self-rooted paths (`lib/common.sh`). Why the repo
  vendors DeepPBS and 3DNA under `lib/` — reproducibility as an architectural choice.
- 4.2 **Stage 1 — BioEmu + HPacker.** Backbone sampling per protein chain; side-chain
  reconstruction. Frame counts (~100/pilot). Why side chains must be rebuilt.
- 4.3 **Stage 2 — interface-aligned Kabsch re-dock.** Each frame superposed onto the crystal
  DNA using *interface* Cα atoms; DNA and structural metals carried across. The monomer
  guard. The `all` / `per_domain` alignment modes that exist but were never used in
  production — and the one experiment that tested this choice (`analysis/align_compare/`,
  Ch. 10).
- 4.4 **Stage 2g — the fnat gate.** Fraction of native contacts, floor 0.5. *Note the
  history:* the gate originally ran after Stage 2 and was moved to after Stage 3
  (commit 2026-06-25) — a fix worth narrating, because it changes what "pass" means.
- 4.5 **Stage 3 — OpenMM minimization** with a metal-coordination cage (sidechain-pair
  restraints, k=20) and a restrained DNA backbone (P, C1'). GBSA gbn2 implicit solvent.
  Phase-0 zero-mass freezing.
- 4.6 **Stage 4 — featurization** via the vendored DeepPBS `process_co_crystal.py` + 3DNA.
- 4.7 **Stage 5 — fold construction.** How an augmented training fold is assembled;
  where the shared PWM label enters the data.
- 4.8 **Stages 6–7 — paired training and benchmark eval.** Baseline and augmented arms are
  always trained as a pair; `id_benchmark_<tf>.json` is the unit of result.
- 4.9 Reproducibility furniture: `wrappers/`, `run_multiseed_pilot.sh`, the
  dependency-gated launcher pattern, `recover_failed_states.sh` array recovery.

**Existing figures:** deck slides "The augmentation pipeline (Stage 1 → 7)", "Stage detail —
sampling, docking, minimization", "Stage detail — the fnat gate, folds, and eval".
**Figures to create:**
- **N4.1** ⭐ **The pipeline plate** — full-page Stage 1→7 diagram with real per-stage
  artifacts (a BioEmu frame, a docked frame, a minimized frame, an `.npz` schematic),
  annotated with typical counts and the gate. This is the book's reference diagram; every
  later chapter can cite "Stage 3" against it.
- **N4.2** Sankey / attrition diagram: 100 sampled frames → docked → gate-passed → featurized,
  one band per pilot. *Data: `analysis/data/perstate_metrics.csv` + `RESULTS_INVENTORY.md`.*

---

## Chapter 5 · The Pilots
*Goal: introduce the cast. Each TF is a character with a biophysical personality.*

- 5.1 Selection criteria: monomeric DBD, crystal co-structure, JASPAR/HOCOMOCO PWM available,
  family spread. The PDB→PWM reference table (2026-06-24 commit).
- 5.2 The 13-pilot roster with structure, family, motif, and *what makes it interesting*:

  | Pilot | PDB | Family | Role in the story |
  |---|---|---|---|
  | ets1 | 1k79 | ETS | rigid winged-helix; the protagonist; 100% gate pass |
  | tbp | 1tgh | TBP β-saddle | the ~80° kink; the decisive negative control |
  | egr1 | 1aay | C2H2 ZnF | metal-dependent; the cage A/B case |
  | engrailed | 3hdd | Homeodomain | the textbook homeodomain |
  | foxa | 1vtn | Forkhead | floppy C-terminus |
  | lef1 | 2lef | HMG-box | ~110° bend, 19% gate pass — the failure exemplar |
  | csl | 3brg | CSL/RBPJ | later addition; sign-unstable |
  | err | 1lo1 | Nuclear receptor | 2 structural Zn |
  | nfat | 1a66 | Rel/NFAT | |
  | runx | 1hjc | Runt | duplicate duplex in the asymmetric unit |
  | hsf | 5d5u | HSF | added 2026-07; self-complementary duplex |
  | irf | 1if1 | IRF | palindromic ISRE; monomer assembly |
  | dux4 | 5z6z | Homeodomain (dimer) | **0 fnat survivors** — the scope boundary |
- 5.3 **Structure-preparation archaeology.** `structures/source_chains/`, chain selection,
  `BINDING_CHAIN` vs `PROTEIN_CHAIN` and the sequence-match guard that fails loudly.
  Per-pilot idiosyncrasies (runx's duplicated duplex; hsf's self-complementary strand;
  err's two Zn) — real work, and exactly the kind of detail this book should keep.
- 5.4 **dux4 as a boundary marker.** A homodimer run through a monomer pipeline: 94 docked
  frames, zero survive the gate, reachability d_min ≈ 10.5 Å. Not a failure to hide — a
  clean statement of the pipeline's domain of validity.
- 5.5 Family assignment done right: motif-level, not per-crystal. The ETS1–RUNX1 co-crystal
  mislabeling and why it inflated an early ETS result from +0.042 to +0.078.
  *(`analysis/data/family_annotation.csv`; the correction note in `analysis/README.md`.)*

**Existing figures:** `analysis/figures/P1_family_table.png`; `structures/pfam/clan_frequency.png`.
**Figures to create:**
- **N5.1** ⭐ **The pilot plate** — 13 small structural renders in a grid, uniform orientation
  and palette, each labelled family + gate pass-rate. Build from
  `analysis/scripts/render_ensembles.py` / `scripts/viz/_render_ensemble.py`.
- **N5.2** Table-figure: pilot × (PDB, family, motif, n frames, gate %, seeds, dnarelax?, AF3?).
  Merge `RESULTS_INVENTORY.md`, `CONDITIONS_INVENTORY.md` §B, and `results_inventory.csv`.

---

## Chapter 6 · The Label Problem
*Goal: plant the villain properly, in its own chapter, before any result.*

- 6.1 Every augmented frame of a complex inherits **one** crystal-derived PWM.
- 6.2 What that means statistically: augmentation by N is indistinguishable from
  up-weighting the original entry ×N, plus input noise. It is a regularizer with a
  structural flavour.
- 6.3 The two readings and why the experiment cannot currently separate them
  (from `mechanism_and_roadmap.md` §3).
- 6.4 What a per-frame label would have to be: structure-derived per-frame PWM from each
  frame's own contacts; or any label with frame-to-frame diversity.
- 6.5 Why this is *the* decisive control, and why it is expensive. Flag forward to Ch. 17.

**Figures to create:**
- **N6.1** ⭐ Schematic contrasting "N frames → 1 label" with "N frames → N labels", with the
  regularization interpretation drawn explicitly as duplicated training rows.

---

## Chapter 7 · Every Knob We Turned
*Goal: the honest conditions inventory — the chapter that makes the book trustworthy.*

- 7.1 The axes actually exercised: pilot, arm (baseline/augmented), DNA treatment
  (frozen k=10 vs relaxed k=1.5 released at ramp 5), metal cage on/off, seeds (1 vs 5),
  fold (0 only), alignment mode, gate floor (0.5 only), feature set.
- 7.2 The trial matrix — what was run, per pilot, per condition
  (`CONDITIONS_INVENTORY.md` §B).
- 7.3 **The predecessor tree.** `TF_conf_init_outputs/old_results/` (May 2026) holds two
  conditions that exist nowhere else: the metal-cage A/B (dux4/egr1/tbp × 5 seeds) and the
  original **wt1** pilot (two seeds, opposite signs: 0.631/0.556 vs 0.608/0.635).
  A short, valuable chapter on why you keep the old tree.
- 7.4 Capability-but-not-run: alignment baselines, feature ablations, other folds, other
  gate floors, dnarelax eval for 5 pilots. Each with an estimate of what it would cost.

**Existing figures:** none — this is a table chapter.
**Figures to create:**
- **N7.1** Condition-coverage heatmap: pilots (rows) × conditions (columns), cells =
  run / partial / capability-only. Directly from `CONDITIONS_INVENTORY.md` §B.
- **N7.2** Metal-cage A/B bar chart (dux4/egr1/tbp; baseline vs cage-ON vs cage-OFF) —
  TBP 0.627 vs 0.585 is the one place the cage clearly matters.

---

# PART III — WHAT THE ENSEMBLE LOOKS LIKE

## Chapter 8 · Rigidity Has a Number
*Goal: the first real result — the gate is a measuring instrument, not just a filter.*

- 8.1 fnat distributions per pilot, against the 0.5 floor. **F1.**
- 8.2 Pass-rate spans 100% (ets1, tbp, hsf) to 19% (lef1) to 0% (dux4). **F2.**
- 8.3 The two gate criteria co-vary (Spearman ρ = −0.84, n=561): fnat and iRMSD are one
  axis. **F3.**
- 8.4 ⭐ **Interface size does not predict fidelity** (ρ = −0.06, n.s.). TBP (40 interface
  residues, 100% pass) vs LEF1 (39 residues, 19% pass). Fidelity is about *module rigidity*,
  not contact count. **F4.**
- 8.5 Distortion is **localized**, not distributed: every point sits above y=x in segment
  max-vs-mean iRMSD. One segment always moves more than the rest. **I2.**
- 8.6 What a passing and a failing frame look like: LEF1's worst frame has its N-terminal
  helix flung off the DNA. **S1 (ets1, lef1).**

**Existing figures:** `F1_fnat_distributions.png`, `F2_passrate_bars.png`,
`F3_fnat_vs_iRMSD.png`, `F4_interface_size.png`, `I1_iRMSD_distributions.png`,
`I2_iRMSD_seg.png`, `I4_interface_geometry.png`, `S1_bestworst_ets1.png`,
`S1_bestworst_lef1.png`.
**Figures to create:**
- **N8.1** S1 renders for the remaining pilots (currently only ets1 + lef1) — at minimum
  tbp, egr1, foxa, hsf. `make_S.py` already takes a pilot list; `render_S.sbatch`.
- **N8.2** Gate pass-rate vs a *quantitative* rigidity proxy (ensemble RMSF or `spread` from
  `reachability.csv`) — turn the qualitative "rigidity" claim into a regression with a CI.

---

## Chapter 9 · Where the Motion Lives
*Goal: per-residue anatomy of the ensemble; the universal termini signature.*

- 9.1 Whole-protein Cα-RMSD from crystal, Stage 2 vs Stage 3. **R1.**
- 9.2 Per-residue profiles: rigid recognition core, floppy termini, universally.
  LEF1/FOXA C-terminal excursions 17–20 Å; TBP's β-saddle flat at ~1 Å. **R2.**
- 9.3 Minimization is a *local relaxation*: median ΔCα-RMSD −0.00 to −0.03 Å — it nudges
  frames marginally toward the crystal, never away. **R3, I3.**
- 9.4 The deeper `rmsd_analysis/` suite: state trajectories, ECDFs, rank stability,
  sidechain-vs-backbone decomposition, per-monomer and per-family summaries, state survival.
- 9.5 Interpretive point for the book: the ensemble's diversity is concentrated exactly where
  it is *least* relevant to the interface. That is a mechanism-level clue, delivered early.

**Existing figures:** `R1_ca_rmsd_stages.png`, `R2_per_residue_profiles.png`,
`R3_minimization_delta.png`, `I3_stage2_to_stage3.png`, `S2_stage_progression_ets1.png`,
plus `rmsd_analysis/plots/` (state_trajectories, stage_ecdf, rank_stability,
sidechain_vs_backbone, state_survival, summary_by_family, variant_agreement, …).
**Figures to create:**
- **N9.1** ⭐ Per-residue mobility mapped onto structure (B-factor-style colouring on one
  representative per family) — the same information as R2 but spatially legible.
- **N9.2** Interface-restricted vs whole-protein RMSD, side by side, all pilots — the
  `rmsd_analysis/` `_interface` variants already exist as separate plots; merge them into
  one comparison that makes the "core rigid / termini floppy" claim visually decisive.

---

## Chapter 10 · Two Diagnostic Detours
*Goal: short chapter; two self-contained methodological studies that would otherwise be lost.*

- 10.1 **Alignment mode matters enormously** (`analysis/align_compare/`). Interface-aligned
  vs globally-aligned docking: DNA displacement under interface alignment is ~0 by
  construction, whereas global alignment displaces the DNA by ~50 Å RMSD. Why interface
  alignment is the only defensible production choice, shown rather than asserted.
- 10.2 **The DNA-relax smoke test** — conditions A (frozen) / C (k=0) / D (tether + late
  release) on a single TBP state, the parameter-selection experiment that chose k=1.5 and
  release-at-ramp-5. Include the condition (B, k=10) whose output was not kept, as an
  honest note.

**Existing figures:** `F5_align_dna_displacement.png`;
`analysis/align_compare/dna_placement_by_mode.csv`, `dna_displacement_interface_vs_global.csv`.
**Figures to create:**
- **N10.1** Smoke-test panel: the single TBP state under A/C/D, DNA backbone overlay +
  RMSD-vs-restraint-strength curve. Reconstruct from `stage3_minimize/validation/`.

---

## Chapter 11 · How Diverse Is "Diverse"? — BioEmu vs AlphaFold3
*Goal: is the ensemble real, and is it the right kind of ensemble?*

- 11.1 Pairwise-RMSD diversity within each ensemble. BioEmu ets1: median pairwise 1.62 Å
  over 4,950 pairs; AF3 ets1 (2 seeds × 5 samples): 0.12 Å over 45 pairs — **AF3 samples are
  ~13× tighter.** AF3 is a high-confidence single-mode predictor, not an ensemble generator.
- 11.2 Distance to the crystal: AF3 lands much closer (egr1 0.37 Å vs BioEmu 2.15 Å) — AF3
  reproduces the *bound* pose; BioEmu explores the *free* one. This is the apo/holo
  distinction made concrete, and it sets up Ch. 14.
- 11.3 dux4 as the outlier in both (AF3 7.8 Å, BioEmu 15.9 Å) — the dimer problem again.
- 11.4 Why AF3 frames were never fed to training, and whether they should be (open question).
- 11.5 The AF3 input-construction work as a methods sidebar: DBD sequence extraction, duplex
  handling, structural metals as CCD ligands, `modelSeeds=[1,2]`, Apptainer on the GPU
  partition. Real, reusable, and currently undocumented outside the run scripts.

**Existing figures:** `D1_diversity.png`; `analysis/data/ensemble_diversity.csv`,
`ensemble_diversity_pairwise.csv`, `af3_rmsd_to_crystal.csv`,
`rmsd_to_crystal_af3_vs_bioemu.csv`; `analysis/scripts/montage_ensembles.py` output.
**Figures to create:**
- **N11.1** ⭐ Two-panel: (a) pairwise-RMSD violins, BioEmu vs AF3, all pilots; (b)
  RMSD-to-crystal, BioEmu vs AF3, same pilots. Currently split across CSVs and a
  single-purpose script; this deserves one canonical figure.
- **N11.2** Structural montage: crystal + AF3 fan + BioEmu fan for one rigid and one mobile
  pilot, same view — the visual statement of "AF3 = one pose, BioEmu = a cloud."

---

# PART IV — THE VERDICT

## Chapter 12 · It Does Not Work
*Goal: deliver the negative result cleanly, with error bars, without flinching.*

- 12.1 The unit of result: `id_benchmark_<tf>.json`, 130 entries, paired arms, per-seed.
- 12.2 Headline table — per-pilot ΔPearson and ΔMAE. **Two aggregations disagree and the
  book must say so:** the seed-pooled table in `RESULTS_INVENTORY.md` (csl +0.027,
  nfat +0.015, ets1 +0.001) vs the current `perentry_accuracy.csv` fold0 aggregation
  (csl −0.037, nfat −0.031, ets1 +0.001, hsf +0.003, runx +0.000). Reconciling *which
  checkpoints go into a mean* is itself a lesson worth a section.
- 12.3 The distributional truth: against entry-to-entry spread (Pearson ~0.63 ± wide), the
  three arms barely differ. **fig3, fig4.**
- 12.4 The seed-paired effect — the only cross-treatment quantity that is comparable, since
  the frozen and relaxed pipelines each retrained their own baseline. **fig2, fig5.**
- 12.5 Per-family means with CIs. **fig1.**
- 12.6 What "mostly slightly negative" means honestly: the effect sizes are ~0.01–0.04 Pearson
  against a benchmark s.d. that dwarfs them.

**Existing figures:** `analysis/figure_scripts/fig1_three_arm_accuracy.png`,
`fig2_augmentation_delta.png`, `fig3_box_pearson.png`, `fig4_box_mae.png`,
`fig5_mae_delta.png`; `analysis/figures/P2_baseline_by_family.png`.
**Figures to create:**
- **N12.1** ⭐ **The one-page verdict figure**: forest plot of seed-paired ΔPearson per pilot
  (frozen and relaxed), with the zero line, ordered by effect. If the book has one figure a
  reader remembers, this is it.
- **N12.2** Reconciliation figure for §12.2: the same pilots under three aggregation choices
  (fold0 only / seed-mean / pooled), showing how the sign of csl and nfat flips. An honesty
  figure, and a genuinely useful methods contribution.

---

## Chapter 13 · Family Structure in the Failure
*Goal: the one place a positive signal survives — presented as a lead, not a result.*

- 13.1 Motif-level family assignment (again — it changes the answer).
- 13.2 **P3:** ETS (n=10, median ΔPearson +0.042, 60% of entries improve) and IRF (n=4,
  +0.039) are the only net-positive families; every other family is net-negative.
- 13.3 The ETS1 gain is echoed by its ERG/FLI1 paralogs (largest gains on low-baseline ERG
  entries, 0.35→0.51) — so it is not one structure's fluke.
- 13.4 Within-family vs cross-family transfer. **fig6, fig7, fig8.**
- 13.5 ⚠️ **The statistics chapter-within-a-chapter.** `REVIEW_figure_scripts.md` found that
  fig9's mixed model (crossed random *intercepts*) treats entry×seed rows as independent
  replicates and shrinks the arm SE ~2×, producing csl p=0.003 / runx p<0.001. The correct
  seed-level analyses (paired t, n=5: csl p=0.16, runx p=0.20; random-slope mixed model:
  p=0.076 / 0.128) contradict it. Per-seed deltas are a coin flip (csl 2 of 5 negative).
  **fig9 must be refit or dropped.** Dramatize this — an internal reviewer catching the
  project's own most exciting p-value is the book's ethical centre.
- 13.6 Also correct the "rigid families" mislabel: fig2 and fig5 name *different* rigid sets,
  and both include TBP — the panel's strongest DNA bender.

**Existing figures:** `P3_augeffect_by_family.png`, `fig6_within_family_transfer.png`,
`fig7_abs_pearson_groups.png`, `fig8_abs_mae_groups.png`,
`fig9_mixedmodel_effects.png` (**to be refit**),
`analysis/dna_relax/figures/samefamily_augmentation_effect.png`,
`crossfamily_augmentation_effect.png`.
**Figures to create:**
- **N13.1** ⭐ Refit fig9: random slope of arm across seeds, or the plain n=5 seed-level
  paired test, side by side with the original — the "before/after" of a statistical
  correction, shown honestly.
- **N13.2** fig6 redrawn with symmetric jitter on both same-family and other-family groups
  (the reviewer's asymmetry complaint), plus an explicit per-pilot same-vs-other test.

---

## Chapter 14 · Apo, Holo, and the Missing Variable
*Goal: the intellectual climax. The mechanism hypothesis.*

- 14.1 Reframe: stop asking "does augmentation work" and start asking "what does the *sign*
  of the augmentation effect tell us about recognition mechanism?" DeepPBS's blindness to
  DNA deformation makes it a clean instrument.
- 14.2 Two protein-side axes, computed per pilot: free-state **spread** (ensemble diversity)
  and **reachability** `d_min` (how close the free ensemble gets to the bound pose).
  *Data: `analysis/data/reachability.csv`, `mechanism_apo_holo.csv`.*
- 14.3 **No single protein-side axis predicts the sign** (|ρ|<0.5, n.s.).
- 14.4 ⭐ **The TBP argument.** TBP has the *best* reachability in the panel (d_min 0.587 Å —
  its rigid β-saddle reaches the bound backbone essentially perfectly) and is hurt among the
  worst. The protein ensemble is not the problem. The missing variable is the DNA.
- 14.5 The two-corner reading (M1 panel B): augmentation helps only where the protein reaches
  its bound pose **and** recognition needs little DNA deformation (ETS1, EGR1, IRF, HSF) —
  the **conformational-selection** corner. It hurts under **induced fit on the DNA**
  (TBP's ~80° kink, LEF1's ~110° bend), because free-protein conformers carry no information
  about a distortion the model holds fixed.
- 14.6 The falsifiable prediction that Part V tests: *a flexible-DNA pipeline should rescue
  exactly the families augmentation hurts today.*
- 14.7 Caveats stated in full: n=13 with only some multi-seed; the DNA-deformation axis is
  curated from literature, not computed; correlations not significant.

**Existing figures:** `M1_apo_holo_mechanism.png`.
**Figures to create:**
- **N14.1** ⭐⭐ **The mechanism quadrant** — reachability (x) × DNA-deformation (y), points =
  pilots, colour/size = augmentation effect sign and magnitude, with the
  conformational-selection corner shaded. M1 has the ingredients; this should be the
  book's cover-candidate figure.
- **N14.2** Replace the curated `dna_deform` axis with the **measured** one now available
  (pyCurves bend from `analysis/dna_relax/data/pycurves_all_summary.csv`) — turning the
  hypothesis's weakest link into a computed quantity. *This is the single highest-value new
  analysis in the outline.*
- **N14.3** ETS1 vs TBP as a two-panel structural argument: same reachability quality,
  opposite DNA deformation, opposite augmentation sign.

---

# PART V — FIXING THE DNA SIDE

## Chapter 15 · Letting the DNA Move
*Goal: the direct test of Ch. 14's prediction — and its outcome.*

- 15.1 The Stage-3 modification: `--dna-restraint-k 1.5` with `--dna-release-stage 5`
  (soft tether, released late) vs the frozen default k=10. Configs for all 13 pilots.
- 15.2 **Does the DNA actually move?** Yes, ~1.7×: median DNA-backbone RMSD vs docked
  0.72 Å relaxed vs 0.41 Å frozen (TBP, n≈95 per arm).
- 15.3 **Does it fray?** No. Max adjacent P–P gap tails match between conditions; no state
  reaches the ≥9 Å unwind regime. The k=1.5 floor holds the duplex batch-wide.
- 15.4 **Does it move in the right direction?** Yes, for TBP: pyCurves overall bend
  crystal 79.7° / frozen 81.6° / relaxed 86.7° (UU); PP crystal 72.2 / frozen 68.0 /
  relaxed 73.5 — relaxed lands closest to crystal while frozen *under*bends. Displacement
  localizes to the central TATA bases — the kink region — not the ends.
- 15.5 Minor-groove width: all three conditions retain TBP's widened groove (~12 Å vs 5.7 Å
  canonical B-DNA). Per-position minor-groove profiles along the helix.
- 15.6 Cross-family: bend deltas with bootstrap CIs; crystal-convergence scores
  (does relaxation move each feature *toward* the crystal value?).
- 15.7 Minor-groove-width fluctuation (MGW-FL) as an ensemble observable, and the
  AF3 comparison: relaxed ensembles show ~1.6–2× the AF3 MGW-FL.
- 15.8 **But the accuracy did not improve.** tbp_dnarelax ΔPearson −0.034, *worse* than
  standard tbp (−0.015). The prediction of Ch. 14.6 fails at the accuracy level even though
  every structural intermediate behaves as predicted. State this plainly.

**Existing figures:** `analysis/dna_relax/figures/` — `tbp_dna_shape.png` (+ 11 per-pilot
`<tf>_dna_shape.png`), `dna_shape_features.png`, `dna_shape_features_fnatpass.png`,
`crossfamily_bend.png`, `bend_delta_bootstrap` data, `crystal_convergence_bootstrap.png`,
`perposition_minorgroove.png`, `pycurves_ensemble_summary.png`, `mgw_fl_all12_panels.png`,
`mgw_style_engrailed_demo.png`, `af3_vs_ensemble_mgwfl.png`, `iface_mgwfl_vs_accuracy.png`.
**Figures to create:**
- **N15.1** ⭐ The prediction-vs-outcome figure: structural convergence toward crystal (x)
  vs ΔPearson (y), one point per pilot. The clean statement that structural improvement did
  **not** buy accuracy — the chapter's whole argument in one panel.
- **N15.2** Ensemble-level pyCurves (bend/groove *distributions* over all states, not just
  state 2) — flagged as "the natural next step" in `analysis/dna_relax/README.md` and
  partially delivered in `pycurves_ensemble_summary.csv`; finish and standardize it.

---

## Chapter 16 · Where the Model Actually Looks
*Goal: attribution — the autopsy that explains why relaxation didn't help.*

- 16.1 Occlusion-based attribution for DeepPBS: mask a protein atom, measure the MAE shift
  in the predicted PWM. `interpret_tfconf_all.py`.
- 16.2 ETS1 result (99 benchmark complexes, three arms):
  - frozen-DNA augmentation: mean importance shift **+0.000489**, p=0.266, d=+0.11 —
    a marginal, non-significant nudge *toward* contact atoms.
  - relaxed-DNA augmentation: shift **−0.000198**, p=0.026, d=−0.23 — a significant move
    *away* from real contact atoms.
- 16.3 The interpretation, and its limits: relaxed frames may introduce interface artifacts,
  or the model learns to discount contacts that are no longer geometrically consistent.
  Effect sizes are small; this is a directional claim, not a strong one.
- 16.4 ⚠️ **The unfinished experiment.** All 12 pilots were submitted on 2026-07-30
  (jobs 5303457–5303468); `output/interpret_results_all/` contains 10 pilot directories but
  **`all_pilots_importance_summary.csv` does not exist** — hsf and irf are missing and the
  cross-pilot comparison was never compiled. Finishing this is a prerequisite for the
  chapter's central claim to generalize beyond ETS1.

**Existing figures:** `analysis/importance_shift_distributions.png`,
`importance_boxplot_comparison.png`, `importance_scatter_comparison.png`;
`importance_stats.json`, `importance_comparison_table.csv`.
**Figures to create:**
- **N16.1** ⭐ Cross-pilot attribution summary — mean shift ± CI per pilot × arm, once
  `compare_importance_all.py` has been run to completion. *Blocked on finishing §16.4.*
- **N16.2** Attribution mapped onto structure: contact atoms coloured by importance shift
  (baseline → frozen → relaxed) for ETS1. Makes an abstract number physical.

---

## Chapter 17 · Overfitting, Regularization, and the Label Again
*Goal: close the loop on Chapter 6's villain with the training-dynamics evidence.*

- 17.1 The training/validation diagnostic across **177–178 runs**: loss curves, best epoch,
  validation minimum, post-minimum rise, train–val gap, for baseline / frozen-aug /
  relaxed-aug regimes. *Data: `analysis/dna_relax/data/training_overfit_metrics.csv`;
  figure `training_overfit_diagnostic.png`.*
- 17.2 Does augmentation regularize? If it did, the augmented arms should show a smaller
  train–val gap and a later validation minimum. Test this explicitly — the data are on disk
  and the analysis is one script away.
- 17.3 The checkpoint-selection problem (`docs/concerns.txt` #6): `best_state_metric: mae`
  is computed on `valid0.txt` (general TFs) while the TF of interest lives in `id.txt`.
  Fair across arms, but it means absolute numbers are not "the best the augmented model can
  do for TF X."
- 17.4 The missing seed control (`concerns.txt` #5): configs do not set `no_random`, so
  baseline and augmented differ in augmentation *and* init *and* shuffle order. The
  five-seed re-run is the mitigation; the fix is one config flag.
- 17.5 Verdict on the regularization hypothesis, and why only per-frame labels can settle it.

**Existing figures:** `training_overfit_diagnostic.png`.
**Figures to create:**
- **N17.1** ⭐ Train–val gap and validation-minimum epoch, baseline vs augmented, paired by
  seed — the direct regularization test. Data already in `training_overfit_metrics.csv`.

---

# PART VI — WHAT WOULD ACTUALLY WORK

## Chapter 18 · Per-Frame Labels
- 18.1 The design: structure-derived per-frame PWM from each frame's own protein–DNA contacts.
- 18.2 What it decides — the two outcomes and their opposite conclusions (both publishable).
- 18.3 Why it is decisive *for ETS1 specifically*: ETS has the lowest baseline among helped
  families, so "low baseline = most improvable by any regularizer" is the live alternative.
- 18.4 Implementation sketch: a new labeling stage between 4 and 5; cost estimate; the
  validation you would need (does the per-frame label reproduce the crystal PWM for the
  crystal frame?).

## Chapter 19 · A Sequence-Dependent DNA Prior
- 19.1 Why uniform k=1.5 is wrong: DNA stiffness is strongly sequence-dependent — TA steps
  are hinges, AT/GC steps are rigid.
- 19.2 The hexABC table: per-step elastic constants from **380 MD sequences**, summed
  stiffness spanning TA/TA 10.0 → AT/AT 18.2 (~1.8× range), TA floppiest — matching TATA-box
  hinge biology. *(`analysis/dna_relax/stiffness_prior/`.)*
- 19.3 `seq_to_kdna()`: mean-normalized per-step k. Worked example — TBP's `CGTATATATACG`:
  TA steps k≈1.09, AT/GC flanks k≈1.98.
- 19.4 The four-step Stage-3 integration plan and the `STAGE3_DNA_STIFFNESS_PRIOR=hexABC`
  flag (default off → byte-identical to current behaviour). **Prototype only; not yet wired.**
- 19.5 Tier 3: learned deformability (Deep DNAshape) as the restraint centre — the in-house
  route, and how it connects to the lab's DNA-shape lineage.

## Chapter 20 · What a Model That Could Learn This Would Look Like
*Goal: the forward-looking chapter. Speculative but grounded.*
- 20.1 The requirement: DNA geometry as a *variable*, not a fixed input feature.
- 20.2 Joint protein–DNA ensembles rather than protein ensembles on a frozen duplex.
- 20.3 Reframing the target: from a PWM to a conformation-conditioned specificity landscape.
- 20.4 What this project's negative result contributes to that design — specifically, the
  quadrant of Ch. 14 as a *specification* for which systems any such model must handle.

## Chapter 21 · How to Be Wrong Well
*Goal: the methods-culture chapter; short, opinionated, and the reason a reader recommends
the book to a student.*
- 21.1 Paired arms, always. Seeds as the experimental unit, not entries.
- 21.2 The internal reviewer pass — and what it caught (fig9 pseudoreplication; the "rigid
  families" mislabel; fig6's asymmetric jitter; over-stated titles on honest error bars).
- 21.3 Keeping the predecessor tree: `old_results/` held the only metal-cage A/B and the only
  wt1 pilot.
- 21.4 Writing down the confounds *before* the result (`docs/concerns.txt`, dated).
- 21.5 Pilot-agnostic figure code: `fig_common.discover_pilots()` — nothing hardcodes the
  pilot list, so a new TF enters every figure automatically. Why this matters more than it
  sounds.
- 21.6 Environment discipline as science: which conda env, which node, why pyCurves cannot
  run on the login node, why PyMOL needs `ray=1` under sbatch.

---

# APPENDICES

- **A · Reproduction guide.** Every figure in the book → the script that makes it → the env
  it needs → the data it reads. Extend `analysis/figscripts/README.md` and
  `analysis/figure_scripts/README.md` into one table.
- **B · Full results tables.** Per-pilot, per-seed, per-arm, per-metric — the machine-readable
  companion (`perentry_accuracy.csv`, `perseed_summary.csv`, `perseed_perentry.csv`,
  `results_inventory.csv`).
- **C · Per-pilot dossiers.** One page per TF: structure render, family, motif, gate curve,
  DNA-shape panel, ΔPearson, and the one-line "what this pilot taught us."
- **D · The unrun experiments.** Alignment baselines, feature ablations, other folds, other
  gate floors, dnarelax eval for 11 pilots, per-frame labels, stiffness prior — each with
  a cost estimate and the question it would answer.
- **E · Known issues and open concerns.** `docs/concerns.txt` verbatim, annotated with what
  has since been fixed (metal cage added; gate moved after Stage 3) and what has not
  (no fixed seed; loose Cα count check; zero-mass freezing brittleness).
- **F · Software stack.** BioEmu, HPacker, OpenMM, 3DNA/DSSR, pyCurves, DeepPBS, AlphaFold3;
  versions, envs, and the vendoring rationale.
- **G · Timeline.** 2026-04-24 first commit → 2026-07-30, annotated with what each phase
  learned (from the git log — a genuinely readable project history).

---

## Cross-cutting production notes

**Figure inventory as it stands**
| Series | Location | Count | Status |
|---|---|---|---|
| F / I / R / S / P / M / D | `analysis/figures/` | 20 | pilot-agnostic, 13-pilot coverage |
| fig1–fig9 (benchmark) | `analysis/figure_scripts/` | 9 | fig9 needs refit |
| DNA-relax structure suite | `analysis/dna_relax/figures/` | 24 | most complete series |
| RMSD deep-dive | `rmsd_analysis/plots/` | ~28 + 100 per-PDB | largely unused in narrative |
| Attribution | `analysis/*.png` | 3 | ETS1 only; all-pilot run unfinished |
| Deck | `deck/TF_conformation_deck.pptx` | 29 slides | narrative skeleton already exists |

**Highest-value new work identified by this survey** (ranked)
1. **N14.2** — replace the curated DNA-deformation axis with measured pyCurves bend.
   Turns the mechanism hypothesis's weakest link into data. *Everything is on disk.*
2. **§16.4 / N16.1** — finish the all-pilot attribution run and compile
   `all_pilots_importance_summary.csv`. The Ch. 16 claim currently rests on one pilot.
3. **N13.1** — refit fig9 correctly. Required before any of it is publishable.
4. **N12.2** — reconcile the three aggregation choices; the sign of csl/nfat depends on it.
5. **N17.1** — the direct regularization test from `training_overfit_metrics.csv`.
6. **N15.1** — structural convergence vs accuracy gain, the Part-V punchline.
7. **N8.1** — S1 renders for the remaining pilots (currently ets1 + lef1 only).

**Style guide for the book's figures**
Cool pastels throughout — light blues, greens, purples, pinks; grey = baseline, teal =
frozen-DNA augmented, green = relaxed-DNA augmented, as already fixed in
`analysis/figure_scripts/_common.py` and `palette.py` (one hue per entity, repo-wide).
Note: `analysis/figscripts/README.md` references a `COLOR_CONSISTENCY.md` at the repo root
that **does not exist** — it should be written as part of the book's production setup.
