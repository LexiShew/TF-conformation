# Conformational Ensembles for DNA-Binding Specificity
### A book about the TF-conformation project — working outline (v2)

> **Scope.** This project augments **DeepPBS** (protein structure → binding-specificity PWM)
> with **conformational ensembles** of monomeric TF–DNA complexes: per TF, sample a backbone
> ensemble with **BioEmu**, rebuild side chains with **HPacker**, dock each conformation onto the
> crystal DNA, filter by interface fidelity (fnat gate), minimize with OpenMM, featurize for
> DeepPBS, and train paired **baseline vs augmented** models. The intellectual payoff is a reframe:
> the *sign* of the augmentation effect becomes a **functional probe of recognition mechanism**
> (conformational selection vs induced-fit-on-DNA) across TF families.
>
> **Audience — future lab members / successors.** Written for a researcher inheriting or extending
> this pipeline. Assumes structural-biology and ML literacy (PWMs, DBDs, Kabsch, GBSA are used, not
> taught); invests instead in *how to run and extend every stage*, *why each design choice was made*,
> and *what was already tried so you don't repeat it*. Reproducibility and the fix spec are
> first-class, not appendix material.
>
> **Framing — mechanism-first.** The headline empirical result — augmentation's net effect is small
> and mostly negative on the general-130 benchmark — is deliberately staged as **setup, not verdict**.
> The conformational-selection vs induced-fit-on-DNA map (Ch. 17) is the payoff the fidelity and
> effect chapters build toward. The hypothesis is threaded from Ch. 4 as the book's backbone; every
> results chapter is read as a test of it. Where the evidence is thin (n=6, single-seed), the book
> says so and points at the control that would settle it (per-frame PWM labels, Part V).

---

## Part I — Motivation: why conformation should matter for specificity

**Ch. 1 — The specificity problem.**
- What a PWM encodes and what it flattens away — position-independence, no readout of shape or context.
- How protein structure encodes DNA-binding preference: base readout and shape readout at the interface.
- Where DeepPBS sits in the landscape: a structure → specificity model that takes the DNA as fixed.
- The gap this project attacks: a single crystal is one frozen pose of a molecule that lives as an ensemble.
- *Figures → (create) a schematic: crystal PWM ⟶ DeepPBS ⟶ predicted specificity, with the "fixed DNA" assumption called out.*

**Ch. 2 — The rigid-structure blind spot.**
- DeepPBS is constitutively blind to DNA deformation and to induced-fit protein rearrangement.
- Why that blindness, reframed, is a clean *instrument* rather than only a limitation — it isolates the protein-side question.
- BioEmu's sampled ensemble introduced as a **computational apo (free) state**; the crystal as **holo (bound)**.
- The apo/holo contrast as the conceptual lever the rest of the book pulls on.
- *Figures → (create) an apo/holo schematic: free-state ensemble (many faint poses) vs the single bound crystal, with the axes DeepPBS is blind to (DNA deformation, induced fit) annotated.*

**Ch. 3 — Why free-state ensembles, and why BioEmu rather than AlphaFold3.**
- The obvious alternative a successor will ask about first: why not just predict the TF–DNA complex directly with a folding model like AlphaFold3, and skip the sampling + docking machinery?
- **The complex-modeling gap:** current structure predictors don't reliably model protein–DNA complexes — which is why the pipeline docks a sampled protein onto crystal DNA and leans on liberal cutoffs + energy minimization instead of trusting a predicted complex.
- **The diversity gap (the decisive argument):** even as a *free-protein* predictor, AF3 collapses to essentially one pose. Across 2 seeds × 5 samples, AF3's median pairwise Cα-RMSD is ~0.1–0.9 Å for the monomeric pilots (runx 0.10, ets1 0.12, engrailed 0.13, tbp 0.15, egr1 0.17 … lef1 0.93), while **BioEmu spreads 1.5–6 Å — an order of magnitude more** (ets1 1.62, tbp 1.50, foxa 5.62, lef1 5.95).
- **Why that settles the method choice:** the whole premise is to populate the free-state conformational ensemble; a model that returns one high-confidence structure cannot do that. BioEmu emulates the equilibrium ensemble, which is exactly what the apo/holo framing (Ch. 2) needs.
- **The DUX4 caveat that proves the rule:** AF3's *only* large apparent diversity is the DUX4 dimer (median pairwise 13.7 Å), and that is assembly-placement variance between the two copies, not conformational sampling — reinforcing the monomer-only scope.
- Forward link: the raw diversity comparison and the D1 per-pilot figure live here as motivation; the metal-cage on/off result (a genuine result, not a motivation) stays in Part IV (Ch. 19).
- *Figures → **D1_diversity** (AF3 vs BioEmu pairwise-RMSD spread, per pilot); **dna_relax/af3_vs_ensemble_mgwfl** (the DNA-side version of the same argument — AF3 collapses minor-groove-width fluctuation toward zero while the physical ensembles spread wide, 12 pilots; reproduces the paper's Fig 3 flexibility claim); **pymol/ensembles_montage** and the per-TF **pymol/<tf>_af3** vs **<tf>_bioemu** pairs (visual: AF3 one pose, BioEmu a fan). (create) a single summary bar/strip of median pairwise Cα-RMSD, AF3 vs BioEmu, with the DUX4 outlier flagged as assembly variance.*

**Ch. 4 — The core idea and the hypothesis.**
- The augmentation idea: train DeepPBS on the conformational ensemble, not just the single crystal.
- The reframe — *does a free-state protein ensemble improve prediction?* — isolates whether free-state conformers are binding-relevant.
- The central hypothesis stated up front: augmentation should help under **conformational selection** and hurt under **induced-fit-on-DNA**.
- A reader's map: why every subsequent results chapter is a test of this one claim.
- *Figures → (create) the 2×2 conceptual schematic — protein reachability (low/high) × DNA deformation (low/high) — with the predicted augmentation sign in each quadrant; this is the hypothesis figure the mechanism chapter (Ch. 17) later fills with data.*

## Part II — Process: the pipeline, stage by stage

*For the successor: each chapter here ends with a **"how to run it / how to extend it"** box — the
actual command, the config knobs, the failure modes, and what a green result looks like — so this
part doubles as the runbook. The narrative explains* why *each stage is shaped the way it is; the
boxes tell you how to drive it.*

**Ch. 5 — Architecture.**
- The per-TF SLURM DAG: one pilot = `config/pilots/<tf>.sh`, submitted stage-by-stage with `afterok` dependencies via `run_pilot.sh`.
- A self-contained repo: vendored DeepPBS package and the 3DNA toolchain under `lib/`, so nothing is sourced from a separate deployment.
- The stage layout — each stage dir runs its co-located script; `wrappers/` holds one SBATCH wrapper per stage; `lib/common.sh` resolves paths and loads the pilot config.
- The endeavour cluster context: `rohs` partition / `rohs_102` account, three conda envs (`bioemu`, `deeppbs`, `hpacker`), and what lives on-disk vs. in-repo.
- *Runbook box:* `run_pilot.sh <tf> [start] [end]`; the config fields a pilot must set; how to read the DAG's job graph.
- *Figures → (create) the pipeline DAG diagram — the 7 stages + fnat gate as boxes with `afterok` edges, envs and on-disk outputs annotated per stage; this is the book's one orientation figure.*

**Ch. 6 — Stage 1: backbone ensembles (BioEmu) + side chains (HPacker).**
- BioEmu backbone sampling + HPacker full-atom side-chain reconstruction in one stage, for **every** protein chain of the structure.
- The per-chain output library (`output/stage1_bioemu/<PDB>_chain<X>_conformations/`) and its `topology.pdb` / `samples.xtc` / `samples_sidechain_rec.*` contents.
- Why the refactor happened: the old single-monomer `stage1_hpacker/` couldn't handle multi-chain complexes or do sampling and reconstruction together.
- Resumability: a chain whose `samples_sidechain_rec.xtc` exists is skipped, so re-submitting only fills gaps.
- *Runbook box:* single-structure vs whole-library submit; the `--all-chains` engine flag; expected output tree.
- *Figures → **rmsd_analysis/sampling_quality** (how broadly BioEmu samples per TF); **rmsd_analysis/sidechain_vs_backbone** (HPacker side-chain vs backbone contribution). (create) a small before/after showing HPacker rebuilding side chains on a BioEmu backbone-only frame.*

**Ch. 7 — Stage 2: interface-aligned re-dock.**
- Per-frame Kabsch dock of each conformation onto the crystal DNA, carrying DNA + structural metals in the protein's coordinate frame (**frame preservation**, INVARIANT 1).
- Why interface-Cα alignment is canonical: global all-Cα caps fnat at ~0.47 by averaging placement error; interface alignment raises it to ~0.72 (INVARIANT 2).
- The monomer guard and `BINDING_CHAIN` selector — docking the right chain's ensemble against the right `.cif` chain (the B2/B3 fixes).
- `per_domain` as a **diagnostic-only** mode (INVARIANT 3): valid for independent C2H2 subsites, dangerous for cooperative binders.
- fnat scored with `--use_model_dna` — model protein vs the DNA embedded in the same docked file.
- *Runbook box:* `--inspect-only` chain layout; `STAGE2_ALIGN_MODE`; what a sane interface + fnat readout looks like.
- *Figures → **S2_stage_progression_ets1** (crystal vs docked vs minimized, tight overlap for a rigid module — shared with Ch. 9). (create) the interface-vs-global alignment illustration: the same frame docked both ways, showing why global caps fnat ~0.47 and interface reaches ~0.72.*

**Ch. 8 — The fnat gate.**
- The single structural-quality filter: fraction-of-native-contacts vs the model's own DNA, dropping wrong-register states before they become training pairs.
- **Why it runs post-minimization (Stage 3), not at dock time:** minimization moves per-state fnat both ways — it can lift a near-miss over the floor or loosen a good pose.
- How it builds `${STAGE3_DIR}_pass/` as a symlink mirror of only the survivors; Stage 4 reads only that dir.
- Fail-loud design: an empty pass-list exits non-zero and the `afterok` edge halts the DAG — no training on an ungated set.
- The one knob, `FNAT_FLOOR=0.5`, overridable globally or per-pilot; why the floor was never varied.
- *Figures → **F1_fnat_distributions** (per-state fnat with the 0.5 floor, pass/fail colored); **F2_passrate_bars** (survival fraction per pilot, 100%→19%); **F3_fnat_vs_iRMSD** (the two criteria co-vary, ρ≈−0.84); **F4_interface_size** (fidelity is not interface size, n.s.).*

**Ch. 9 — Stage 3: energy minimization (the physics stage).**
This is the pipeline's one physics-based refinement step, and the chapter earns its length: the
docked frames from Stage 2 are geometrically plausible but energetically strained, and how we relax
that strain shapes every downstream feature and label. Because the successor may not have a molecular-
dynamics background, the chapter is broken into sub-chapters that each introduce one MD concept, say
why it matters *in general*, then show the specific choice this pipeline made and why. The engine is
`stage3_minimize/stage3_minimize.py` (OpenMM).

> **MD-concepts box (read first).** A *force field* is a formula for the potential energy of a set of
> atoms given their positions — springs for bonds and angles, periodic terms for torsions, and
> non-bonded electrostatics + van der Waals. *Energy minimization* rolls the structure downhill on
> that energy surface to the nearest low-energy geometry; it is **not** dynamics — no heating, no
> time evolution, no sampling of new conformations. (This pipeline sets up a Langevin integrator only
> because OpenMM requires one to build a simulation; it calls `minimizeEnergy`, never `step`, so
> temperature and friction never actually act.) The goal here is narrow and deliberate: take a
> clashy docked pose and settle it into the nearest physically reasonable one **without letting it
> wander away from the bound geometry we trust.**

**Ch. 9.1 — Why minimize a docked pose at all?**
- What Stage 2 hands over: a BioEmu conformer Kabsch-fit onto crystal DNA — right global placement, but with atomic clashes, stretched bonds, and strained contacts at the new interface.
- Why raw docked poses can't be featurized directly: DeepPBS reads interface geometry, and un-relaxed clashes would encode artifacts, not recognition chemistry.
- What minimization buys: it removes clashes and relieves local strain (the chapter's clash-count diagnostic drops from hundreds to near-zero, minimum heavy-atom protein–DNA distance opens up), producing a structure that is *physically* self-consistent.
- The one-sentence contract: relax the strain, preserve the pose — every design choice below serves that tension.
- *Figures → (create) the clash-count trajectory: heavy-atom protein–DNA clashes and min-distance across Phase 0 → ramp stages → final (the engine already prints these per state; plot one representative state).*

**Ch. 9.2 — The force field and the solvent model.**
- **AMBER ff14SB (`amber14-all.xml`)** — the energy function for protein and DNA; explain the anatomy (bonded springs + torsions + non-bonded terms) so the reader knows what is and isn't captured.
- **Implicit solvent, GB-Neck2 (`implicit/gbn2.xml`)** — instead of explicit water molecules, water is modeled as a continuous dielectric; why this is dramatically cheaper and standard for minimization.
- **Why it matters generally, and the caveat:** implicit solvent has no explicit counter-ions, so GBSA tends to **over-stabilize close ion pairs** — a charged protein–DNA interface can relax artificially compact (the Part III scientific caveat, stated here at its source).
- **`NoCutoff` non-bonded** — with implicit solvent there's no periodic box, so all pairwise interactions are computed exactly rather than truncated; hydrogen bond lengths are constrained (`HBonds`) in the ramped phases.
- Platform detail worth recording: CUDA mixed precision when a GPU is present, CPU fallback otherwise.

**Ch. 9.3 — Restraints: hold what you trust, relax what you don't.**
- **The core idea of a harmonic restraint:** a virtual spring, energy `½k(x−x₀)²`, that pulls an atom back toward a reference position `x₀`; the stiffness `k` (kcal/mol/Å²) sets how hard it resists moving.
- **What this pipeline pins:** protein backbone atoms (N, Cα, C) and — by default — DNA backbone atoms (P, C1′) are restrained to their docked positions at **k = 10 kcal/mol/Å²**; side chains and everything else are free.
- **Why restrain the backbone but not the side chains:** the backbone *is* the docked pose we trust from Stage 2; letting it float would discard the frame's identity. Side chains carry much of the recognition chemistry and genuinely need to repack, so they're left free to relax.
- **Why it matters generally:** without restraints, minimization of a strained complex can slide the whole protein off the DNA to a lower-energy but biologically meaningless pose — restraints are how you relax *locally* without losing *global* placement.

**Ch. 9.4 — The metal-coordination cage.**
- **The problem:** many DBDs (C2H2 zinc fingers, GATA, homeodomain-Zn) depend on a structural Zn²⁺/Mg²⁺; but `removeHeterogens` strips all non-standard residues, and AMBER ff14SB can't parameterize a bare metal ion anyway — so a naive minimization collapses the coordination shell into an unphysical apo-Zn fold.
- **The trick (no metal in the simulated system):** before stripping heterogens, detect each structural metal (ZN, MG, MN, FE, CA, CO, NI, CU) and record which side-chain heavy atoms coordinate it (ligand atoms SG/ND1/NE2/OD1/OD2/OE1/OE2/OG/OG1/OH within a **3.0 Å** cutoff; clusters of ≥2 ligands kept).
- **How the geometry is preserved:** add pairwise harmonic bonds between every pair of coordinating side-chain atoms, pinned at their original separations, at **k = 20 kcal/mol/Å²** — a "cage" that holds the coordination shell in place without ever simulating the ion.
- **Why it matters generally:** it's a reusable pattern for any cofactor a force field can't describe — restrain the *geometry the cofactor imposes* rather than the cofactor itself.
- **Diagnostics and the A/B:** the engine reports mean/max cage drift after minimization; `--ignore-metals` reproduces the pre-cage behavior for the metal-cage on/off comparison (Ch. 19 — the cage matters most for TBP, where minimization is otherwise least constrained).
- *Figures → (create) a zinc-finger coordination-shell render, cage-ON vs cage-OFF after minimization (intact tetrahedral Zn geometry vs collapsed apo-Zn), with the cage-drift number annotated.*

**Ch. 9.5 — The three-phase ramp (why not just minimize once).**
- **Phase 0 — hydrogens only:** freeze every heavy atom (set mass = 0), minimize just the newly added H's (≤2000 iters); fixes hydrogen placement before it can perturb the heavy-atom frame.
- **Phase 1 — van der Waals ramp:** the key move for clashy starts. Atoms that overlap feel enormous repulsive forces that can explode a minimization; so their vdW radii (σ) are scaled down and grown back over five stages (**σ = 0.1 → 0.3 → 0.5 → 0.7 → 1.0**, 500 steps each), letting clashes resolve gently as atoms "swell" to full size — with backbone and cage restraints active throughout.
- **Phase 2 — final minimization:** at full vdW, a long final relaxation (up to 10,000 iters) settles the structure into its minimum.
- **Why it matters generally:** softening-then-growing (a soft-core-like schedule) is the standard way to bring a badly clashing model to a stable state without numerical blow-ups — the alternative, one hard minimization, frequently diverges on docked poses.
- *Figures → (create) the ramp schedule as a plot: potential energy and clash count vs ramp stage (σ = 0.1 → 1.0), showing the gentle descent versus the divergence of a single hard minimization.*

**Ch. 9.6 — DNA treatment: frozen default vs relaxed variant.**
- **Frozen DNA (default):** DNA backbone shares the protein's k = 10 restraint, so the double helix is effectively held at its crystal pose while the protein relaxes around it — the pipeline's baseline, and the source of the "rigid-DNA ceiling" the mechanism story diagnoses.
- **The relaxation variant (`_dnarelax`):** split the DNA onto its own restraint force (`k_dna`) so it can be tuned independently — a **soft tether at k = 1.5 kcal/mol/Å²** lets the duplex deform toward the protein-preferred geometry while paying an elastic cost that stops it fraying; `k_dna = 0` releases it entirely.
- **Late release (`--dna-release-stage`):** hold the DNA at the stiff protein k for the early ramp stages so the protein settles first, then relax it to the soft `k_dna` at a chosen stage (5 in the `_dnarelax` config) and through the final minimization — co-relaxation without letting the DNA move before the protein is stable.
- **Why it matters generally:** the DNA restraint stiffness is the single knob that decides whether the pipeline augments "protein ensemble around a fixed binding pose" or "ensemble of genuine binding modes" — it's the physical embodiment of the book's central question.
- **Forward links:** the structural consequences (DNA bends toward the bound pose without fraying) are Ch. 18; the principled stiffness priors that should replace a flat k = 1.5 are Ch. 20.
- **Runbook box:** `--restraint-k` / `--dna-restraint-k` / `--dna-release-stage` / `--metal-cage-k` / `--ramp-stages`; the `_dnarelax` config pattern (`STAGE3_DNA_RESTRAINT_K=1.5`, `STAGE3_DNA_RELEASE_STAGE=5`); array-parallel Stage-3 recovery for states that fail an individual minimization.
- *Figures → **R1_ca_rmsd_stages** / **R3_minimization_delta** / **rmsd_analysis/minimization_motion** (minimization is a local relaxation, nudges toward crystal — shared with Ch. 15); the full DNA-shape consequence figures are Ch. 18. (create) a schematic of the release schedule: DNA restraint k held stiff for early stages, dropped to 1.5 at stage 5.*

**Ch. 10 — Stages 4–5: featurization and the augmented training set.**
- Stage 4: DeepPBS `process_co_crystal.py` turns each passing state into per-state `.npz` features (3DNA via vendored `lib/`).
- The helical-segment requirement and the self-complementary-DNA failure ("helix count problem 0") that Stage 4 exposes.
- Stage 5: building the augmented training fold + combined assembly + paired train configs.
- How augmented frames are labeled — **all frames of a complex inherit the one crystal PWM** (the confound that Part III and Part V both return to).
- *Figures → (create) a feature-tensor schematic (what a per-state `.npz` contains — protein-shape + DNA groove/shape channels DeepPBS reads); (create) an augmented-fold diagram: 1 crystal + N passing frames all pointing at one shared PWM label (visually foreshadowing the confound).*

**Ch. 11 — Stages 6–7: paired training and benchmark evaluation.**
- Baseline (crystal only) vs augmented (crystal + ensemble), always trained as a matched pair.
- Multi-seed design: single-seed vs 5 paired seeds (s1–s5), and which pilots got which.
- The general-130 held-out benchmark and the metrics reported (Pearson, Spearman, AUROC, IC-weighted PCC, MAE).
- Paired statistics (bootstrap / t-test) and the ΔPearson>0 = "augmentation helps" convention.
- *Runbook box:* `run_multiseed_pilot.sh`; where trained models and eval JSONs land; how to read `id_benchmark_<tf>.json`.
- *Figures → **figure_scripts/fig3_box_pearson** and **fig4_box_mae** (per-entry metric distributions over the 130-entry benchmark, three arms) as the "what the eval produces" reference; the interpretation figures (Δ by family etc.) come in Part IV.*

## Part III — Caveats: what the pipeline assumes and where it can mislead

**Ch. 12 — Scientific caveats.**
- **Rigid DNA carried from the crystal:** the protein's Kabsch transform is applied to the reference DNA and the backbone is then restrained — the interface is perturbed, not augmented (closer to "ensemble around a fixed binding pose" than "ensemble of binding modes").
- **Metal/cofactor stripping:** `removeHeterogens` can silently drop structural Zn²⁺/Mg²⁺; without the cage, zinc-finger / GATA / homeodomain-Zn folds minimize into an unphysical apo-Zn state.
- **GBSA implicit solvent** over-stabilizes close ion pairs at a charged protein–DNA interface, so the relaxed interface can be artificially compact.
- **The shared-PWM-label confound (confound #1):** every augmented frame shares one crystal PWM, so "does the ensemble help" is entangled with "does averaging N frames onto one label just regularize" — signal vs smoothing, currently inseparable.
- *Figures → (create) a "perturbed not augmented" schematic — protein frames varying around a frozen DNA backbone — making the rigid-DNA caveat visual; reuse the Ch. 10 shared-label diagram to anchor confound #1.*

**Ch. 13 — Statistical caveats.**
- No fixed random seed in generated configs — baseline vs augmented differ in augmentation *and* init *and* shuffle order, so with few seeds the augmentation effect isn't cleanly separable from noise.
- Best-state checkpoint selected on the wrong distribution: MAE on a general-TF validation set, not on the target TF — fair across arms, but the absolute numbers aren't "the best for TF X."
- Per-condition denominator drift: an npz that fails or returns NaN for one condition only shrinks that condition's n, so the Δ compares means over different subsets.
- Effect-size pairing is by regex; conditions not matching `baseline|augmented` are silently dropped from the summary with no warning.
- **The fig9 cautionary case:** the mixed-model "power-recovered" result (csl p=0.003, runx p<0.001) is **pseudoreplication** — it treats entry×seed rows as independent when augmentation is applied once per seed (n=5); correct seed-level tests give csl p=0.16, runx p=0.20. A worked example of how a plausible model overstates significance.
- *Figures → **figure_scripts/fig2_augmentation_delta** (the correct seed-paired ΔPearson with honest CIs) shown against **fig9_mixedmodel_effects** (the overconfident forest plot) as a side-by-side teaching contrast; **dna_relax/training_overfit_diagnostic** (loss curves, early-stopping histogram, late-training overfit, and generalization gap across the three arms — a caveats figure, not a results figure); (create) the per-seed delta dot-strip for csl/runx (2 of 5 negative) that exposes the coin-flip.*

**Ch. 14 — Engineering debt and the fix spec.**
- The root problem behind `PIPELINE_FIXES.md`: Stage 1 was refactored to a per-chain library, but `common.sh`, the Stage-2 wrapper, and pilot configs still encoded the old single-monomer layout — the two halves were disconnected.
- The invariants that must never regress: frame preservation, interface-canonical alignment, `per_domain` diagnostic-only, monomer-only scope.
- The blockers B1–B7 as a worked debugging narrative — from the Stage-1/Stage-2 path/naming break through the fnat-gate placement.
- Special case: self-complementary / symmetry-mate DNA (only one strand in the ASU) triggers "helix count problem 0"; the fix regenerates the biological-assembly duplex before featurization.
- Special case: the DUX4 dimer exclusion — the monomer guard and why a dimer has no place in a monomer benchmark.
- *Figures → (create) a small before/after of the self-complementary-DNA fix (single ASU strand → regenerated base-paired duplex); mostly a text/table chapter (B1–B7 as a table), few figures needed.*

**Ch. 14½ — Dead ends and things not to repeat (successor's chapter).**
- The two-docker foot-gun: the legacy PyMOL `cmd.align` global dock caps fnat ~0.47; `stage2_redock.py` is the sole production docker.
- Stage-2 fnat rejection removed: filtering at dock time shrank the denominator before minimization could rescue near-misses — the gate moved to post-minimization.
- The DUX4 dimer (0 fnat survivors) and `dux4_dnarelax` (still 0 survivors) — tandem/multi-chain pilots that never reached training.
- The sign-unstable `wt1` two-seed pilot in the predecessor tree: seed fold0 says augmentation hurts, seed1 says it helps — a caution against single-seed conclusions.
- Conditions that exist only as capability, never run: alignment baselines (`all` / `per_domain`), feature ablations (shape-only / groove / dnaseqInfo), folds ≠ 0, fnat floors ≠ 0.5.
- Format for each entry: what was tried, why it failed or stalled, and whether it is worth revisiting.
- *Figures → **pymol/dux4_bioemu** (the excluded dimer) as an illustrative "why it fails"; (create) the wt1 two-seed sign-flip dot plot (seed0 hurts, seed1 helps) as the single most memorable dead-end figure. Otherwise a narrative chapter.*

## Part IV — Results: building toward the mechanism map

*Read this part as one argument. Ch. 15 establishes that the ensembles are structurally faithful in
a rigidity-ordered way; Ch. 16 shows the net augmentation effect is small and mostly negative — the
apparent null that motivates the reframe; Ch. 17 resolves the null into a mechanism map and is the
climax; Ch. 18–19 stress-test that map.*

**Ch. 15 — Structural fidelity (the ensembles are trustworthy, unequally).**
- fnat, interface-RMSD, and Cα-RMSD across pilots, with the recurring **rigidity ordering**: ETS1, TBP (100% pass) → EGR1, engrailed (~91%) → FOXA (~75%) → LEF1 (19%).
- Fidelity tracks recognition-module **rigidity, not interface size**: TBP (40 interface residues) and LEF1 (39) are near-identical in size but opposite in fidelity.
- fnat and iRMSD co-vary tightly (Spearman ρ ≈ −0.84); pilots cluster along the curve by rigidity.
- Distortion is **localized** (one interface segment always moves more than the rest) and the per-residue signature is universal — rigid cores, floppy termini.
- Minimization is a local relaxation: global backbone barely moves, and per-state changes nudge marginally *toward* the crystal, never away. *Sets up the rigidity axis Ch. 17 splits on.*
- *Figures → **I1_iRMSD_distributions**, **I2_iRMSD_seg** (localized vs distributed distortion), **I4_interface_geometry**; **R1_ca_rmsd_stages**, **R2_per_residue_profiles** (rigid cores/floppy termini), **R3_minimization_delta**; **S1_bestworst_ets1** and **S1_bestworst_lef1** (best/worst frame renders — rigid hugs crystal, mobile flings off); **S2_stage_progression_ets1**. Optionally the `rmsd_analysis/` deep-dive set (state_trajectories, stagewise_progression) for a successor appendix.*

**Ch. 16 — The augmentation effect (the productive null).**
- Baseline vs augmented across 11 pilots: net effect small and **mostly negative** on the general-130 benchmark (most-negative egr1 −0.027, err −0.026, tbp −0.015; net-positive csl +0.027, nfat +0.015, ets1 +0.001).
- The **motif-level family map** and its correction of an earlier per-PDB mislabeling — ETS1–RUNX1 co-crystals were misfiled as Runt, inflating the ETS median.
- ETS (median ΔPearson +0.042, 60% of entries improve) and IRF (+0.039) are the **only** net-positive families; every other family is net-negative.
- The ETS1 gain is echoed by its ERG/FLI1 paralogs (biggest gains on low-baseline ERG entries), not confined to one structure.
- *Framed explicitly as setup:* a flat/negative mean is not "it doesn't work" — it is the signal that a hidden variable is splitting the families.
- *Figures → **P3_augeffect_by_family** (the key figure — ΔPearson by motif-level family, ETS/IRF the only net-positive); **P2_baseline_by_family** (baseline accuracy by family); **P1_family_table** (pilot→family map); **dna_relax/iface_mgwfl_vs_accuracy** (a second, DNA-side axis: more flexible interface DNA → augmentation hurts, ρ≈−0.55 — a companion to M1's protein-reachability axis, hinting the hidden variable is the DNA); **figure_scripts/fig1_three_arm_accuracy** and **fig6_within_family_transfer** (same- vs other-family transfer). Note the earlier per-PDB-mislabel version must NOT be used.*

**Ch. 17 — The mechanism story (the climax).**
- The apo/holo test: BioEmu ensemble = computational apo, crystal = holo, scored on two axes per pilot (free-state spread, reachability d_min) against qualitative DNA deformation.
- **No single protein-side axis predicts the sign** (all |ρ|<0.5, n.s., n=6) — the protein ensemble is not the problem.
- **Conformational selection** regime (ETS1, EGR1 helped): the free protein already samples binding-competent geometry and needs little DNA deformation.
- **Induced-fit-on-DNA** regime (TBP kink, LEF1 bend, hurt): free-protein frames carry no information about the DNA distortion DeepPBS holds fixed, so they add noise.
- **TBP as the decisive negative control:** best protein reachability (0.59 Å) yet hurt most (−0.26) → the missing variable is the DNA, not the protein sampling.
- Honest limits: n=6, single-seed, correlations not significant, DNA-deformation axis still qualitative — a hypothesis with a clean story and a decisive control, not a settled result.
- *Figures → **M1_apo_holo_mechanism** (the climax figure — reachability/spread vs DNA deformation, aug sign per pilot, TBP flagged). (create) the "hypothesis figure filled with data" — the Ch. 4 2×2 quadrant schematic overlaid with the six pilots at their measured positions, so motivation and result share one visual.*

**Ch. 18 — DNA relaxation.**
- The `_dnarelax` variant (soft tether k=1.5, released at ramp stage 5) lets protein and DNA co-relax into each frame.
- TBP DNA moves ~1.7× more under relaxation (median backbone RMSD 0.72 Å vs 0.41 Å frozen) — **without excess fraying**; the stiffness floor holds the duplex batch-wide.
- Displacement localizes to the central TATA bases (the kink region), and relaxed DNA **bends toward** the bound conformation (pyCurves: relaxed PP-bend lands closest to crystal; frozen underbends).
- The falsifiable prediction that ties it back to the mechanism map: flexible DNA should **recover exactly the families augmentation hurts today** (TBP, LEF1, homeodomain).
- Current status honestly stated: data for 7 pilots, trained for 2, evaluated for 1 (tbp_dnarelax, itself slightly worse than frozen tbp) — the prediction is not yet tested at scale.
- *Figures → **dna_relax/dna_shape_features** (the lead — a 6-panel violin overview across 7 TFs: global bend, axis shortening, minor/major groove width, BP inclination/tilt, crystal vs frozen vs relaxed); **dna_relax/tbp_dna_shape** (the 4-panel TBP ensemble: backbone RMSD, per-residue P displacement, P–P gap, Δbend); **dna_relax/crossfamily_bend** (Δ axis bend under relaxation per pilot, bootstrap CIs — most cross zero) and **crystal_convergence_bootstrap** (bend toward bound, with CIs); the per-TF `<tf>_dna_shape` set (11 pilots, missing hsf/irf) as a family appendix; interactive **pycurves_viz/*.html** viewers (crystal/frozen/relaxed) for the digital edition. (create) the falsifiable-prediction figure: predicted vs (eventual) actual aug-sign flip for the induced-fit families.*
- *Coverage caveat (applies to this whole chapter): the DNA-shape analysis has not been run on all 13 pilots — the overview spans 7, the per-TF panels 11 (no hsf/irf), the bend bootstrap 7, and the AF3/interface comparisons 12. Caption every DNA-relaxation figure as preliminary / partial-coverage: the story is directional, not yet complete across the full pilot set.*

**Ch. 19 — The metal-cage A/B (does the coordination cage change the result?).**
- The control this answers: the metal-coordination cage (Ch. 9.4) is a strong intervention on zinc-dependent folds — does turning it on/off move the downstream augmentation result, or is it invisible to DeepPBS?
- The A/B from the predecessor `old_results` tree (dux4, egr1, tbp; 5 seeds each; `--ignore-metals` = cage OFF).
- The cage matters most where minimization is least constrained: for TBP, cage-ON augmented (0.627) clearly beats cage-OFF (0.585); for dux4/egr1 the two are comparable — all still ≈ baseline.
- Why this lives only in the May-2026 predecessor tree and what that means for reproducing it.
- *Note:* the AF3-vs-BioEmu ensemble-diversity comparison, which used to sit here as a "reference comparison," moved to Ch. 3 — it is a motivation for the method (why BioEmu at all), not a result of it.
- *Figures → (create) the cage on/off bar chart — baseline vs aug(cage-ON) vs aug(cage-OFF) mean Pearson for dux4/egr1/tbp, TBP's 0.627-vs-0.585 gap highlighted (the numbers are in the predecessor `id_benchmark_<tf>_legacy_ab.json`; no plotted version exists yet).*

## Part V — Next directions

**Ch. 20 — Flexible DNA, done principledly.**
- The thesis from the mechanism result: the rigid-DNA assumption is the ceiling, not the protein sampling.
- **Tier 1** — relax DNA in the existing Stage-3 minimization (cheapest; the `_dnarelax` variant already built) — with the risk that unrestrained B-DNA frays under GBSA.
- **Tier 2** — a cgDNA+/cgNA+ sequence-dependent stiffness prior as a harmonic restraint: DNA free to deform toward the protein-preferred geometry but paying a sequence-dependent elastic cost (the principled hexABC connection; a per-hexamer stiffness lookup).
- **Tier 3** — Deep DNAshape as a learned prior: predicted sequence→shape used as the relaxation target so frames carry DNA geometry consistent with intrinsic preferences.
- Why not use hexABC trajectories directly: they characterize *free* B-DNA, whereas recognition is *protein-induced* deformation.
- *Figures → **dna_relax/stiffness_prior/stiffness_prior_demo** (the sequence-dependent stiffness concept); **dna_relax/perposition_minorgroove** and **mgw_fl_all12_panels** (per-position groove profiles the prior would target). (create) a three-tier schematic (existing minimization → cgDNA+ prior → Deep DNAshape prior) with cost/principledness axes.*

**Ch. 21 — Per-frame PWM labels: the decisive control.**
- The confound restated: all frames of a complex currently share one crystal PWM, so signal and smoothing are inseparable.
- The idea: give each frame its **own** label — a structure-derived per-frame PWM from that frame's protein–DNA contacts (even predicted per-frame labels suffice for the control).
- What it decides: if per-frame labels help where shared labels didn't → real conformation→specificity signal; if they wash out the gain → the benefit was regularization.
- Why it is decisive **for ETS specifically:** ETS has the lowest baseline of the helped families, so per-frame labels separate "genuinely conformationally selected" from "low baseline = most improvable by any regularizer."
- Sequencing: run the 5-seed control first (error bars), then per-frame PWM (the sharper question), then DNA relaxation.
- *Figures → (create) the decisive-control schematic — shared-label vs per-frame-label training, with the two predicted outcomes (signal survives / washes out) as branches; (create) a mock-up of the expected ETS result under each branch.*

**Ch. 22 — Scaling and completeness.**
- More seeds: put error bars on the per-family effects (the 5-seed re-run for 6 pilots, single-seed for the rest).
- More TFs and families: extend beyond the 11 pilots to fill out the family map.
- Coverage gaps to close: folds beyond 0, feature ablations (shape-only / groove / dnaseqInfo), and fnat-floor sensitivity — all present as capability, none exercised.
- DNA-relax eval for the 5 pilots that have data/configs but never reached evaluation.
- *Figures → (create) a coverage matrix (pilots × conditions, cells = run / data-only / capability-only) — the same grid as Appendix C but read as a to-do map of what remains.*

**Ch. 23 — The broader vision.**
- The reframe taken to its conclusion: from "does augmentation help" to "**augmentation sign as a mechanism map** of the TF universe."
- What a flexible-DNA, per-frame-labeled DeepPBS would enable — recognition-mechanism classification as a first-class output.
- The falsifiable through-line: if co-relaxing DNA flips TBP from negative toward neutral/positive, the mechanism hypothesis is confirmed and the method is fixed by the same change.
- Where this connects to the broader DNA-recognition and structural-ML programs.
- *Figures → (create) the closing "mechanism map of the TF universe" concept figure — families placed on the conformational-selection ↔ induced-fit axis, with the pilots as measured anchors and the rest of the universe as the frontier.*

## Appendices
- **A. Reproducibility.** The three conda envs (`bioemu` / `deeppbs` / `hpacker`) and their stage assignments; endeavour host/partition/account notes; the login-node OpenBLAS thread-cap gotcha; PyMOL edu-license install and headless ray-traced rendering.
- **B. Config reference.** A pilot config field-by-field (`PDB_ID`, `BINDING_CHAIN`, `PROTEIN_CHAIN`, `DNA_CHAINS`, `PWM_LABEL`, `FOLD`, `FNAT_FLOOR`, DNA-relax knobs); the worked "adding a new TF" walkthrough and its acceptance test.
- **C. Data & results inventory.** The 11-pilot × conditions matrix (std / 5-seed / dnarelax / AF3); the full per-pilot results table; the predecessor `old_results` tree; where every output lives on the cluster.
- **D. Figure atlas + production plan.** Master index of every figure, its home chapter, and its status. See the two tables below.
- **E. Glossary.** fnat, iRMSD, apo/holo, conformational selection, induced fit, PWM, DBD families, and the pipeline's stage vocabulary.

---

## Figure plan

> **How to read the per-chapter `Figures →` lines.** A **bold name** (e.g. **P3_augeffect_by_family**) is a figure that **already exists** on the cluster and just needs a caption + placement. A `(create)` item does **not exist yet** and must be produced. Existing files live under `analysis/figures/` (D/F/I/R/S/P/M series + `pymol/`), `analysis/dna_relax/figures/` (DNA-shape suite) and `pycurves_viz/` (interactive), `analysis/figure_scripts/` (three-arm benchmark fig1–9), and `rmsd_analysis/plots/` (the structural deep-dive set).

### D.1 — Existing figures → home chapter

| Figure (on disk) | Home | Role |
|---|---|---|
| `D1_diversity`, `pymol/ensembles_montage`, `pymol/<tf>_{af3,bioemu}` | **Ch. 3** | AF3-collapses / BioEmu-spreads motivation |
| `F1_fnat_distributions`, `F2_passrate_bars`, `F3_fnat_vs_iRMSD`, `F4_interface_size` | **Ch. 8** | fnat gate behaviour |
| `S2_stage_progression_ets1` | **Ch. 7** (+ reused Ch. 15) | crystal/docked/minimized overlap |
| `rmsd_analysis/sampling_quality`, `sidechain_vs_backbone` | **Ch. 6** | BioEmu sampling + HPacker |
| `R1_ca_rmsd_stages`, `R3_minimization_delta`, `rmsd_analysis/minimization_motion` | **Ch. 9.6** (+ reused Ch. 15) | minimization is a local relaxation |
| `figure_scripts/fig3_box_pearson`, `fig4_box_mae` | **Ch. 11** | what the eval produces |
| `I1_iRMSD_distributions`, `I2_iRMSD_seg`, `I4_interface_geometry`, `R2_per_residue_profiles`, `S1_bestworst_{ets1,lef1}` | **Ch. 15** | structural fidelity / rigidity ordering |
| `P3_augeffect_by_family` ⭐, `P2_baseline_by_family`, `P1_family_table`, `figure_scripts/fig1_three_arm_accuracy`, `fig6_within_family_transfer` | **Ch. 16** | augmentation effect by family |
| `M1_apo_holo_mechanism` ⭐ | **Ch. 17** | the mechanism climax |
| `dna_relax/tbp_dna_shape`, `crossfamily_bend`, `crystal_convergence_bootstrap`, `pycurves_ensemble_summary`, per-TF `<tf>_dna_shape`, `pycurves_viz/*.html` | **Ch. 18** | DNA relaxation consequences |
| `figure_scripts/fig2_augmentation_delta` vs `fig9_mixedmodel_effects` | **Ch. 13** | correct-vs-overconfident stats contrast |
| `dna_relax/stiffness_prior_demo`, `perposition_minorgroove`, `mgw_fl_all12_panels` | **Ch. 20** | stiffness-prior groundwork |
| `pymol/dux4_bioemu` | **Ch. 14½** | the excluded dimer |
| `rmsd_analysis/` deep-dive set (`state_trajectories`, `stagewise_progression`, `state_heatmap`, …) | **Appendix (successor)** | full structural drill-down |

### D.2 — Figures still to create

**Conceptual / schematic (no new computation — highest priority, the book has almost none of these):**
1. **Ch. 5** — the pipeline DAG diagram (7 stages + fnat gate, `afterok` edges, envs/outputs annotated). *The one orientation figure.*
2. **Ch. 4 → Ch. 17** — the 2×2 hypothesis quadrant (reachability × DNA deformation → predicted aug sign), reused in Ch. 17 overlaid with the six pilots' measured positions. *Ties motivation to result.*
3. **Ch. 1 / Ch. 2** — the DeepPBS "fixed-DNA" schematic and the apo/holo schematic.
4. **Ch. 10 / Ch. 12** — the shared-PWM-label diagram (N frames → one label), reused to anchor confound #1.
5. **Ch. 21** — the decisive-control schematic (shared vs per-frame label, two predicted branches).
6. **Ch. 20** — the three-tier flexible-DNA schematic (existing min → cgDNA+ → Deep DNAshape).
7. **Ch. 23** — the closing "mechanism map of the TF universe" concept figure.

**Data-backed (numbers already on disk; just need plotting):**
8. **Ch. 3** — median pairwise Cα-RMSD bar/strip, AF3 vs BioEmu, DUX4 flagged (from `ensemble_diversity.csv`).
9. **Ch. 19** — the metal-cage on/off bar chart (from the predecessor `id_benchmark_<tf>_legacy_ab.json`; no plot exists yet).
10. **Ch. 9.1 / 9.5** — clash-count + potential-energy trajectory across the ramp (the engine already prints these per state; parse a log and plot).
11. **Ch. 13** — the per-seed ΔPearson dot-strip for csl/runx (2 of 5 negative) exposing the coin-flip.
12. **Ch. 14½** — the wt1 two-seed sign-flip dot plot.
13. **Ch. 22** — the coverage matrix (pilots × conditions, run / data-only / capability-only).

**New computation / rendering required:**
14. **Ch. 6** — HPacker before/after (side-chain rebuild on a backbone-only BioEmu frame).
15. **Ch. 7** — interface-vs-global alignment illustration (same frame docked both ways, fnat 0.47 vs 0.72).
16. **Ch. 9.4** — zinc-finger coordination shell, cage-ON vs cage-OFF after minimization (render + cage-drift annotation).
17. **Ch. 10** — feature-tensor schematic (what a per-state `.npz` contains).
18. **Ch. 14** — self-complementary-DNA fix before/after (single ASU strand → regenerated duplex).
19. **Ch. 18** — the falsifiable-prediction figure (predicted vs eventual aug-sign flip for induced-fit families) — *pending the at-scale dnarelax eval.*

> **Statistical caution carried into the figure plan:** `fig9_mixedmodel_effects` is retained **only** as the negative example in Ch. 13 (its p-values are pseudoreplicated — seed is the experimental unit, n=5). Do not reuse it anywhere as a positive significance result. The honest counterpart is `fig2_augmentation_delta`.
