# Figure plan — Chapters 10–11 and 15–18

> **Status of this plan.** The existing-figure inventory below was **verified against the cluster** on
> 2026-08-04: all 43 analysis PNGs under `analysis/figures/`, `analysis/mechanism/figures/`,
> `analysis/figure_scripts/`, `analysis/dna_relax/figures/`, and `rmsd_analysis/plots/` were listed, and
> the six load-bearing figures (M1, M2, M3, M4, P3, fig3) were opened and read so the captions and the
> keep/swap decisions reflect what the figures actually show — not what the outline named. A figure marked
> **[exists ✓]** is confirmed on disk; **[create]** does not exist and must be produced.
>
> **How to read the columns.** Each chapter lists: figures to **reuse** (exist, just caption + place),
> figures to **create**, and any **caption caveat** where the existing figure's obvious caption would
> overstate what the data support.

---

## The lopsided summary (read this first)

| chapter | exists on disk | to create | state |
|---|---|---|---|
| 10 — featurization / augmented set | **0** | 3 (all schematic) | **unillustrated — highest build need** |
| 11 — training / evaluation | 2 (`fig3`, `fig4`) | 3 (1 data, 2 schematic) | thin |
| 15 — structural fidelity | 8+ (I/R/S series) | 0 | **fully covered** |
| 16 — augmentation effect | 5 + `M2` | 1–2 (data) | well covered; centerpiece is `M2` |
| 17 — mechanism (revised) | 3 (`M1`,`M3`,`M4`) | 1–2 (data + schematic) | covered, but **feature `M4`, not `M1`** |
| 18 — DNA relaxation | 7+ (dna_relax suite) | 1 (data) | well covered; one caption caveat |

**Two swaps that matter more than any new figure**, both using figures that already exist:
1. **Ch. 17: make `M4_dnaflex_vs_effect` the climax, not `M1_apo_holo_mechanism`.** `M1` plots the
   *protein-side* axes and hard-codes the "TBP best reach, hurt most" reading — on its own it tells the
   *static* story the data overturned. `M4` is the actual result (dynamic bend-IQR ρ = −0.76 relaxed vs
   flat static control ρ = +0.04). Caption `M1` as *the hypothesis*; let `M4` carry *the result*.
2. **Ch. 16: make `M2_where_the_signal_lives` the centerpiece.** Its four panels already *are* the chapter's
   argument (own vs other family, seed-level contrast, headroom control, dynamic-vs-static axis).

---

## Chapter 10 — Stages 4–5: featurization & the augmented set

**Existing: none.** No figure anywhere in the analysis tree covers featurization or set assembly. This chapter
is unillustrated in the current material — the top build priority of the six.

**To create (all conceptual, no computation):**
1. **Feature-tensor schematic** — what one `.npz` holds: the protein graph (charge, radius, SASA, Atchley
   factors, circular variance) + the DNA point cloud (`X_dna_point`) + the 14-column shape matrix (6 intra-bp,
   6 inter-bp, 2 groove). *The reader currently has no mental model of what DeepPBS ingests.* **Priority 1.**
2. **Augmented-fold diagram** — 1 crystal + N passing frames → one shared PWM label, drawn to foreshadow the
   confound; reused in Ch. 12. **Priority 1.**
3. **Stage-4 yield funnel** — fnat-pass states → helix-guard survivors → written `.npz` per pilot, making the
   silent-skip attrition visible. *Data-backed* (Stage-4 log line counts vs fnat-pass counts). **Priority 2.**

---

## Chapter 11 — Stages 6–7: training & evaluation

**Reuse [exists ✓]:**
- `figure_scripts/fig3_box_pearson` and `fig4_box_mae` — per-entry metric distributions over the 130-entry
  benchmark, three arms (baseline / augmented·frozen / augmented·relaxed). The "what the eval produces"
  reference; `fig3`'s own title already makes the honest point ("the three arms barely differ against the
  entry-to-entry spread"). Verified: 12 pilots, three-arm boxes.

**Caption caveat — do NOT use here:**
- `figure_scripts/fig9_mixedmodel_effects` **exists on disk** and looks like a significance result. It is the
  **pseudoreplicated** forest plot (seed is the experimental unit, n = 5) and is reserved **only** as the
  negative example in Ch. 13. Never place it in Ch. 11 as a positive result. Its honest counterpart is
  `fig2_augmentation_delta`.

**To create:**
1. **The two-p-values contrast** (data) — the csl **own-family** example drawn twice: pooled entry×seed t-test
   (ΔP +0.199, p = 0.001) beside the five per-seed deltas (`−0.024, −0.043, +0.159, +0.162, +0.164`; 2/5
   negative). The single clearest pseudoreplication picture; shared with Ch. 13. **Priority 1.**
2. **Paired-training schematic** — one seed → shared init/shuffle → arms differ only by `data_dir`; the
   `--array=0-1` pair and the `_sN` seed stack as the experimental design. **Priority 2.**
3. **Stage-6/7 flow diagram** — config pair → array train → checkpoint-on-validation-MAE → id.txt eval →
   denominator-intersection → paired Δ. **Priority 2.**

---

## Chapter 15 — Structural fidelity

**Fully covered by existing figures [exists ✓] — reuse and caption, nothing to create:**
- `I1_iRMSD_distributions`, `I2_iRMSD_seg` (localized vs distributed distortion), `I4_interface_geometry`,
  `R1_ca_rmsd_stages`, `R2_per_residue_profiles` (rigid cores / floppy termini), `R3_minimization_delta`,
  `S1_bestworst_ets1`, `S1_bestworst_lef1`, `S2_stage_progression_ets1`.
- Also available and on-topic: `I3_stage2_to_stage3`, and the `rmsd_analysis/plots/` deep-dive set
  (`sampling_quality`, `stage_ecdf`, `state_trajectories`, `summary_by_family`, `variant_agreement`, …) for a
  successor appendix.
- **Reuse from the Ch. 9 ensemble suite** (already built this session): fnat violins / pass-rate ladder,
  backbone-vs-sidechain, interface-tightening, fnat-vs-iRMSD — same per-state metrics, rigidity-ordered.

---

## Chapter 16 — The augmentation effect (the productive null)

**Reuse [exists ✓]:**
- `M2_where_the_signal_lives` — **the centerpiece.** Four panels: own vs other family, seed-level within-pilot
  contrast, the headroom control (ρ = −0.27), and the dynamic-vs-static axis. Viewed and confirmed.
- `P3_augeffect_by_family` ⭐ — ΔPearson by **motif-level** family; viewed, correctly labeled "NOT per-pilot,"
  IRF/ETS the only net-positive. `P2_baseline_by_family`, `P1_family_table`,
  `figure_scripts/fig1_three_arm_accuracy`, `fig6_within_family_transfer`.

**Caption caveat:** do **not** use the earlier per-PDB-mislabeled family figure (ETS1–RUNX1 co-crystals misfiled
as Runt, inflating the ETS median). `P3` on disk is the corrected motif-level version.

**To create:**
1. **Cross-benchmark-vs-own-family scatter** (data) — each pilot's pooled cross-benchmark ΔP (≈0) against its
   own-family ΔP, ETS1 flagged. The one picture that shows the signal the average hides. **Priority 1.**
2. *(optional)* the **pseudoreplication contrast** if not already placed in Ch. 11 (shared figure).

---

## Chapter 17 — The mechanism story (revised)

**Reuse [exists ✓] — but mind the ordering:**
- `M4_dnaflex_vs_effect` — **the climax.** Viewed: three panels, bend-IQR ρ = −0.76 (relaxed, the only CI
  excluding zero) vs the flat static-crystal-bend control (ρ = +0.04). Its subtitle already states the finding.
  *This carries the result; feature it.*
- `M3_samefamily_both_arms` — own-family effect both arms + the frozen→relaxed shift (8/11 up). Viewed.
- `M1_apo_holo_mechanism` ⭐ — visually strong but plots the **protein-side** axes and hard-codes the static
  "TBP best reach, hurt most" reading. **Caption it as the hypothesis it replaced, not the result** — on its
  own it implies the story the data overturned.

**To create:**
1. **Static-fails / dynamic-works pair** (data) — own-family effect vs static crystal deformation (ρ = +0.04,
   flat) beside vs ensemble bend IQR (ρ = −0.76, relaxed). A tighter, talk-ready distillation of `M4`'s message.
   **Priority 1.**
2. **Hypothesis-filled-with-data quadrant** (schematic) — the Ch. 4 2×2 (reachability × DNA *fluctuation*) with
   pilots at measured positions, TBP flagged as the decisive control. Ties Ch. 4 motivation to Ch. 17 result.

---

## Chapter 18 — DNA relaxation

**Reuse [exists ✓]:**
- `dna_relax/dna_shape_features` (6-panel overview, 7 TFs), `tbp_dna_shape` (4-panel TBP), `crossfamily_bend`,
  the per-TF `<tf>_dna_shape` set (11 pilots, no hsf/irf), `mgw_fl_all12_panels`, `perposition_minorgroove`,
  `pycurves_ensemble_summary`, and the interactive `pycurves_viz/*.html` viewers.
- `M3_samefamily_both_arms` doubles as the frozen→relaxed effect figure here.

**Caption caveat — important, from the data check:**
- `crystal_convergence_bootstrap` must **not** be captioned "DNA bends toward the bound conformation." The
  whole-duplex convergence data do not support it: relaxation moves DNA *closer* to crystal in only 7/42
  feature×pilot cases, 2/7 on the bend median, and pyCurves flags only 3/13 pilots as `closer_to_crystal`
  (frozen *over*bends 6/7). Use this figure only if its **per-position** direction is checked directly; caption
  it as "per-position bend vs crystal, with CIs," not as whole-duplex convergence. (See Ch. 18.1.)
- Every DNA-relaxation figure carries the **partial-coverage** caveat: shape overview 7 pilots, per-TF panels 11
  (no hsf/irf), bend bootstrap 7, AF3/interface 12. Caption as preliminary.

**To create:**
1. **Falsifiable-prediction figure** (data) — predicted vs actual aug-sign shift for the induced-fit families
   (TBP, LEF1, homeodomain), frozen → relaxed, with the 8/11-upward result and TBP's zero-crossing flagged.
   **Priority 1.**

---

## Build order (if producing new figures)

**Tier 1 — the chapters that are unillustrated or whose argument needs a figure that doesn't exist:**
1. Ch. 10 feature-tensor schematic + augmented-fold diagram (2 schematics — Ch. 10 has nothing).
2. Ch. 11 two-p-values contrast (data — the pseudoreplication picture the whole stats story rests on).
3. Ch. 16 cross-benchmark-vs-own-family scatter + Ch. 17 static-fails/dynamic-works pair (data — both carry
   the *corrected* arguments and both come straight from CSVs already local).

**Tier 2 — schematics that improve orientation but aren't load-bearing:**
4. Ch. 10 yield funnel, Ch. 11 paired-training + flow schematics, Ch. 17 hypothesis quadrant, Ch. 18
   falsifiable-prediction figure.

**Tier 0 — no build, just decisions already made:** the two swaps (feature `M4` not `M1` in Ch. 17; `M2`
centerpiece in Ch. 16), the two "do not use here" caveats (`fig9` in Ch. 11; whole-duplex convergence caption in
Ch. 18), and captioning `M1` as hypothesis.
