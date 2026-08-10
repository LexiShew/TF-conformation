#!/bin/bash
# =============================================================================
# tools/git_commit_batches_2026-08-10.sh
#
# One-time history grouping: turns the large uncommitted working tree that
# accumulated between 2026-07-17 (ff9cc14) and 2026-08-10 into a sequence of
# topically-grouped commits on feature/dna-restraint-k.
#
# Run from the repo root:
#     cd /project2/rohs_102/shewchuk/TF-conformation
#     bash tools/git_commit_batches_2026-08-10.sh
#
# Idempotent-ish: every step is "git add <paths> && commit only if staged".
# A step whose paths are already committed becomes a no-op.
#
# Author identity comes from the repo's own git config
# (user.name=LexiShew, user.email=lexishew@live.com) -- not overridden here.
#
# NOTHING IS PUSHED. Review with `git log --stat` before `git push`.
#
# WHAT IS DELIBERATELY NOT COMMITTED (see step 1, .gitignore):
#   analysis/figures/pymol/prepost_min/   70 MB  regenerate: figscripts/render_prepost.sbatch
#   analysis/figures/pymol/_viewtest/      5 MB  scratch view tests
#   analysis/align_compare/{_work,logs}   82 MB  regenerate: align_compare/align_array_min.sbatch
#   docs/*.pptx                           15 MB  re-export from committed figures
#   *_auto*.inp (generated pyCurves decks) rebuilt by pycurves_array.sh
#   *.bak_* rollback snapshots, *.out slurm logs
# =============================================================================
set -uo pipefail
cd "$(git rev-parse --show-toplevel)"

STEP=0
# commit <<subject+body>> -- commits whatever is staged; skips cleanly if nothing is.
commit () {
    STEP=$((STEP+1))
    if git diff --cached --quiet; then
        printf '[%02d] SKIP (nothing staged): %s\n' "$STEP" "$(echo "$1" | head -1)"
        return 0
    fi
    local n; n=$(git diff --cached --name-only | wc -l | tr -d ' ')
    git commit -q -m "$1" && printf '[%02d] OK  %5s files  %s\n' "$STEP" "$n" "$(echo "$1" | head -1)"
}

# A -- add paths, tolerating ones that do not exist
A () { for p in "$@"; do [ -e "$p" ] && git add -A -- "$p"; done; return 0; }


# ---------------------------------------------------------------- 1. hygiene --
# Guarded so a re-run does not append the block twice.
if grep -q 'Added 2026-08-10' .gitignore; then
  echo "[01] .gitignore block already present -- not re-appending"
else
cat >> .gitignore <<'IGN'

# ---------------------------------------------------------------------------
# Added 2026-08-10 -- generated / regenerable output that should not be tracked
# ---------------------------------------------------------------------------

# pyCurves input decks: written next to the caller by pycurves_array.sh and
# af3_pycurves.sh, one per structure. Regenerated on every run.
*_auto.inp
*_auto_1.inp
*_auto_2.inp
# ...except the 20 hand-picked visualization decks already tracked:
!analysis/dna_relax/pycurves_viz/*_auto.inp

# Rollback snapshots taken by hand before an edit (.bak_precolor, .bak_10pilot,
# .bak_pre12, .bak_csl_nd, ...). Git history is the rollback mechanism.
*.bak_*

# Slurm stdout logs written next to the code that submitted them (*.err was
# already ignored). Vendored libs ship legitimate .out fixtures -- keep those.
*.out
!lib/**/*.out

# Scratch minimization/docking intermediates. Regenerate with
# analysis/align_compare/align_array_min.sbatch (+ align_compare_job.sbatch);
# the two summary CSVs at align_compare/ root ARE tracked.
analysis/align_compare/_work/
analysis/align_compare/logs/

# Bulk PyMOL renders. Per-pilot af3/bioemu panels ARE tracked; these two trees
# are 75 MB of overlay frames and camera tests, regenerated with
# analysis/figscripts/render_prepost.sbatch -> render_prepost_min.py.
analysis/figures/pymol/prepost_min/
analysis/figures/pymol/_viewtest/

# Exported slide decks: binary, non-delta-compressing, and rebuildable from the
# committed figures. (deck/TF_conformation_deck.pptx predates this rule and
# stays tracked.)
*.pptx

# Scratch structure dumps
*.tmp.pdb
IGN
fi
git add .gitignore
commit "Repo hygiene: ignore generated pyCurves decks, .bak_* snapshots, slurm logs, bulk renders, and exported decks

Between 2026-07-17 and 2026-08-10 the working tree accumulated ~1680 uncommitted
paths, most of them regenerable output. This establishes what is and is not
tracked before committing the rest.

Newly ignored, with the regeneration path for each:
  *_auto*.inp (1366 files)      rebuilt by pycurves_array.sh / af3_pycurves.sh
  *.bak_* (18)                  hand rollback copies; git history replaces them
  *.out                         slurm logs written next to submitting scripts
                                (lib/**/*.out negated -- vendored x3dna fixtures)
  align_compare/{_work,logs}/    82 MB  align_array_min.sbatch
  figures/pymol/prepost_min/     70 MB  render_prepost.sbatch
  figures/pymol/_viewtest/        5 MB  camera-angle scratch
  *.pptx                         15 MB  re-export from committed figures
  *.tmp.pdb

Kept tracked despite bulk, because regeneration costs compute-node hours rather
than a render call: analysis/dna_relax/pycurves/ JSON (147 MB new, joining 2570
files already tracked) and the per-pilot af3/bioemu render panels."


# ------------------------------------------------------- 2. canonical palette --
A palette.py analysis/scripts/plot_diversity.py rmsd_analysis/plot_rmsd_progression.py \
  rmsd_analysis/plots
commit "palette.py: single canonical colour spec, and rewire existing figure scripts to it

One hue per entity, so colour is the cross-reference between figures rather than
a per-script choice:
  GREY     baseline model / crystal reference   the fixed reference
  TEAL     BioEmu ensemble / augmented-FROZEN   one narrative thread (BioEmu is
                                                the source of the frozen arm)
  GREEN    augmented-RELAXED DNA                the intervention on that thread
  AF3      AlphaFold3 comparator (lavender)     a different tool, Ch. 3 only
  ALARM    rose                                 annotation only, never a series

Rewired to import from it instead of carrying local hex literals:
  analysis/scripts/plot_diversity.py        C_BIO/C_AF3 + median label colours
  rmsd_analysis/plot_rmsd_progression.py    TF_PALETTE, VARIANT_COLOR
rmsd_analysis/plots/ (22 PNGs) regenerated under the new palette."


# --------------------------------------------------- 3. pilot expansion 12/13 --
A config/pilots wrappers
commit "Add hsf/irf pilot configs (frozen + dnarelax) and raise stage-3 walltimes

config/pilots/{hsf,irf}{,_dnarelax}.sh bring the pilot set to 13 structures.
hsf.sh documents a non-obvious structural fix: 5D5U's asymmetric unit deposits
one HSF1 protomer + one 12-nt strand (GGTTCTAGAACC, exactly self-complementary);
the biological assembly is tetrameric with the duplex generated by the
crystallographic 2-fold (oper 2, x,-y,-z). The bare ASU strand is not a duplex,
so DeepPBS featurization failed with 'helix count problem 0' on every frame.
5d5u.cif now holds the regenerated monomer assembly (one protein chain over the
full 12-bp duplex, chain A + 2-fold mate D); the original ASU is preserved as
5d5u_asu_original.cif.

Walltime bumps -- the larger pilots were hitting the wall:
  wrappers/stage3_array.sh    00:30:00 -> 01:30:00
  wrappers/stage3_recover.sh  04:00:00 -> 08:00:00"


# ---------------------------------------------------------- 4. pyCurves batch --
A analysis/dna_relax/pycurves_batch analysis/dna_relax/pycurves_worklist.tsv \
  analysis/dna_relax/pycurves analysis/pycurves
commit "pyCurves batch: extend worklist to all 13 pilots, add helical parameters for 6 new TFs

build_worklist.sh PID map grows from 7 to 13 entries (adds csl=3brg, err=1lo1,
nfat=1a66, runx=1hjc, hsf=5d5u, irf=1if1); pycurves_worklist.tsv regenerated
(+1136 rows) covering crystal + frozen + relaxed states per pilot.

New pyCurves output (curvesplus + legacy axis conventions, one JSON pair per
structure) under analysis/dna_relax/pycurves/:
  csl 348 - err 366 - hsf 392 - irf 370 - nfat 350 - runx 388 files
plus the relaxed-state JSONs that completed for dux4 (14), egr1 (6), foxa (10).

These are tracked rather than regenerated on demand: the array job runs on a
compute node (pyCurves' JAX thread pool exceeds the login-node RLIMIT_NPROC) and
the runner skips rows whose output JSON already exists, so the committed files
are also the resume state."


# ------------------------------------------------- 5. DNA shape generalization --
A analysis/dna_relax/scripts analysis/dna_relax/README.md analysis/dna_relax/data \
  analysis/dna_relax/figures
commit "DNA-shape analysis: auto-detect base-pair register so it runs on any pilot

batch_dna_shape_v2.py generalizes batch_dna_shape.py, which hardcoded the 1tgh
register. The antiparallel base-pair register is now detected from the docked
reference: strands are the first two chains sorted, strand B is paired to strand
C positionally, and duplex orientation is resolved geometrically from the C1'
endpoints. Metric definitions are unchanged -- validated byte-for-byte against
the original TBP CSVs (tbp_dna_{perstate,perres,perbp}, 0 diffs). Also emits
<prefix>_register.json (n_bp + strand ids) for the plotter.

plot_dna_shape.py rebuilds the 4-panel figure for any pilot from the CSVs:
(a) DNA backbone RMSD vs docked, (b) per-residue mean P displacement split by
strand, (c) max adjacent P-P gap with the 9 A unwind line, (d) delta global bend
vs docked. Panel titles are computed per TF rather than written for TBP.

Ran for 11 pilots: data/<tf>_dna_{perstate,perres,perbp}.csv + _register.json
and figures/<tf>_dna_shape.png. Relaxed-state counts reflect the pipeline
snapshot at run time -- csl and err were still being populated.

Also here: pycurves_ensemble_summary, crystal_vs_ensemble_mgw, bend/shape
fidelity panels, and the aggregate CSVs (pycurves_all_summary,
pycurves_bend_summary, mgw_fl_summary)."


# --------------------------------------------------------- 6. interface MGW-FL --
A analysis/dna_relax/interface_mgwfl.py analysis/dna_relax/interface_mgwfl.README.py
commit "Interface MGW-FL: minor-groove width restricted to contacted base pairs

interface_mgwfl.py computes minor-groove-width fluctuation over only the base
pairs the protein actually contacts, rather than the whole duplex, and relates
it to per-pilot augmentation effect. interface_mgwfl.README.py carries the
method notes (contact definition, FL statistic, how the register is reused from
the pyCurves _register.json). Runs end-to-end from the repo root under the
deeppbs env. Outputs iface_mgwfl_vs_accuracy.csv and the companion figures
committed in the preceding DNA-shape commit."


# ------------------------------------------------------------ 7. AF3 expansion --
A af3/input af3/run_af3_batch2.sh af3/af3_pycurves.sh af3/af3_dna \
  analysis/scripts analysis/data/af3_rmsd_to_crystal.csv \
  analysis/data/ensemble_diversity.csv analysis/data/ensemble_diversity_pairwise.csv \
  analysis/data/rmsd_to_crystal_af3_vs_bioemu.csv \
  analysis/data/rmsd_crystal_and_pairwise_summary.csv
commit "AF3 comparator: 7 new pilot inputs, DNA-duplex extraction, and 12-pilot diversity comparison

Input JSONs (dialect alphafold3 v1, modelSeeds=[1,2]) for csl/dux4/err/hsf/irf/
nfat/runx, built by pulling the DBD protein sequence from
structures/source_chains/<pid>_chains/*_protein.pdb, the duplex strands from
<pid>_dna.pdb, and structural metals from <pid>.cif as CCD ligands. Per-pilot
care that is easy to get wrong: runx deposits two identical duplex copies (keep
one, B/C); hsf's self-complementary duplex is supplied as the same strand twice;
err carries 2 Zn ligands.

run_af3_batch2.sh runs the batch under Apptainer on the qcbgpu l40s partition
against /project2/rohs_102/share/alphafold3. af3_pycurves.sh then extracts each
model's DNA and runs pyCurves in legacy axis convention -- af3/af3_dna/ holds
the per-seed/per-sample duplex PDBs + legacy JSON for 13 pilots (254 files).

Diversity comparison extended 6 -> 12 pilots (af3_rmsd_diversity.py,
plot_diversity.py, montage_ensembles.py PILOTS/ORDER/LABELS) with refreshed
ensemble_diversity{,_pairwise}.csv and af3 vs bioemu RMSD-to-crystal tables.

render_ensembles.py: crystal DNA switched from gold ring-mode cartoon to a
white 0.55-transparency surface, and the camera now looks perpendicular to both
the protein->DNA centroid axis and the DNA helical axis (SVD of the duplex
coordinates), so the contact face is toward the viewer. Plain orient() was
landing end-on down the duplex or putting the protein in front of the
interface."


# ----------------------------------------------------- 8. figscripts F/I/R/S/P/M --
A analysis/figscripts analysis/figures/F5_align_dna_displacement.png \
  analysis/figures/D1_diversity.png analysis/figures/F1_fnat_distributions.png \
  analysis/figures/F2_passrate_bars.png analysis/figures/F3_fnat_vs_iRMSD.png \
  analysis/figures/F4_interface_size.png analysis/figures/I1_iRMSD_distributions.png \
  analysis/figures/I2_iRMSD_seg.png analysis/figures/I4_interface_geometry.png \
  analysis/figures/M1_apo_holo_mechanism.png analysis/figures/P1_family_table.png \
  analysis/figures/P2_baseline_by_family.png analysis/figures/P3_augeffect_by_family.png \
  analysis/figures/R1_ca_rmsd_stages.png analysis/figures/R2_per_residue_profiles.png \
  analysis/figures/R3_minimization_delta.png analysis/figures/pymol \
  analysis/data/perstate_metrics.csv analysis/data/perentry_accuracy.csv \
  analysis/data/mechanism_apo_holo.csv analysis/data/reachability.csv \
  analysis/data/ca_rmsd_perresidue.csv
commit "figscripts/: pilot-agnostic F/I/R/S/P/M figure suite -- pilots discovered from disk

New analysis/figscripts/ rebuilds the structural and accuracy figure series from
pipeline output with no hardcoded pilot list: fig_common.discover_pilots() walks
the output tree, so a newly-run TF appears in every figure once its pipeline
output and eval JSON exist. Adding a pilot means editing PILOT_META only.

  fig_common.py             discovery, PILOT_META, pass-rate ordering, savefig
  extract_metrics.py        pipeline output -> perstate_metrics.csv,
                            perentry_accuracy.csv, mechanism_apo_holo.csv
  compute_reachability.py   coordinate pass for M1's reachability axis
                            (d_min/reach_ratio/rmsf_mean) and R2's per-residue
                            profile -- runs in the pycurves env
  compute_align_displacement.py  interface- vs global-alignment DNA displacement
  make_F/I/R/S/P/M.py       the figure series
  make_af3flex/bend/aligncompare/clash_trajectory.py  companion panels
  render_prepost_min.py + render_prepost.sbatch / render_S.sbatch  PyMOL renders

Regenerated F/I/R/P/M/D panels and the metric CSVs behind them (perentry
+2342 rows, perstate +1951, ca_rmsd_perresidue +2353 as pilots 7-12 land).
PyMOL: per-pilot af3/bioemu panels for csl/err/nfat/runx/dux4/hsf complete the
already-tracked series; the 70 MB prepost_min/ overlay tree and _viewtest/
camera scratch are ignored (regenerate via render_prepost.sbatch)."


# --------------------------------------------------------- 9. alignment compare --
A analysis/align_compare
commit "Alignment-mode comparison: interface vs global docking, DNA displacement

align_compare_job.sbatch + align_array_min.sbatch re-dock and minimize each
pilot's ensemble under interface-aligned vs globally-aligned superposition, to
test whether the production choice (interface alignment) moves the DNA relative
to its crystal pose more or less than the alternative.

Tracked results: dna_placement_by_mode.csv and
dna_displacement_interface_vs_global.csv (the F5 figure is in the figscripts
commit). The 718 intermediate docked/minimized PDBs and 75 array logs under
_work/ and logs/ are ignored -- rerun align_array_min.sbatch to rebuild them."


# --------------------------------------------------------------- 10. mechanism --
A analysis/mechanism
commit "Mechanism analysis: the augmentation signal is dynamic, and lives in each pilot's own family

Reanalysis of the induced-fit vs conformational-selection hypothesis using the
measured DNA-shape data now on disk, replacing the curated literature-derived
deformation labels. Everything is computed from pipeline output; no numbers are
carried forward from prose. FINDINGS.md is the write-up (status: provisional).

Three results:
1. The 130-entry cross-benchmark delta-Pearson averages the effect away. It
   lives in each pilot's OWN family: ETS1 own-family delta-Pearson = +0.111
   (p = 0.013, 0 of 6 seeds negative) against a cross-benchmark -0.002.
2. Not explained by baseline headroom (rho = -0.27, p = 0.42) or own-family
   subset size (rho = -0.32, p = 0.34).
3. The discriminating DNA variable is dynamic, not static: ensemble DNA-bend IQR
   tracks the own-family effect at rho = -0.53, while crystal bend -- the static
   deformation the original hypothesis was written on -- gives rho = +0.04
   (p = 0.92).

Scripts: mechanism_analysis.py (main), mechanism_probe.py,
mechanism_confound.py (headroom/size/leakage controls), pilot_audit.py,
make_mechanism_fig.py, make_samefamily.py, make_dnaflex{,_labeled}.py.
Figures M2-M9 + correlation matrix + summary diagram + per-pilot samefam PWM
panels; 18 result tables incl. bootstrap_ci.json, confound_tests.csv,
regime_contrast{,_test}, static_vs_dynamic.csv."


# ---------------------------------------------------- 11. figure_scripts 1-9 --
A analysis/figure_scripts
commit "figure_scripts/: three-arm benchmark figures 1-9, plus a review flagging fig9 as pseudoreplication

Regenerates the baseline vs augmented-frozen vs augmented-relaxed comparison
broken out by Pfam family. _common.py holds the shared load/colour/label layer;
each make_figN_*.py writes one PNG.

The core design decision, documented in README.md: because the frozen and
relaxed pipelines each retrained their own baseline, absolute Pearson is not
comparable across pipelines, so the cross-treatment quantity is the
within-pipeline seed-paired delta-Pearson = aug_sN - base_sN. figs 2 and 5 use
exactly that. Data inputs: perentry_condition.csv, perseed_summary.csv,
perseed_perentry.csv (12 pilots x 130-entry general benchmark x 5 seeds).

REVIEW_figure_scripts.md is a reviewer pass that re-fit the models
independently, and it finds a statistical error worth reading before anyone
cites fig9: make_fig9_mixedmodel_effects.py fits metric ~ C(arm)*C(dna) with
crossed random INTERCEPTS on entry and seed, estimates seed Var ~ 0, and
therefore treats all 50 (csl) / 100 (runx) entry x seed rows as independent
replicates of the treatment contrast -- shrinking the arm SE about 2x (csl 0.028
vs 0.048 seed-level) and yielding csl p=0.003, runx p<0.001. But augmentation is
applied once per seed (5 retrainings), so the replication unit is n=5 seeds. The
correct analyses -- seed-level paired t, or a random SLOPE of arm across seeds
(re_formula='~arm', groups=seed) -- give csl p=0.076-0.16 and runx p=0.128-0.20,
i.e. not significant. The script is committed as-is for provenance; fig9 should
not be reported as written. The other figures' error bars are honest."


# ------------------------------------------------------- 12. atom importance --
A analysis/interpret_tfconf.py analysis/interpret_tfconf_all.py \
  analysis/compare_importance.py analysis/compare_importance_all.py \
  analysis/make_importance_figures.py analysis/INTERPRET_SUITE_README.md \
  analysis/ATOM_IMPORTANCE_DELIVERABLES.md analysis/IMPORTANCE_ANALYSIS_RESULTS.md \
  analysis/importance_boxplot_comparison.png analysis/importance_scatter_comparison.png \
  analysis/importance_shift_distributions.png analysis/importance_comparison_table.csv \
  analysis/importance_stats.json ALL_PILOTS_SUBMISSION_STATUS.md \
  batch_interpret_all.sh run_all_interpretations.sh submit_interpret.sh submit_figs.sh
commit "Atom-importance attribution suite: occlusion-based DeepPBS interpretability across 3 arms

Measures how much each protein atom contributes to the predicted binding
specificity by masking it and recording the MAE shift in the predicted PWM, then
compares the three arms per pilot (baseline / augmented-frozen /
augmented-relaxed) to ask whether augmentation redistributes which atoms the
model reads.

  interpret_tfconf.py / _all.py       single-pilot and 12-pilot occlusion runs
  compare_importance.py / _all.py     cross-arm comparison and statistics
  make_importance_figures.py          boxplot / scatter / shift-distribution
  batch_interpret_all.sh              submits one GPU job per pilot
  run_all_interpretations.sh          submit, wait, then compile
  INTERPRET_SUITE_README.md           how to run it end to end
  ALL_PILOTS_SUBMISSION_STATUS.md     2026-07-30 submission record (job ids)

Results: importance_comparison_table.csv, importance_stats.json, the three
comparison figures, and the ATOM_IMPORTANCE_DELIVERABLES /
IMPORTANCE_ANALYSIS_RESULTS write-ups."


# -------------------------------------------------- 13. stage-3 minimization figs --
A analysis/stage3_figs
commit "Stage-3 minimization figure set: 20 panels + log parser

parse_stage3_logs.py turns the OpenMM minimization logs into
data/{summary,quantiles,recovery,traj_sample}.csv; figures fig01-fig20 cover
clash-ensemble and energy trajectories, initial-vs-final clash, collapse stage,
runtime, backbone vs sidechain motion, interface stage2->stage3, strain vs fnat,
cage drift, metals vs clusters, recovery, frozen vs relaxed DNA, the cage
schematic and two-force timeline, plus the fnat violin / pass-rate ladder /
ensemble coverage / fnat-vs-iRMSD / fnat ECDF summaries."


# --------------------------------------------------------- 14. stage-7 A/B PWM --
A stage7_eval
commit "Stage-7 A/B PWM eval: per-pilot id files and GPU runner

run_ab_pwms.sbatch scores the baseline and augmented models on matched entry
lists (ab_idfiles/id_<pilot>_ab.txt, 9 pilots) under the deeppbs env on an
rtx5000, so the two arms are compared on identical entries."


# ------------------------------------------------------- 15. rmsd regeneration --
A rmsd_analysis
commit "rmsd_analysis: batch regeneration job and refreshed stage-progression tables

rmsd_regen.sbatch recomputes the stage-wise RMSD tables on a compute node under
the pycurves env with BLAS threads pinned to 1 (the login-node thread pool
segfaults otherwise). per_state_rmsds.csv (+11178 rows) and summary_stats.csv
now cover the expanded pilot set."


# --------------------------------------------------- 16. frozen vs relaxed grid --
A analysis/frozen_vs_relaxed_grid.csv analysis/frozen_vs_relaxed_slope.png \
  analysis/augmentation_paired_stats.csv analysis/trial_matrix.csv
commit "Frozen vs relaxed DNA: full paired grid and seed-level augmentation statistics

frozen_vs_relaxed_grid.csv is the completed frozen x relaxed comparison across
pilots, with frozen_vs_relaxed_slope.png as the per-pilot summary.
augmentation_paired_stats.csv carries the seed-paired tests (n = 5 seeds, the
correct replication unit -- see REVIEW_figure_scripts.md), and trial_matrix.csv
indexes which arm/DNA/seed combinations exist on disk."


# ----------------------------------------------------- 17. narrative figures --
A analysis/figures/struct_pwm analysis/figures/chapter_figs
commit "Structure+PWM exemplars and chapter schematics

struct_pwm/: BioEmu ensemble render beside crystal-target / baseline / augmented
PWM logos for ETS1 (3wty, delta-r +0.06), TBP (2ko0, the induced-fit case,
delta-r +0.33) and CSL/RBPJ (6qhd, +0.10) -- plus struct_pwm_context.png, which
exists to keep the exemplars honest: it shows the per-entry delta-r distribution
(104 entries/pilot) with the exemplars marked as the high-gain tail and the
population medians near zero. Also the positive/negative variants and a README
mapping each file to its entry.

chapter_figs/: schematics for the write-up -- augmented-fold construction,
feature tensor, stage-4 yield, paired training, stage 6/7 flow, the hypothesis
quadrant, and the falsifiable prediction."


# ----------------------------------------------------------- 18. inventories --
A analysis/RESULTS_INVENTORY.md analysis/CONDITIONS_INVENTORY.md \
  analysis/results_inventory.csv analysis/conditions_axes.csv \
  docs/TF_conformation_Ch9_minimization.md \
  docs/TF_conformation_Ch10-11_featurization_training.md \
  docs/TF_conformation_Ch15-18_results.md \
  docs/TF_conformation_Ch10-18_figure_plan.md \
  docs/TF_conformation_book_outline.md docs/dna_shape_recognition_story.md
commit "Inventories and chapter drafts: what was actually run vs what is only capability

CONDITIONS_INVENTORY.md enumerates every experimental knob with the distinction
that matters when reading the results -- RUN (output on disk) vs CAPABILITY ONLY
(script/flag exists, no output tree). Notably: DNA-relax is run for all pilots
at the data stage but trained for 2 and evaluated for 1 at survey time; the
metal-cage A/B exists only in the May 2026 old_results tree; FOLD=0 exclusively;
the fnat gate was never varied off 0.5; and 4 DeepPBS pred_configs exist but the
benchmark used the full prot_shape model only.

RESULTS_INVENTORY.md maps the stage-1..6 output tree, and the two CSVs give the
same information in machine-readable form.

docs/: chapter drafts for minimization (Ch9), featurization and training
(Ch10-11), results (Ch15-18), the figure plan, and the DNA-shape recognition
narrative. Exported .pptx decks are ignored -- rebuild them from the committed
figures."


# ---------------------------------------------------------------- 19. book/ --
A book
commit "book/: chapter-by-chapter outline built from a read-only repo survey

OUTLINE.md v0.1 -- 6 parts, 21 chapters, 7 appendices. Each chapter lists its
argument, subpoints with concrete numbers, the figures that already exist (with
repo paths), and the figures still to be made (IDs N<chapter>.<n>).

README.md records how it was built, which is what makes it auditable: a survey
of all 20 project markdown docs, the four figure inventories, the result tables
in analysis/data/ and analysis/dna_relax/data/, deck_spec.json's 29 slide
titles, and git log from 2026-04-24 to 2026-07-30 for chronology. Numbers quoted
in the outline were re-derived from the CSVs on disk rather than copied from
prose, so a few disagree with older write-ups -- the CSVs win."


# ------------------------------------------------------- 20. this script itself --
A tools
commit "tools/: record the script that grouped the 2026-07-17..2026-08-10 backlog into commits

git_commit_batches_2026-08-10.sh is the exact sequence used to turn ~1680
uncommitted paths into the preceding commits, including the .gitignore rules and
the reasoning for what was excluded. Kept so the grouping is reproducible and
auditable rather than a shell-history artifact."


# ------------------------------------------------------------------- report --
echo
echo "=============== remaining uncommitted (should be 0) ==============="
git status --porcelain=v1 | head -20
echo "count: $(git status --porcelain=v1 | wc -l)"
echo
echo "=============== commits created ==============="
git log --oneline ff9cc14..HEAD
