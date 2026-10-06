# analysis/analyses/nmr_rmsd — NMR bundles vs BioEmu ensembles

Two layers live here.

## 1. Per-entry (original)

`compute_nmr_state_rmsd.py` + `plot_best_state.py` — every BioEmu conformation
against every NMR state of its own reference, and how often each state is the
closest one. One entry at a time; each entry uses whatever residues its own
sequence alignment happens to match.

## 2. Apo/holo pairs (`pairs.csv`, `pair_core.py`, `compute_pair_rmsd.py`)

`pairs.csv` is the single source of truth for pairs. Each row carries a
**`status`**:

- **`ready`** — both members have a `<pdb>_chain<X>_conformations/` BioEmu
  ensemble under `structures/nmr/`. These are what every script analyses by
  default. Five of them:

  | pair | apo | holo | family |
  |---|---|---|---|
  | trf1  | 1ity A | 1iv6 A | Telobox/Myb |
  | nhp6a | 1lwm A | 1j5n A | HMG-box |
  | vnd   | 1vnd A | 1nk2 P | Homeodomain |
  | maze  | 2mrn A | 2mru A | AbrB/MazE |
  | rok   | 5zuz A | 5zux A | Rok |

- **`needs_ensembles`** — 33 further pairs found by PDB search (solution-NMR
  protein+DNA entries matched to protein-only solution-NMR entries in the same
  95% sequence cluster), filtered to >=10 models on both sides, >=84% coverage
  of the shorter sequence, and DNA-binding regulators. BioEmu has not been run
  on them. They are skipped by default and do not clutter any run; name one with
  `--pair-ids` to analyse it the moment its ensembles land, then flip the status.

  `c02` (ETS1, `1r36`/`2stt`) is the cheapest: 25/25 models, full containment,
  and the 2stt ensemble already exists — only the apo side needs generating.
  Several others map onto existing pilots (c09 forkhead/foxa, c32 Rel/nfat,
  c12 HMG-box/lef1, c08 Telobox/trf1), which ties this validation to the
  systems in the benchmark.

  They are **not structurally vetted**: "apo" means the entry holds no nucleic
  acid, but it may carry a peptide or ligand, differ in oligomeric state, or
  cover different domain boundaries. Some share an apo member (1lqc appears 4x)
  — pick one holo before generating anything. Each pair costs **two** BioEmu
  runs (`tools/bioemu_from_pdb.sh`), which is the real expense. Expect the
  step-(1) gate to reject a substantial fraction: 2 of the current 5 failed it.

### Why a separate script

The per-entry script aligns each ensemble to its own reference, so the apo and
holo numbers of a pair land on different atom sets (57 CA for 1iv6, 67 residues
in 1ity). They are not comparable: a histogram built from them measures
construct length as much as conformation. `pair_core.py` builds **one** atom set
per pair shared by all four members — apo bundle, holo bundle, BioEmu(apo seq),
BioEmu(holo seq) — as

1. the 4-way sequence intersection (difflib, so ragged termini and numbering
   offsets do not matter), then
2. minus expression tags (poly-His, enterokinase, TEV) matched by regex on each
   member's own sequence. The intersection alone does not remove these when both
   members carry the tag — 5zux/5zuz share a C-term `LEHHHHHH`, and it would
   otherwise contribute pure noise to every RMSD, then
3. minus disordered termini: contiguous positions trimmed inward from each end
   while the intra-bundle spread in *either* NMR bundle exceeds
   `--spread-cutoff` (default 3.0 Å). Only the ends are trimmed, so interior
   mobile segments stay in the core (NHP6A keeps an internal loop reaching ~7 Å).

The core is defined by the two NMR bundles only — the BioEmu ensembles are what
is being tested and must not define the yardstick.

### Run

```bash
# the `pymol` conda env has mdtraj / scipy
~/miniconda3/envs/pymol/bin/python compute_pair_rmsd.py               # CA
~/miniconda3/envs/pymol/bin/python compute_pair_rmsd.py --atoms backbone
```

Outputs, all suffixed `_<atoms>.csv`, all tidy/long so plotting stays separate:

| file | contents |
|---|---|
| `pair_core_summary`     | per pair: core size, positions trimmed, states, frames |
| `pair_core_residues`    | per core position: residue ids in all four members, bundle spread, kept flag |
| `pair_bundle_rmsd`      | long: `apo_apo` / `holo_holo` / `apo_holo` state-pair RMSDs |
| `pair_bundle_summary`   | (1) medians, effect size, bootstrap CI, informative gate |
| `pair_frame_state_rmsd` | long: every (ensemble, bundle, state, frame) RMSD |
| `pair_frame_summary`    | (3) per BioEmu frame: d_min to each bundle, `delta_A` |
| `pair_2x2_summary`      | (2) per (pair, ensemble, bundle) cell: median d_min, coverage, state recall |

### (1) Is the pair informative?

Within-bundle heterogeneity (apo–apo, holo–holo) vs between-bundle (apo–holo)
separation, on the shared core. NMR bundle spread is a precision artefact as
much as dynamics, so a pair whose apo/holo separation does not clear its own
bundle noise cannot support any claim about BioEmu. Gate: the bootstrap CI on
`median(between) − median(within)` must exclude 0 **and** the ratio must reach
`--min-ratio` (default 1.25). States, not state-pairs, are the bootstrap
sampling unit; self-comparisons from resampling are dropped, not counted as 0.

CA results (backbone agrees to within 0.05 Å throughout):

| pair | within | between | separation [95% CI] | ratio | verdict |
|---|---|---|---|---|---|
| nhp6a | 0.50 | 1.55 | +1.05 [+0.96, +1.12] | 3.08 | **informative** |
| rok   | 0.75 | 1.75 | +1.00 [+0.89, +1.11] | 2.33 | **informative** |
| vnd   | 0.41 | 1.15 | +0.75 [+0.68, +0.81] | 2.83 | **informative** |
| trf1  | 1.11 | 1.28 | +0.17 [+0.09, +0.23] | 1.15 | gated out |
| maze  | 0.90 | 1.07 | +0.17 [−0.03, +0.31] | 1.19 | gated out |

trf1 and maze are gated out: their apo/holo difference is real but tiny relative
to bundle noise (trf1's CI excludes 0 only because n is large — the effect is
0.17 Å on a 1.1 Å noise floor). maze additionally has only 7 holo states.
**Do not build BioEmu conclusions on these two.**

### (2) The 2×2 and (3) delta

Both ensembles scored against both bundles; `delta_A = d_min(holo) − d_min(apo)`
per frame, positive meaning the frame sits closer to the free state.

| pair | ensemble | med d_min apo | med d_min holo | med delta | frames closer to apo |
|---|---|---|---|---|---|
| nhp6a | apo seq  | 0.95 | 1.28 | +0.34 | 95% |
| nhp6a | holo seq | 0.94 | 1.27 | +0.32 | 93% |
| rok   | apo seq  | 1.61 | 2.14 | +0.55 | 99% |
| rok   | holo seq | 1.34 | 2.04 | +0.62 | 99% |
| vnd   | apo seq  | 0.96 | 0.86 | −0.10 | 24% |
| vnd   | holo seq | 0.91 | 0.85 | −0.08 | 27% |

(trf1/maze omitted — gated out at step 1.)

Reading (this table covers the original five pairs only; it predates the
expansion to 16 and has not been regenerated), restricted to the three
informative pairs of that set: 93-99% of nhp6a and rok
frames sit closer to the apo bundle, vnd's sit closer to holo, and the apo-seq
and holo-seq ensembles agree throughout (the construct control behaves).

**But the sign of delta is not the whole story — see the projection below, which
overturns the obvious reading of this table.** delta counts which bundle is
nearer without asking whether the frame lies anywhere near the line between
them.

## 3. Where the ensembles actually sit (`sweep_trim.py`, `plot_pair_pca.py`)

### Trim sensitivity

`sweep_trim.py` re-runs everything across `--spread-cutoff` 2-5 A plus a no-trim
control (`trim_sweep_<atoms>.csv`, `plots/trim_sweep_<atoms>.png`).

- **2.5-4.0 A is a stable plateau** for every pair; the default 3.0 A sits in it.
  vnd's negative delta holds at -0.07 to -0.18 A right across it, so its sign is
  *not* a trimming artefact.
- **2.0 A is degenerate** — the core collapses to 7 positions (nhp6a) and 4
  (vnd). Do not use it.
- **No-trim is meaningless**, and flips signs (vnd -0.10 -> +1.59 A, nhp6a
  +0.34 -> -0.41 A). With 12-20 A of disordered tail left in, the RMSD measures
  tails. This is the control that justifies trimming at all.

### The deposited-subspace PCA (`plot_pair_pca.py`) — the main view

Superpose all four members on the shared core, fit a Cartesian PCA on the
**deposited states only** (apo + holo NMR), and project the BioEmu ensembles in
as test data. The bundles define the coordinate system because they are the
conformations that were actually observed; the ensembles never influence the
axes. Two numbers per structure, both in Angstrom and directly comparable to an
RMSD:

- **PC1 / PC2** — position inside the deposited subspace. For an informative
  pair PC1 comes out as the apo->holo direction on its own.
- **residual** — the part of the structure the subspace cannot reproduce,
  `||x - mean - P_k(x - mean)|| / sqrt(n_atoms)`, with k the number of PCs
  holding 90% of the deposited variance. Bundle residuals are computed
  **leave-one-out** (subspace refitted without that state), so the states that
  defined the space do not score a free zero.

This replaced a classical-MDS embedding, which kept only 35-67% of the positive
eigenvalue mass (RMSD after superposition is not a Euclidean metric) and had to
be read as a sketch. The PCA is metric and its variance fractions are honest.

**The subspace validates itself, and independently reproduces the step-(1) gate
a third time.** For the three informative pairs a single PC dominates and it is
the apo/holo axis; for the two gated-out pairs no direction dominates (table
covers the original five pairs only):

| pair | PC1 share of deposited variance | k for 90% | apo/holo split on PC1 |
|---|---|---|---|
| nhp6a | 78% | 4 | −5.60 / +5.25 |
| vnd | 75% | 4 | −2.57 / +3.93 |
| rok | 70% | 7 | −6.33 / +6.64 |
| trf1 | 25% | 15 | −2.01 / +2.01 |
| maze | 32% | 10 | +0.23 / +1.93 |

Where the ensembles land (CA; backbone agrees):

| pair | BioEmu PC1 (apo-seq / holo-seq) | % of the way apo→holo | residual | bundle residual (LOO) |
|---|---|---|---|---|
| rok | −5.88 / −6.32 | **3%** | 1.46 / 1.24 | 0.34 / 0.39 |
| nhp6a | −1.13 / −0.96 | **41%** | 0.80 / 0.80 | 0.31 / 0.22 |
| vnd | +1.66 / +1.47 | **65%** | 0.70 / 0.67 | 0.25 / 0.14 |

### What this shows

**Gate.** 14 of 16 pairs are informative; only trf1 and maze fail. The new pairs
also bring a much wider dynamic range — THAP1 and TrpR separate their bundles by
4.3-4.4 A, against 1.1-1.8 A for the original nhp6a/rok/vnd — which is what was
missing for judging whether BioEmu tracks a transition at all.

**Position along the apo->holo axis, all 14 informative pairs.** `%across` is
0% at the apo bundle, 100% at holo; `resid ratio` is the ensemble residual over
the deposited bundles' own leave-one-out residual, i.e. how far outside the
observed subspace the ensemble sits (values >4 mean `%across` is not worth
reading closely).

| pair | protein | %across | resid ratio |
|---|---|---|---|
| rok | Rok C-term DBD | 2% | 4 |
| c09 | FOXD3 / Genesis | 18% | 4 |
| c10 | SarA | 25% | 10 |
| nhp6a | NHP6A | 42% | 3 |
| c17 | YB-1 | 46% | 7 |
| c03 | c-Myb | 47% | 2 |
| c20 | THAP1 | 51% | 9 |
| c13 | Lac repressor | 58% | 3 |
| c04 | MeCP2 MBD | 61% | 2 |
| vnd | VND/NK-2 | 64% | 4 |
| c22 | Trp repressor | 65% | 3 |
| c01 | Mu repressor C | 84% | 4 |
| c11 | Dead ringer (ARID) | 85% | 5 |
| c05 | WRKY4 | 94% | 3 |

**This retires the conclusion drawn from the first three pairs.** On nhp6a, rok
and vnd alone the reading was "BioEmu reaches neither basin and never gets to
holo". At n=14 that does not survive: `%across` spans 2% to 94%, median 55%, and
restricting to the 8 pairs whose ensembles stay near the deposited subspace
(resid ratio < 4) changes nothing — median 60%, same 2-94% range. Two pairs sit
at the apo end, three at the holo end, nine in between. **BioEmu is not
systematically apo-biased; where its ensemble lands on the apo->holo coordinate
varies by system across the entire range.** c05 (WRKY4) at 94% is essentially at
the bound state — the case the three-pair set said did not occur.

What does survive is that the ensembles are consistently *wider* than the
deposited bundles and displaced outside their subspace, by 2-10x the bundles'
own leave-one-out residual. That displacement, not a directional apo bias, is
the reproducible BioEmu signature here.

The obvious follow-up is what predicts `%across`: it is not DBD family (two
homeodomain-adjacent pairs sit at 42% and 64%), and the two-ensemble control
rules out construct identity. Bundle separation, protein size, and the residual
ratio are all in `pair_pca_summary_*.csv` and `pair_bundle_summary_*.csv` if you
want to test those directly — with 14 points, expect to be able to reject strong
effects only.

### The two-ensemble control, and what the input sequences actually differ by

Both members of a pair are the same protein and BioEmu never sees DNA, so the
two ensembles should be draws from the same distribution. They are — and two
pairs make that a quantitative statement rather than an assumption.

Every construct difference is **flanking only**; no pair has an interior
substitution:

| pair | apo len | holo len | difference |
|---|---|---|---|
| nhp6a | 93 | 93 | **identical input sequence** |
| vnd | 77 | 77 | **identical input sequence** |
| trf1 | 67 | 57 | apo has `EKHRA` (N) + `SDSED` (C) |
| rok | 103 | 92 | apo has `AGIPD` (N) + `AESANE` (C) |
| maze | 67 | 50 | apo has `NHKVHHHHHHMSDDDDK` (N-term His-tag) |

For **nhp6a and vnd the two ensembles are independent BioEmu runs on a
byte-identical sequence**, so their disagreement is pure sampling noise. It is
small: PC1 medians differ by 0.17 A and 0.19 A, KS cannot separate them
(p = 0.64, 0.33), and the IQRs match (vnd 1.02 vs 0.96; the larger apo-seq *sd*
is two outlier frames, not a wider bulk).

That gives a measured noise floor of **~2-3% of the apo-holo axis**, which is
what the displacements in the previous section should be read against:

| pair | ensemble disagreement / apo-holo sep | displacement from apo |
|---|---|---|
| nhp6a | 0.016 | 0.41 |
| vnd | 0.029 | 0.65 |
| rok | 0.034 | 0.03 |

So nhp6a's and vnd's failure to reach holo is 15-25x the sampling noise, and
rok's parking at apo is within ~1 noise unit of the apo bundle. trf1's two
ensembles are indistinguishable (KS p = 0.84) despite a 10-residue difference,
and rok's are borderline (p = 0.057) despite 11 — flanking residues cost little.

**maze is the exception, and not for the reason the tag suggests.** Its two
ensembles match on PC1 (p = 0.23) but differ 1.7x in residual (2.12 vs 3.62 A).
The worse one is `holo_seq` — the *untagged* 50-mer — so this is truncation, not
the His-tag: both maze ensembles are poorly converged (PC1 sd 8-10 A against
1-5 A everywhere else), consistent with EcMazE being a partly disordered
antitoxin, and the bare 50-residue fragment is the worse of the two. maze is
gated out at step (1) regardless.

One caution for any new pair: a tag can be dropped from the measurement core
(`pair_core.py` does) but not from the *generation*. If BioEmu is run on a
tagged FASTA, it sampled a protein with a disordered tag attached and the core
conformation carries that. Strip tags from the input sequence, not just the
analysis.

### Known limits

- `state_recall_2A` saturates at 1.00 nearly everywhere: after terminal
  trimming, 2 Å is too loose a threshold on these cores. Use a sweep (item 4)
  rather than a fixed τ.
- Deposited NMR states are not Boltzmann-weighted. Everything here is *coverage
  of the deposited bundle*, never population agreement.
- Sixteen pairs analysed, fourteen informative (trf1 and maze gated out). The
  cross-pair spread is descriptive; with n = 14 only strong effects are
  testable. Twenty-two further pairs are staged in `pairs.csv` at
  `status: needs_ensembles`.
- PC2 holds only 5-10% of the deposited variance even for the informative
  pairs, so the vertical spread of the PCA panel carries little information —
  PC1 and the residual are what to read.
- The PCA subspace is fitted on 27-50 deposited states in a 3N-dimensional
  space, so k=4-18 components is a generous fit. The leave-one-out bundle
  residuals are the guard against reading too much into that.
- `delta_A` in `pair_frame_summary` is a sign statistic and ignores magnitude and
  direction. It is kept because it is what the per-frame tables naturally hold,
  but conclusions should come from t and p.
