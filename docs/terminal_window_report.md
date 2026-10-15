# §1.1 — Terminal-window bend. The gate resolves, and not to any of the four rows.

**Run:** 2026-09-29, from pyCurves output on disk. **Outcome: closest to row 3 — neither the ends nor
the end-to-end orientation reproduces −0.88 — but with a finding that changes what the conflict means.
Recommendation: do not spend the 24 arms.**

The decisive result is not any of the three terminal variants. It is the **control**.

---

## What ran

Twist-aware accumulation (cumulative Ω carried through every step) applied to:

1. **each end separately** — first W steps, last W steps
2. **mean of the two ends**
3. **ends in a common frame** — all steps accumulated with the *central* W wedges zeroed but their
   **twist retained**, so the two ends stay in one reference frame. This is the variant-3 comparison
   the plan asked for: relative orientation of the two ends, centre removed.
4. **all steps** — the control. The twist-aware reconstruction of whole-duplex bending from the same
   local parameters.

**Coverage note.** W = 6 needs ≥ 15 bp for two disjoint end windows plus a centre, leaving 2 pilots
after the ΔP join — uninterpretable, and the ρ = ±1.000 entries in the raw table are that artifact, not
a result. **W = 4, n = 9 pilots, is the usable test.** On these same 9, the benchmark to beat is
whole-duplex IQR × √L at **ρ = −0.883 (p = 0.0016)**.

## Result: nothing terminal reproduces it

| variant | relaxed ρ | p | frozen ρ | p | ρ vs length |
|---|---|---|---|---|---|
| 5′ end (W=4) | -0.267 | 0.488 | +0.033 | 0.932 | -0.769 |
| 3′ end (W=4) | -0.367 | 0.332 | -0.317 | 0.406 | -0.077 |
| mean of ends | -0.667 | 0.050 | -0.033 | 0.932 | -0.522 |
| ends, common frame | -0.233 | 0.546 | -0.533 | 0.139 | -0.094 |
| all steps (control) | +0.083 | 0.831 | +0.133 | 0.732 | +0.128 |

The best terminal variant is the mean of ends at **−0.667 (p = 0.050)**, which (a) does not reach the
−0.88 benchmark, (b) is **not reproduced in the frozen arm** (−0.033), and (c) weakens to −0.563
(p = 0.146) controlling duplex length — it is itself 52% length. Variant 3, the end-to-end relative
orientation, gives −0.233 relaxed / −0.533 frozen, neither significant.

**Per the plan's own rule, the signal is in neither the centre nor the ends as constructed.**

## The control is the real finding

`full` accumulates **every** step, twist-aware — the local-parameter reconstruction of exactly what
`bend_uu` is supposed to measure. If `bend_uu` IQR tracked real molecular bending, these two must agree.

> **full-reconstruction IQR vs `bend_uu` IQR: ρ = +0.050 (p = 0.90).**

They are unrelated. The `full` axis is reliable (arm-to-arm +0.900) and its effect correlation is
+0.083 relaxed / +0.133 frozen — flat, like the central window.

So the position taken in the v7 permutation report — *"two instruments looking at different parts of
the molecule"* — is **wrong, and I should correct it**. The central window and the terminal windows are
partial views, but `full` is not: it is the same molecule, the same states, the same local parameters,
completely accumulated. **Every reconstruction from local base-pair-step geometry is flat against ΔP.
Only `bend_uu` itself correlates.** That is no longer a disagreement between physical regions; it is a
disagreement between `bend_uu` and the local geometry it is nominally a summary of.

## What `bend_uu` IQR is actually tracking

The conditional footprint follow-up, run anyway because the control made it more informative:

| axis | vs footprint width (`n_iface`) | vs contact count (`n_iface_res`) |
|---|---|---|
| **whole-duplex `bend_uu` IQR** | +0.210 (p = 0.59) | **+0.849 (p = 0.004)** |
| terminal IQR (mean of ends) | +0.297 (p = 0.44) | +0.672 (p = 0.047) |
| full reconstruction | +0.358 (p = 0.34) | +0.185 (p = 0.63) |

**`bend_uu` IQR tracks the number of protein residues at the interface better than it tracks any DNA
geometry** — including a complete reconstruction of DNA geometry from local parameters. Contact count
also correlates with duplex length (−0.772), so the three are entangled.

This does **not** rescue the axis as a protein-side property, and the mediation test says so plainly:

| | unadjusted | controlling contact count |
|---|---|---|
| `bend_uu` IQR vs ΔP | −0.733 (p = 0.025) | −0.609 (p = 0.109) |
| IQR × √L vs ΔP | −0.883 (p = 0.002) | **−0.873 (p = 0.005)** |

Contact count does not mediate the ×√L correlation. And the reverse partial — contacts vs ΔP
controlling the axis — is +0.212 / +0.493, both n.s. So contact count is associated with the *axis* but
is not the path from axis to effect. **The quantity remains unidentified.** What we now know is that
it is not local DNA bending, not terminal bending, not end-to-end relative orientation, and not
footprint extent.

## Recommendation: do not spend the 24 arms

The plan gated 12 new pilots on this result. The honest reading:

- **Row 3 is the outcome** ("neither reproduces it") → *present as unresolved; do not spend.*
- The row-2 branch (interpretable, protein-side, worth 24 arms) **required** contact count to explain
  the correlation. It does not — the ×√L axis is untouched by it.
- Replication would confirm that a correlation exists. It would not say what the correlation is *of*,
  and after this run there is no remaining candidate for the quantity among the measurable DNA or
  interface properties. Spending weeks of compute to replicate an unidentified quantity is the wrong
  order, exactly as the plan's own prepared answer argues.

**What replaces it as the next question.** `bend_uu` summarises a duplex through terminal base-pair
axis *vectors*, and it disagrees with a complete local reconstruction at ρ = +0.05. That is a
measurement-definition problem, not a biology problem, and it is cheap to pin down: compare `bend_uu`
against the reconstruction on a synthetic duplex with known imposed curvature. If they diverge there
too, `bend_uu` is not measuring what its name implies and the correlation is an artifact of the
summary — which would close the mechanism line properly rather than leaving it open. **That is an hour,
and it is a better use of the slot than new pilots.**

## Effect on the talk

**Slide 6 (the retraction) gets stronger and simpler.** The current script says two orthogonal
reproducible axes disagree and it is unresolved at n = 10. Replace with: *every* reconstruction from
local DNA geometry — central window, both ends, end-to-end orientation, and the complete accumulation —
is flat against ΔP; only the whole-duplex terminal-vector summary correlates, and it tracks protein
contact count more strongly than it tracks DNA shape. The claim is withdrawn because the quantity
cannot be identified, and we enumerate what it is not.

**Slide 7 (F5) needs the same amendment flagged in the permutation report**, and it is now better
supported: the enumeration should say *no reconstructed DNA-geometry axis* predicts ΔP, with
`bend_uu` × √L named as the single unidentified exception.

**Prepared answer to update.** "Why not just run more pilots?" currently says 10–12 would suffice for
the whole-duplex axis. Keep the power arithmetic but add: we chose not to, because this run removed the
last candidate interpretation, and replicating an unidentified quantity does not advance the claim.

**Part 7's decision can be made now rather than after the talk, and the answer is no.** That is worth
saying out loud — the gate was defined in advance, it returned row 3, and the commitment was not made.

## Reproduce

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh && conda activate deeppbs
cd /scratch1/shewchuk/deepdnashape_seqshape
python extract_terminal.py      # -> terminal_bend.csv, all five accumulations per pilot/cond/W
```

Footprint inputs: `analysis/analyses/dna_relax/data/iface_mgwfl_vs_accuracy.csv` (`n_iface`, DNA levels
contacted) and `output/stage3_min/<tf>/<pdb>_fnat.csv` (`n_iface_res`, protein residues at interface).
Use W = 4 only; W = 6 leaves 2 pilots after the ΔP join. irf excluded throughout (no base-paired
duplex); lef1 has no own-family ΔP.
