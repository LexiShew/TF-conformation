# Sox2 validation of the apo-fnat read-out — guide v2

**Supersedes** `sox2_validation_guide.md`. Revised 2026-10-05 against the pre-flight addendum, which
resolved all three gating checks from deposited files. Two of them changed the plan.

**What this is.** A half-day test that uses a published NMR result to validate the Chapter 19½ mechanism
read-out against a case with both a **known right answer and a known wrong one**. Most validation asks
"does the method agree with the literature." This asks: **can the read-out tell a good apo structure
from a misleading one?**

**What changed from v1.**

| v1 | v2 |
|---|---|
| "Search for a binary Sox2–DNA structure, prefer it if one exists" | **No binary structure exists.** 18 Sox2–DNA entries: 5 ternary, 13 nucleosomal. Go straight to 1GT0 chain D. |
| Risk table did not include numbering frames | **New blocking assertion (§3.2).** Three structures, three frames, offsets +7 and +8. A number-based intersection returns a confident wrong answer with no error raised — and it looks exactly like the result the test is designed to find. |
| "Shared core may trim the α3 interface" listed as the largest risk | **Closed, negative.** 79 contiguous residues covering the complete HMG box, all three helices in scope. |
| BioEmu ~1 GPU-hour | ~13 min wall for an 80-residue construct at 100 samples. |

---

## 1 · The setup

Orsetti *et al.*, *Nucleic Acids Research* **53**(20), 2025 (doi:10.1093/nar/gkaf1121) determined the
solution structure of the free human Sox2 DNA-binding domain and overturned the prevailing induced-fit
model.

| PDB | provenance | conformation | implied mechanism |
|---|---|---|---|
| **2LE4** | unpublished CESG deposit (2011); 100 conformers calculated, 20 submitted, selection by target function; 81 residues | *"the relative orientation of the three helices is significantly different compared to the DNA-bound state"* | **induced fit** — the model the paper says was *"seemingly supported"* by it |
| **9QPF** | Orsetti *et al.* 2025; solution NMR, 200 conformers calculated, lowest-energy selection | *"a three-helical bundle conformation identical to its DNA-bound state, indicating that folding is already encoded in the free state"* | **conformational selection** |

**Conclusion:** *"Binding of Sox2 to DNA thus involves conformational selection, rather than exclusively
induced fit, as was previously proposed."*

Established by NOESY and RDC analysis, ¹⁵N backbone relaxation, thermal stability assays and MD — none
of which is interface geometry. If the apo-fnat read-out separates 9QPF from 2LE4, it recovered the same
conclusion by an independent route.

**The dynamics result that matters for §6:** the HMG core has *"very similar rigidity in the free and
bound states, while the N- and C-terminal tails become much more rigid upon binding."*

---

## 2 · The test

The Chapter 19½ Layer-2 procedure with one new pair and one extra arm:

```
dock a free-state structure onto the holo DNA  →  score fnat against the holo native contact set
```

Four sources through the **same unmodified dock**, all trimmed to one shared residue set:

| source | expectation | role |
|---|---|---|
| 1GT0 chain D, DNA stripped, re-docked | passes (~0.8+) | positive control **and alignment canary** (§3.2) |
| **9QPF** (correct apo) | **passes** | |
| **2LE4** (misleading apo) | **fails** | |
| BioEmu ensemble of the Sox2 construct | — | tests Finding 3 against two apo references |

---

## 3 · Procedure

### 3.1 — The holo reference: 1GT0 chain D

**There is no binary Sox2–DNA structure.** Of 18 RCSB entries containing Sox2 with a DNA entity:

| class | entries | usable |
|---|---|---|
| ternary, Sox2 + a POU partner | 1GT0, 1O4X, 6HT5, 8BX1, 8BX2 | yes, with chain surgery |
| nucleosomal (4–15 protein chains) | 6T79, 6T7A, 6T7C, 6T7D, 6T90, 6T93, 6YOV, 9RL4, 9RMC, 9RN1, 9RN2 | no |

So the irf failure mode is **a condition to be handled, not a risk to be avoided.** The two-DNA-chain
verification below is mandatory, not precautionary.

> ### ⚠ The 6T78 trap
> `6T78` is X-ray, 2.50 Å, **one** protein entity plus two DNA chains — it looks exactly like the binary
> hit you want, and a text search for "SOX2" returns it. It is **Sox11**: *"Structure of human Sox11
> transcription factor in complex with a short DNA fragment."*
>
> This is the same class of error as the TBP `PWM_LABEL` mistake: a correct-looking identifier pointing
> at the wrong protein, invisible to any check on the file itself. **Assert the deposited title and
> UniProt ID match the intended protein**, not just the chain composition.

**Use 1GT0 chain D** (2.60 Å). It is the structure the paper itself aligns against, so a disagreement is
attributable to your pipeline rather than to a different choice of holo. 8BX1 chain A is marginally
higher resolution (2.50 Å) but pairs Sox2 with Oct4 rather than Oct1 — keep it as the replicate (§5).

**Sox2 is a different chain in every entry:** D in 1GT0 and 1O4X, C in 6HT5, A in 8BX1 and 8BX2.
**Specify it explicitly; do not let the prep script pick "the first protein chain."**

**Mandatory before scoring:** verify both DNA chains survive into `structures/source_chains/1gt0_chains/`.
RCSB metadata will not predict this — 1IF1 deposits both strands and local prep still dropped one, and
that went undetected for months.

**Restrict the native contact set to Sox2-chain contacts only.** Treat the POU chain as absent for
scoring.

### 3.2 — Sequence alignment ⭐ **blocking assertion; new in v2**

The three structures use three different auth numbering frames:

| structure | auth numbering | residues | offset to 9QPF | identity at that offset |
|---|---|---|---|---|
| 9QPF (correct apo) | 1–97 | 97 | — | — |
| 2LE4 (misleading apo) | 1–81 | 81 | **+7** | 80 / 81 |
| 1GT0 chain D (holo) | 1–79 | 79 | **+8** | 79 / 79 |

**By raw residue number, 2LE4 and 9QPF agree at only 7 of 81 positions.**

**Why this is the most dangerous step in the procedure.** A number-based intersection superposes residue
*i* of one structure onto a different residue of the other, returns a large RMSD, and **raises no
error**. The result looks exactly like the three-helix rearrangement the test is designed to detect —
2LE4 would "fail," and you would conclude the method works. Worse, §3.4's per-residue decomposition
would also be wrong, and could look localized. **The safeguard built into v1 does not catch this.**

**Assert all three, before the core build:**

1. the core builder aligns **by sequence**, not by residue number;
2. recovered offsets are exactly **+7** (2LE4) and **+8** (1GT0.D);
3. identity fractions are **≥ 80/81** and **≥ 79/79**.

An offset check alone passes if the aligner lands on a plausible-but-wrong frame. The identity fraction
is what makes it safe.

**Then one post-hoc canary.** If the holo-redocked arm does **not** score ~0.8+, the alignment is
suspect — a structure that *is* the bound state must survive the rigid-body round trip. That is the
cheapest independent signal that something upstream broke.

### 3.3 — Shared core: no trimming problem

v1's largest listed risk is closed, negative:

> **79 contiguous residues, span 9–87 in 9QPF numbering, no gaps.**
>
> `DRVKRPMNAFMVWSRGQRRKMAQENPKMHNSEISKRLGAEWKLLSETEKRPFIDEAKRLRALHMKEHPDYKYRPRRKTK`

That is the complete HMG box — the `RPMNAFMVW` signature through the basic `RPRRKTK` tail — so **all
three helices including α3 are inside the shared set** and §3.4 has the positions it needs. Equivalent
ranges: 2LE4 2–80, 1GT0 chain D 1–79.

Two consequences worth stating:

- The intersection is limited by the **holo** chain (79 residues), **not** by 2LE4. So trimming is not
  discarding 2LE4's problematic regions — if 2LE4 fails, it fails on the part of the molecule actually
  being compared. This removes an objection v1 left open.
- 9QPF's extra N- and C-terminal residues — the tails carrying the paper's rigidity result — fall
  outside the intersection. Consistent with §6: this test is about the fold, and the tails are
  deliberately out of scope.

### 3.4 — Dock, score, and decompose per residue

Run all four sources through the unmodified Stage-2 dock and `fnat_gate` at `FNAT_FLOOR` = 0.5. Report
median fnat and pass rate per source.

**Two documented traps from the Chapter 19½ README:**

- Scoring the deposited holo *with* its own DNA is near-circular — fnat's native contact set comes from
  that same structure. **Strip the DNA and re-dock.**
- The gate writes `NA` when the reference has no native contacts (`fnat = rec/len(nat) if nat else nan`).
  An all-`NA` file means the **reference chain selection is wrong**, not that the ensemble scored zero.
  With a ternary complex this is the likeliest first failure.

**The per-residue decomposition is what makes this a test rather than a coincidence.** 2LE4 is a 2011
structural-genomics entry selected by target function, and a lower-quality structure docks worse
regardless of conformation — a **quality artifact**, not a conformational finding.

**The discriminator:** is 2LE4's loss **localized to helix α3**, the specific rearrangement the paper
identifies?

| pattern | reading |
|---|---|
| Loss concentrated in α3, rest of the interface scoring like 9QPF | The read-out measured the actual structural difference. **This is the result.** |
| Loss diffuse across the interface | Measuring structure quality, not conformation. Test is inconclusive; report with the caveat. |

Run the same decomposition on 9QPF as the contrast.

**Partial control if the result is ambiguous:** check whether 2LE4's per-residue loss overlaps the POU
interface or the DNA interface. If α3's loss is on the DNA-facing side, the ternary confound (§5) is
less load-bearing.

### 3.5 — BioEmu arm

One BioEmu run on the Sox2 construct — **~13 min wall** for ~80 residues at 100 samples — scored through
the same dock.

The Finding 3 question: **does BioEmu's best frame beat the deposited apo median**, now with two apo
references at different quality?

| result | reading |
|---|---|
| Lands nearer 9QPF | Better than the 5-of-5 pilot result predicts. Worth explaining. |
| Lands nearer 2LE4 | **A concrete, named failure mode** rather than an aggregate statistic: the generator reproduces a conformation the field has since shown to be wrong. Strong and quotable. |
| Below both | Consistent with the pilots. Extends Finding 3 to a 15th system. |

---

## 4 · Reading the whole result

| 9QPF | 2LE4 | α3 localized | reading |
|---|---|---|---|
| pass | fail | yes | **The strongest validation available.** The read-out independently recovers a published NMR conclusion *and* correctly rejects the misleading deposit, from interface geometry alone. |
| pass | fail | no | Suggestive, confounded with structure quality. Report with the caveat. |
| pass | pass | — | The read-out is insensitive to the α3 rearrangement. **A real scope limitation** — worth reporting, since it bounds what the mechanism split can resolve. |
| fail | fail | — | Check §3.2's assertions and the holo canary first, then §3.1's chain handling. Do not interpret. |
| fail | pass | — | Requires explanation before any other conclusion. Dock and contact set first. |

---

## 5 · The ternary confound, and the replicate

In **every** available reference, Sox2 is bound to DNA *and* to a POU partner, and the POU/HMG interface
is the point of those structures. So the native contact set is defined on a Sox2 conformation that may
be partly shaped by the protein–protein interface, not by DNA alone.

The paper has the same limitation — it also aligns to 1GT0 — so this does not invalidate the test. But
both measurements share the confound rather than it disappearing. **State it as:** *free Sox2 versus
Sox2 in a POU-partnered DNA complex; no binary reference exists for either analysis.*

**The replicate.** If the primary result is clean, **8BX1 chain A** is nearly free: different partner
(Oct4 rather than Oct1), different crystal form, same procedure, one more dock. Two independent holo
references agreeing would strengthen the claim considerably. Queue it.

---

## 6 · What the test does and doesn't see

The paper reports the HMG **core** has similar rigidity free and bound, while the **N- and C-terminal
tails** rigidify on binding. The read-out trims disordered termini by design.

Chapter 19½ already flags this with **ets1** as the case where trimming hid genuine induced fit: Ets-1's
autoinhibitory module unfolds on DNA binding, but it is flanking and was trimmed, so the 103 residues
scored were the preformed ETS domain.

**Sox2 is the mirror image:** core preformed, tails ordered on binding. "Selection-like for the fold" is
the *correct* answer here, and the literature agrees.

> The read-out answers **"is the DBD fold preformed?"** — reliably. It does **not** answer "does this
> protein undergo induced fit," because the regions where induced fit most often lives are the ones it
> trims. Ets-1 is where that distinction bites; Sox2 is where it doesn't.

**Report them as a pair.** It converts a limitations paragraph into a characterization of the method's
scope.

---

## 7 · What this buys beyond validation

- **It fills a gap already flagged.** Chapter 19½'s limits section notes the HMG-box split — NHP6A
  induced-fit, Sox2 reported as conformational selection — and **lef1, the project's HMG-box pilot, has
  no apo structure.** Sox2 gives the family an anchor. The paper itself notes that *"analysis of Sox2
  orthologue structures in free and bound state presents conflicting evidence"* for the structural-change
  model, which is exactly the disagreement a geometric read-out can arbitrate.
- **It upgrades Finding 3 from a statistic to a named case.** "BioEmu's best frame is below the apo
  median in 5 of 5 pilots" is strong; "and on Sox2 it reproduces the deposited structure the field has
  since shown to be wrong" is memorable.
- **It validates a method rather than adding a datapoint** — which is why it ranks above the 2euv NDT80
  pilot.

---

## 8 · Risk table

| risk | check | status |
|---|---|---|
| **Numbering frames mismatch** | §3.2 — sequence alignment, offsets +7/+8, identity ≥ 80/81 and 79/79 | **blocking, new in v2** |
| Wrong protein (6T78 = Sox11) | Assert deposited title **and** UniProt ID match the intended protein | **blocking** |
| Ternary complex breaks chain extraction (irf mode) | Two DNA chains present in `source_chains/` before scoring; assert at extraction | **mandatory** |
| Native contact set includes POU contacts | Restrict to Sox2-chain contacts; inspect the contact list | mandatory |
| All-`NA` fnat output | Wrong reference chain selection, not a zero score. Scripts should raise `BrokenReference` | — |
| Alignment silently broken | Holo-redocked arm must score ~0.8+ | canary |
| 2LE4 fails on quality, not conformation | §3.4 per-residue decomposition. **Do not skip** | — |
| Shared core trims α3 | — | **closed, negative** (§3.3) |
| Construct differences 2LE4 vs 9QPF | — | **closed** — intersection covers the full HMG box |

---

## 9 · Cost

| step | time |
|---|---|
| Fetch 1GT0 / 2LE4 / 9QPF; title + UniProt assertion; verify two DNA chains | 30 min |
| **Sequence-alignment assertion** | 20 min |
| Shared core across four sources | 45 min |
| Dock + `fnat_gate`, three deposited sources | 1 h |
| Per-residue decomposition | 1 h |
| BioEmu run | ~13 min, unattended |
| BioEmu scored | 30 min |
| *(optional)* 8BX1 replicate | 30 min |

**~4 hours hands-on.** It is `fnat_on_pairs.py` with a new row, one extra arm, and one new assertion.

---

## Reference

Orsetti, A., Slejfer, J., Ha, S., Kevelam, D.I., Tekkelenburg, J., Van Duijn, T., Leppäkoski, A.,
Sedrakyan, A., Szilagyi, Á., Schellevis, R.D., Soufi, A., Cojocaru, V. & van Ingen, H.
**Solution structure of the Sox2 DNA-binding domain reveals conformational selection in DNA binding.**
*Nucleic Acids Research* **53**(20) (2025). doi:10.1093/nar/gkaf1121.

Key quotes:
- *"the relative orientation of the three helices is significantly different compared to the DNA-bound
  state"* — on 2LE4
- *"a three-helical bundle conformation identical to its DNA-bound state, indicating that folding is
  already encoded in the free state"* — on the new structure
- *"Binding of Sox2 to DNA thus involves conformational selection, rather than exclusively induced fit,
  as was previously proposed."*
- *"the HMG has very similar rigidity in the free and bound states, while the N- and C-terminal tails
  become much more rigid upon binding"*

Structures: **2LE4** (Sahu, Markley, Tonelli, Bahrami, Eghbalnia; CESG; deposited 2011, unpublished),
**9QPF** (Orsetti *et al.*, solution NMR, 200 conformers calculated), **1GT0** chain D (Oct-1 POU +
Sox2 on DNA, 2.60 Å), **8BX1** chain A (Oct4 + Sox2, 2.50 Å — replicate).
