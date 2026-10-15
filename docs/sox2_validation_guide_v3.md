# Sox2 validation of the apo-fnat read-out — guide v3

**Supersedes** v1 and v2. Revised 2026-10-05 after reading Orsetti *et al.* in full rather than from the
abstract.

**What this is.** A half-day test that uses a published NMR result to validate the Chapter 19½ mechanism
read-out against a case with both a **known right answer and a known wrong one**. Most validation asks
"does the method agree with the literature." This asks: **can the read-out tell a good apo structure
from a misleading one?**

**What changed from v2.**

| v2 | v3 |
|---|---|
| Three auth numbering frames, offsets +7 / +8 | **A fourth frame.** The paper reports everything in UniProt P48431 numbering, and that is the frame α3 is defined in. Without the bridge (+30 / +37 / +38), §3.5's decomposition runs over the wrong residues. **Blocking.** |
| Reading table is qualitative (pass/fail) | **Pre-registered magnitudes.** 9QPF–1GT0 backbone RMSD 0.8 Å, 2LE4–9QPF 3.7 Å, interhelical angle 107° vs 77°, RDC Q 0.08 vs 0.76. |
| "1GT0 is the paper's own alignment target" — removed as unsourced | **Restored**, now sourced. |
| BioEmu "below both apo" = extends Finding 3 | **Ambiguous.** The free state genuinely contains partially unfolded members. Needs the folded-subpopulation split — and so does Finding 3 itself. |
| 9QPF "200 conformers calculated" | Softened — the paper describes 20 lowest-energy conformers; the deposited entry has 20 models. |

---

## 1 · The setup

Orsetti *et al.*, *Nucleic Acids Research* **53**(20), 2025 (doi:10.1093/nar/gkaf1121) determined the
solution structure of the free human Sox2 DNA-binding domain and overturned the prevailing induced-fit
model.

| PDB | provenance | conformation | implied mechanism |
|---|---|---|---|
| **2LE4** | unpublished CESG deposit (2011); 20 submitted conformers, selection by target function; 81 residues | *"the relative orientation of the three helices is significantly different compared to the DNA-bound state"* | **induced fit** — the model the paper says was *"seemingly supported"* by it |
| **9QPF** | Orsetti *et al.* 2025; solution NMR, final ensemble of the 20 lowest-energy conformers | *"a three-helical bundle conformation identical to its DNA-bound state, indicating that folding is already encoded in the free state"* | **conformational selection** |

**Conclusion:** *"Binding of Sox2 to DNA thus involves conformational selection, rather than exclusively
induced fit, as was previously proposed."*

Established by NOESY and RDC analysis, ¹⁵N backbone relaxation, thermal stability assays and MD — none
of which is interface geometry. If the apo-fnat read-out separates 9QPF from 2LE4, it recovered the same
conclusion by an independent route.

**Two further results from the paper that shape the test:**

- the HMG core has *"very similar rigidity in the free and bound states, while the N- and C-terminal
  tails become much more rigid upon binding"* → §7
- the folded free protein *"coexists in dynamic equilibrium with partially unfolded states, which are
  quenched upon specific DNA binding"* → §3.6, and it has consequences for Finding 3 generally

---

## 2 · The test

The Chapter 19½ Layer-2 procedure with one new pair and one extra arm:

```
dock a free-state structure onto the holo DNA  →  score fnat against the holo native contact set
```

| source | expectation | role |
|---|---|---|
| 1GT0 chain D, DNA stripped, re-docked | ~0.8+ | positive control **and alignment canary** |
| **9QPF** (correct apo) | **passes — and should land near the holo control** | |
| **2LE4** (misleading apo) | **fails** | |
| BioEmu ensemble, split folded / partially unfolded | — | tests Finding 3, with the §3.6 caveat |

---

## 3 · Procedure

### 3.1 — The holo reference: 1GT0 chain D

**There is no binary Sox2–DNA structure.** Of 18 RCSB entries containing Sox2 with a DNA entity:

| class | entries | usable |
|---|---|---|
| ternary, Sox2 + a POU partner | 1GT0, 1O4X, 6HT5, 8BX1, 8BX2 | yes, with chain surgery |
| nucleosomal (4–15 protein chains) | 6T79, 6T7A, 6T7C, 6T7D, 6T90, 6T93, 6YOV, 9RL4, 9RMC, 9RN1, 9RN2 | no |

So the irf failure mode is **a condition to be handled, not a risk to be avoided.**

> ### ⚠ The 6T78 trap
> `6T78` is X-ray, 2.50 Å, **one** protein entity plus two DNA chains — it looks exactly like the binary
> hit you want, and a text search for "SOX2" returns it. It is **Sox11**: *"Structure of human Sox11
> transcription factor in complex with a short DNA fragment."*
>
> Same class as the TBP `PWM_LABEL` error: a correct-looking identifier pointing at the wrong protein,
> invisible to any check on the file itself. **Assert the deposited title and UniProt ID match the
> intended protein**, not just chain composition.

**Use 1GT0 chain D** (2.60 Å). It is the structure the paper itself aligns against —
*"Structural alignment of the RDC-refined Sox2 31–127 structure (HMG residues 47–104) with an
unpublished free state structure (PDB: 2LE4) and the DNA-bound structure (PDB: 1GT0)"* — so a
disagreement is attributable to your pipeline rather than to a different choice of holo.

*Caveat from the paper:* 1GT0 is **missing the N- and C-terminal tails**, which the authors rebuilt with
MODELLER. Consistent with chain D covering only UniProt 39–117. Irrelevant to this test, which is about
the fold (§7), but worth knowing.

**Sox2 is a different chain in every entry:** D in 1GT0 and 1O4X, C in 6HT5, A in 8BX1 and 8BX2.
**Specify it explicitly; do not let the prep script pick "the first protein chain."**

**Mandatory before scoring:** verify both DNA chains survive into `structures/source_chains/1gt0_chains/`.
RCSB metadata will not predict this — 1IF1 deposits both strands and local prep still dropped one.

**Restrict the native contact set to Sox2-chain contacts only.**

### 3.2 — Numbering: four frames, not three ⭐ **blocking**

v2 tabulated three auth frames and asserted pairwise offsets +7 and +8. Correct but **insufficient**:
the paper reports everything in **UniProt P48431** numbering, and that is the frame in which α3 is
defined. Without the bridge, §3.5's decomposition runs over the wrong residues.

| structure | auth span | **UniProt = auth +** | identity | covers (UniProt) |
|---|---|---|---|---|
| 9QPF chain A | 1–97 | **+30** | 96/97 | 31–127 |
| 2LE4 chain A | 1–81 | **+37** | 80/81 | 38–118 |
| 1GT0 chain D | 1–79 | **+38** | 79/79 | 39–117 |

Verified against the paper's hydrophobic-core landmarks (A47, F48, V50, W51, I71, L75, W79, F90, A94,
L97 — all ten match). Consistent with the pairwise offsets: 37 − 30 = +7, 38 − 30 = +8.

**Why a number-based intersection is the most dangerous step in the procedure.** By raw residue number,
2LE4 and 9QPF agree at only **7 of 81** positions. A number-based intersection superposes residue *i* of
one structure onto a different residue of the other, returns a large RMSD, and **raises no error**. The
result looks exactly like the three-helix rearrangement the test is designed to detect — 2LE4 would
"fail," and you would conclude the method works. §3.5's decomposition would also be wrong, and could
look localized.

**Four blocking assertions, before the core build:**

1. the core builder aligns **by sequence**, not by residue number;
2. pairwise offsets are exactly **+7** (2LE4→9QPF) and **+8** (1GT0.D→9QPF);
3. identity fractions are **≥ 80/81** and **≥ 79/79**;
4. the **UniProt bridge** is **+30 / +37 / +38**, and α3 (UniProt 84–103) maps into the shared core.

An offset check alone passes if the aligner lands on a plausible-but-wrong frame. Identity fraction is
what makes it safe.

**Secondary structure, from the paper** — *"The HMG domain spans from residue 47 to 103, with helix α1
spanning between residue 47 and 62, α2 from residue 68 to 81 and α3 from residue 84 to 103"* — in each
file's own numbering. **All elements are fully present in all three structures:**

| element | UniProt | 9QPF auth | 2LE4 auth | 1GT0.D auth |
|---|---|---|---|---|
| HMG domain | 47–103 | 17–73 | 10–66 | 9–65 |
| α1 | 47–62 | 17–32 | 10–25 | 9–24 |
| α2 | 68–81 | 38–51 | 31–44 | 30–43 |
| **α3** | **84–103** | **54–73** | **47–66** | **46–65** |

**Report the per-residue decomposition in UniProt numbering** so it can be compared to the paper
directly.

### 3.3 — Shared core: no trimming problem

> **79 contiguous residues, UniProt 39–117, no gaps.**
>
> `DRVKRPMNAFMVWSRGQRRKMAQENPKMHNSEISKRLGAEWKLLSETEKRPFIDEAKRLRALHMKEHPDYKYRPRRKTK`

The complete HMG box — `RPMNAFMVW` signature through the basic `RPRRKTK` tail. The HMG domain (47–103)
sits inside with **8 residues of margin at the N-terminus and 14 at the C-terminus**, and the paper's own
comparison region for 2LE4 (47–103) is entirely contained.

Two consequences:

- The intersection is limited by the **holo** chain, **not** by 2LE4. Trimming is not discarding 2LE4's
  problematic regions — if 2LE4 fails, it fails on the part of the molecule actually being compared.
- 9QPF's extra terminal residues — the tails carrying the paper's rigidity result — fall outside.
  Consistent with §7: this test is about the fold.

### 3.4 — Pre-registered magnitudes ⭐ new in v3

The paper supplies numbers that convert §4's qualitative table into magnitude predictions:

| quantity | value | use |
|---|---|---|
| 9QPF vs 1GT0, backbone RMSD over the HMG | **0.8 Å** | 9QPF should land **within noise of the holo-redocked control**, not merely above the gate |
| 2LE4 vs 9QPF, backbone RMSD (47–103) | **3.7 Å** | the deviation the dock must be sensitive to |
| α1–α3 interhelical angle | **107°** refined vs **77°** in 2LE4 | the ~30° hinge §3.5 is testing for |
| RDC Q-factor | **0.08** refined, 0.35 NOESY-only, **0.76** for 2LE4 | independent confirmation that 2LE4 is the wrong structure |

**The 0.8 Å figure is the most useful.** It turns §4's top row from "9QPF passes" into "**9QPF scores
within noise of the holo control**" — much harder to pass by accident, and a much stronger statement if
it does.

**And a pre-registered caution:** a 30° hinge in one helix may or may not drop fnat below 0.5 after
rigid-body docking. **The pass/pass row of §4 is a live possibility, not a remote one.** It is a scope
limitation, not a failure, and it should be stated as an expected outcome rather than discovered as a
surprise.

### 3.5 — Dock, score, decompose

Run all sources through the unmodified Stage-2 dock and `fnat_gate` at `FNAT_FLOOR` = 0.5.

**Two documented traps from the Chapter 19½ README:**

- Scoring the deposited holo *with* its own DNA is near-circular — fnat's native contact set comes from
  that structure. **Strip the DNA and re-dock.**
- The gate writes `NA` when the reference has no native contacts (`fnat = rec/len(nat) if nat else nan`).
  An all-`NA` file means the **reference chain selection is wrong**, not a zero score. With a ternary
  complex this is the likeliest first failure.

**The per-residue decomposition is what makes this a test rather than a coincidence.** 2LE4 is a 2011
structural-genomics entry selected by target function; a lower-quality structure docks worse regardless
of conformation — a **quality artifact**, not a conformational finding.

**The discriminator:** is 2LE4's loss localized to **α3 (UniProt 84–103)**?

| pattern | reading |
|---|---|
| Loss concentrated in α3, rest scoring like 9QPF | The read-out measured the actual structural difference. **This is the result.** |
| Loss diffuse | Measuring structure quality. Inconclusive; report with the caveat. |

Run the same decomposition on 9QPF as the contrast. **Partial control if ambiguous:** check whether
2LE4's loss overlaps the POU interface or the DNA interface. If α3's loss is DNA-facing, the ternary
confound (§6) is less load-bearing.

### 3.6 — BioEmu arm, with the folded/unfolded split ⭐ revised in v3

One BioEmu run on the Sox2 construct — **~13 min wall** for ~80 residues at 100 samples.

**v2's reading of "lands below both apo references" was wrong.** The paper reports the folded free
protein *"coexists in dynamic equilibrium with partially unfolded states, which are quenched upon
specific DNA binding."* So the free-state ensemble genuinely contains partially unfolded members, and a
generator sampling them is behaving **correctly** while scoring low fnat for the right reason.

"Below both" is therefore ambiguous between *the generator is bad* and *the generator is right and fnat
is measuring foldedness*.

**Discriminator:** split the ensemble by whether the HMG core (UniProt 47–103) is folded — the paper's
MD shows this region stable across 15 µs — and score the folded subpopulation separately.

| result | reading |
|---|---|
| Folded frames reach 9QPF-like fnat; overall median dragged down by unfolded ones | **The generator is reproducing a real equilibrium.** A very different result from generator failure, and it reframes Finding 3. |
| Folded frames still well below 9QPF | The generator never reaches the bound-competent fold even when folded. **Finding 3 is substantially stronger** — the deficit is not a foldedness artifact. |
| Lands nearer 2LE4 than 9QPF | **A concrete, named failure mode:** the generator reproduces a conformation the field has since shown to be wrong. |

> ### This exposure is not confined to Sox2
> Finding 3 across all five pilots has the same ambiguity. Deposited apo bundles are biased toward
> folded states — NMR structure determination selects for a determinable fold — so if BioEmu samples a
> real partially-unfolded tail, "best frame below the apo median" partly measures that bias rather than
> generator failure.
>
> **Run the folded-subpopulation split on the five pilots too.** If folded BioEmu frames reach apo-like
> fnat, Finding 3 needs restating; if they don't, it gets substantially stronger. Worth settling before
> the generator-benchmark paper is written.

---

## 4 · Reading the whole result

| 9QPF | 2LE4 | α3 localized | reading |
|---|---|---|---|
| **near the holo control** | fail | yes | **The strongest validation available.** Independently recovers a published NMR conclusion *and* correctly rejects the misleading deposit, from interface geometry alone. |
| passes but well below holo | fail | yes | Still a success; note the gap to the 0.8 Å expectation. |
| pass | fail | no | Suggestive, confounded with structure quality. Report with the caveat. |
| pass | pass | — | **Expected possibility, not a surprise** (§3.4). The read-out is insensitive to a 30° single-helix hinge after rigid-body docking. A real scope limitation that bounds what the mechanism split can resolve. |
| fail | fail | — | Check §3.2's four assertions and the holo canary, then §3.1's chain handling. Do not interpret. |
| fail | pass | — | Requires explanation before any other conclusion. |

---

## 5 · Cost

| step | time |
|---|---|
| Fetch 1GT0 / 2LE4 / 9QPF; title + UniProt assertion; verify two DNA chains | 30 min |
| **Sequence alignment + UniProt bridge assertions (4 checks)** | 30 min |
| Shared core across four sources | 45 min |
| Dock + `fnat_gate`, three deposited sources | 1 h |
| Per-residue decomposition, reported in UniProt numbering | 1 h |
| BioEmu run | ~13 min, unattended |
| BioEmu scored, with folded/unfolded split | 45 min |
| *(optional)* 8BX1 replicate | 30 min |

**~4.5 hours hands-on.**

---

## 6 · The ternary confound, and the replicate

In **every** available reference, Sox2 is bound to DNA *and* to a POU partner, and the POU/HMG interface
is the point of those structures. The native contact set is therefore defined on a Sox2 conformation
possibly shaped in part by the protein–protein interface.

The paper has the same limitation — it also aligns to 1GT0 — so this does not invalidate the test, but
both measurements share the confound rather than it disappearing. **State it as:** *free Sox2 versus
Sox2 in a POU-partnered DNA complex; no binary reference exists for either analysis.*

**The replicate.** If the primary result is clean, **8BX1 chain A** is nearly free: different partner
(Oct4 rather than Oct1), different crystal form, same procedure, one more dock. Two independent holo
references agreeing would strengthen the claim considerably.

---

## 7 · What the test does and doesn't see

The paper reports the HMG **core** has similar rigidity free and bound, while the **N- and C-terminal
tails** rigidify on binding. The read-out trims disordered termini by design.

Chapter 19½ flags this with **ets1** as the case where trimming hid genuine induced fit: Ets-1's
autoinhibitory module unfolds on DNA binding, but it is flanking and was trimmed, so the 103 residues
scored were the preformed ETS domain.

**Sox2 is the mirror image:** core preformed, tails ordered on binding. "Selection-like for the fold" is
the *correct* answer here, and the literature agrees.

> The read-out answers **"is the DBD fold preformed?"** — reliably. It does **not** answer "does this
> protein undergo induced fit," because the regions where induced fit most often lives are the ones it
> trims. Ets-1 is where that distinction bites; Sox2 is where it doesn't.

**Report them as a pair.** It converts a limitations paragraph into a characterization of scope.

---

## 8 · What this buys beyond validation

- **Fills a flagged gap.** Chapter 19½ notes the HMG-box split — NHP6A induced-fit, Sox2 conformational
  selection — and **lef1, the HMG-box pilot, has no apo structure.** Sox2 gives the family an anchor. The
  paper itself notes that *"analysis of Sox2 orthologue structures in free and bound state presents
  conflicting evidence"* for the structural-change model — exactly the disagreement a geometric read-out
  can arbitrate.
- **Upgrades Finding 3 from a statistic to a named case** — and, via §3.6, may force a sharper version
  of it.
- **Validates a method rather than adding a datapoint** — which is why it ranks above the 2euv pilot.

---

## 9 · Risk table

| risk | check | status |
|---|---|---|
| **UniProt bridge missing — α3 decomposed over wrong residues** | §3.2 assertion 4: +30 / +37 / +38 | **blocking, new in v3** |
| Numbering frames mismatch | §3.2 assertions 1–3 | **blocking** |
| Wrong protein (6T78 = Sox11) | Assert title **and** UniProt ID | **blocking** |
| Ternary complex breaks chain extraction (irf mode) | Two DNA chains in `source_chains/` before scoring | **mandatory** |
| Native contact set includes POU contacts | Restrict to Sox2-chain contacts; inspect the list | mandatory |
| All-`NA` fnat output | Wrong reference chain selection, not a zero score | — |
| Alignment silently broken | Holo-redocked arm must score ~0.8+ | canary |
| 2LE4 fails on quality, not conformation | §3.5 per-residue decomposition. **Do not skip** | — |
| BioEmu low fnat from real partial unfolding | §3.6 folded-subpopulation split | **new in v3** |
| Shared core trims α3 | — | **closed, negative** |
| Construct differences 2LE4 vs 9QPF | — | **closed** |

---

## Reference

Orsetti, A., Slejfer, J., Ha, S., Kevelam, D.I., Tekkelenburg, J., Van Duijn, T., Leppäkoski, A.,
Sedrakyan, A., Szilagyi, Á., Schellevis, R.D., Soufi, A., Cojocaru, V. & van Ingen, H.
**Solution structure of the Sox2 DNA-binding domain reveals conformational selection in DNA binding.**
*Nucleic Acids Research* **53**(20) (2025). doi:10.1093/nar/gkaf1121. (PMC12585912)

Key quotes:
- *"the relative orientation of the three helices is significantly different compared to the DNA-bound
  state"* — on 2LE4
- *"a three-helical bundle conformation identical to its DNA-bound state, indicating that folding is
  already encoded in the free state"*
- *"Binding of Sox2 to DNA thus involves conformational selection, rather than exclusively induced fit,
  as was previously proposed."*
- *"the HMG has very similar rigidity in the free and bound states, while the N- and C-terminal tails
  become much more rigid upon binding"*
- *"coexists in dynamic equilibrium with partially unfolded states, which are quenched upon specific DNA
  binding"*
- *"the DNA-bound structure of Sox2 HMG domain (PDB 1GT0) resembles the free state Sox2 31–137 closely,
  with heavy atom backbone RMSD of 0.8 Å"*
- *"The HMG domain spans from residue 47 to 103, with helix α1 spanning between residue 47 and 62, α2
  from residue 68 to 81 and α3 from residue 84 to 103"*

Structures: **2LE4** (Sahu, Markley, Tonelli, Bahrami, Eghbalnia; CESG; 2011, unpublished), **9QPF**
(Orsetti *et al.*), **1GT0** chain D (Oct-1 POU + Sox2 on DNA, 2.60 Å), **8BX1** chain A (Oct4 + Sox2,
2.50 Å — replicate).
