# Sox2 validation of the apo-fnat read-out — guide v4

**Supersedes** v1–v3. Revised 2026-10-05 after the interface and bundle-spread computations on 1GT0,
2LE4 and 9QPF.

**What this is.** A half-day test using a published NMR result to validate the Chapter 19½ mechanism
read-out against a case with both a **known right answer and a known wrong one**: can the read-out tell
a good apo structure from a misleading one?

**Prerequisite, new in v4:** run the **trim-survival audit** (`trim_survival_audit_spec.md`) first. It is
one hour, it covers the existing 14 pairs and 5 pilots, and it tells you what the read-out has been
measuring all along. Sox2's answer — 50% of the native contact set survives trimming — is the reason
this guide changed.

**What changed from v3.**

| v3 | v4 |
|---|---|
| UniProt bridge **+30** / +37 / +38 | **+0** / +37 / +38. 9QPF is deposited *in* UniProt numbering; the +30 came from reading `label_seq_id` instead of `auth_seq_id`. v3's 9QPF helix column was wrong throughout. |
| Pairwise offsets +7 / +8 | Auth-to-auth they are **+37 / +38**. The +7/+8 figures were 2LE4/1GT0 auth against 9QPF's *label* numbering — not what an aligner on deposited files sees. |
| Discriminator: is 2LE4's loss **localized to α3**? | **Inverted.** α3 carries 2 of 72 contact pairs (3%) — a hinge there *cannot* produce localized loss. Replaced by **major-wing vs whole-HMG superposition**. |
| pass/pass "a live possibility" | **The modal expectation**, pre-registered. The trim removes 50% of the contact set, and what remains is where the two apo structures agree. |
| Tails "deliberately out of scope" | The C-terminal arm carries **33% of the DNA contacts** and is where binding-coupled ordering happens. The scope statement is narrower and more honest. |
| Ternary confound needs a caption | **Smaller than feared** — after trimming, the scored interface contains no Sox2–POU contact residues at all. |

---

## 1 · The setup

Orsetti *et al.*, *Nucleic Acids Research* **53**(20), 2025 (doi:10.1093/nar/gkaf1121; PMC12585912).

| PDB | provenance | conformation | implied mechanism |
|---|---|---|---|
| **2LE4** | unpublished CESG deposit (2011); 20 conformers, selection by target function | *"the relative orientation of the three helices is significantly different compared to the DNA-bound state"* | **induced fit** — *"seemingly supported"* by it |
| **9QPF** | Orsetti *et al.* 2025; final ensemble of the 20 lowest-energy conformers | *"a three-helical bundle conformation identical to its DNA-bound state, indicating that folding is already encoded in the free state"* | **conformational selection** |

**Conclusion:** *"Binding of Sox2 to DNA thus involves conformational selection, rather than exclusively
induced fit, as was previously proposed."* — established by NOESY and RDC analysis, ¹⁵N relaxation,
thermal stability and MD. None of it interface geometry.

**Two further results that shape the test:**

- the HMG core has *"very similar rigidity in the free and bound states, while the N- and C-terminal
  tails become much more rigid upon binding"* → §7
- the folded free protein *"coexists in dynamic equilibrium with partially unfolded states, which are
  quenched upon specific DNA binding"* → §3.6

---

## 2 · The test

```
dock a free-state structure onto the holo DNA  →  score fnat against the holo native contact set
```

| source | expectation | role |
|---|---|---|
| 1GT0 chain D, DNA stripped, re-docked | ~0.8+ | positive control **and alignment canary** |
| **9QPF** (correct apo) | passes, and should land **near the holo control** | |
| **2LE4** (misleading apo) | **likely also passes** — see §3.4 | |
| BioEmu ensemble, split folded / partially unfolded | — | tests Finding 3, with the §3.6 caveat |

---

## 3 · Procedure

### 3.1 — The holo reference: 1GT0 chain D

**There is no binary Sox2–DNA structure.** Of 18 RCSB entries with Sox2 and a DNA entity: 5 ternary
(1GT0, 1O4X, 6HT5, 8BX1, 8BX2), 13 nucleosomal (4–15 protein chains, unusable). The irf failure mode is
**a condition to be handled, not a risk to be avoided.**

> ### ⚠ The 6T78 trap
> X-ray, 2.50 Å, one protein entity, two DNA chains — it looks exactly like the binary hit you want, and
> a text search for "SOX2" returns it. It is **Sox11**. Same class as the TBP `PWM_LABEL` error: a
> correct-looking identifier pointing at the wrong protein, invisible to any check on the file itself.
> **Assert the deposited title and UniProt ID**, not just chain composition.

**Use 1GT0 chain D** (2.60 Å) — the paper's own alignment target: *"Structural alignment of the
RDC-refined Sox2 31–127 structure (HMG residues 47–104) with an unpublished free state structure
(PDB: 2LE4) and the DNA-bound structure (PDB: 1GT0)."* A disagreement is then attributable to your
pipeline, not to a different choice of holo.

*Caveat:* 1GT0 is missing the N- and C-terminal tails, rebuilt by the authors with MODELLER. Consistent
with chain D covering UniProt 39–117.

**Sox2 is a different chain in every entry:** D in 1GT0 and 1O4X, C in 6HT5, A in 8BX1/8BX2. Specify it;
do not let the prep script take "the first protein chain."

**Mandatory before scoring:** both DNA chains present in `structures/source_chains/1gt0_chains/`.
**Restrict the native contact set to Sox2-chain contacts.**

### 3.2 — Numbering: four frames ⭐ **blocking**

| structure | auth span | **UniProt = auth +** | identity | covers (UniProt) |
|---|---|---|---|---|
| 9QPF chain A | **31–127** | **+0** | 96/97 | 31–127 |
| 2LE4 chain A | 1–81 | **+37** | 80/81 | 38–118 |
| 1GT0 chain D | 1–79 | **+38** | 79/79 | 39–117 |

Verified against the paper's hydrophobic-core landmarks (A47, F48, V50, W51, I71, L75, W79, F90, A94,
L97 — all ten match). **9QPF is already deposited in UniProt numbering.**

> **Parse `auth_seq_id` via the mmCIF loop header, never by fixed column index.** v3 reported +30 for
> 9QPF because it read `label_seq_id` (which does run 1–97). That is exactly the plausible-but-wrong
> frame this section warns about, and only the identity-fraction check catches it.

**Why a number-based intersection is the most dangerous step.** By raw auth number, 2LE4 and 9QPF agree
at almost no positions. A number-based intersection superposes residue *i* onto a different residue,
returns a large RMSD, and **raises no error** — looking exactly like the rearrangement the test is
designed to detect.

**Four blocking assertions before the core build:**

1. the core builder aligns **by sequence**, not residue number;
2. auth-to-auth offsets are **+37** (2LE4→9QPF) and **+38** (1GT0.D→9QPF);
3. identity fractions **≥ 80/81** and **≥ 79/79**;
4. the UniProt bridge is **+0 / +37 / +38**, and α3 (UniProt 84–103) maps into the shared core.

**Secondary structure** — *"The HMG domain spans from residue 47 to 103, with helix α1 spanning between
residue 47 and 62, α2 from residue 68 to 81 and α3 from residue 84 to 103"* — in each file's numbering:

| element | UniProt | 9QPF auth | 2LE4 auth | 1GT0.D auth |
|---|---|---|---|---|
| HMG domain | 47–103 | **47–103** | 10–66 | 9–65 |
| α1 | 47–62 | **47–62** | 10–25 | 9–24 |
| α2 | 68–81 | **68–81** | 31–44 | 30–43 |
| **α3** | **84–103** | **84–103** | **47–66** | **46–65** |

All elements complete in all three structures. **Report everything in UniProt numbering.**

### 3.3 — Shared core, and what it costs

> **79 contiguous residues, UniProt 39–117, no gaps.**
>
> `DRVKRPMNAFMVWSRGQRRKMAQENPKMHNSEISKRLGAEWKLLSETEKRPFIDEAKRLRALHMKEHPDYKYRPRRKTK`

The complete HMG box. The intersection is limited by the **holo** chain, not by 2LE4 — so trimming is
not discarding 2LE4's problematic regions.

**But the 3 Å disorder trim is not free.** Per-residue CA spread across each 20-model bundle:

| segment | 2LE4 median | 9QPF median | residues above 3 Å |
|---|---|---|---|
| α1 47–62 | 0.65 Å | 0.27 Å | 0 / 0 |
| α2 68–81 | 0.64 Å | 0.28 Å | 0 / 0 |
| α3 84–103 | 0.94 Å | 0.35 Å | 0 / 0 |
| N-term 39–46 | 2.80 Å | 3.37 Å | 4 / 4 |
| C-term 104–117 | **13.26 Å** | **14.81 Å** | 13 / 13 |

**Both apo structures lose the same 10 interface residues** — 40, 41, 42, 108, 109, 110, 112, 113, 114,
115 — which is **36 of 72 contact pairs, exactly 50%**. What survives is dominated by α1 + α2, **the
region where 2LE4 and 9QPF agree**.

The test trims the regions that differ and scores the region that doesn't. This is the basis for §3.4.

### 3.4 — The interface, and the inverted discriminator ⭐ **new in v4**

Sox2–DNA contacts from 1GT0 chain D (4.5 Å heavy-atom, waters excluded): **25 residues, 72 contact
pairs.**

| segment | interface residues | contact pairs | share |
|---|---|---|---|
| N-term 39–46 | 6 | 20 | 28% |
| α1 47–62 | 4 | 11 | 15% |
| α2 68–81 | 7 | 15 | 21% |
| **α3 84–103** | **1** (R98) | **2** | **3%** |
| C-term 104–117 | 7 | 24 | 33% |

**α3 — the helix whose orientation distinguishes 2LE4 from the refined structure — makes two DNA
contacts.**

This inverts v3's discriminator. v3 read *localized α3 loss* as the true positive and *diffuse loss* as a
quality artifact. But a genuine 30° α3 hinge **cannot** produce localized α3 loss: there are only 2
contact pairs to lose. α3 is well ordered in both bundles (0.94 Å, 0.35 Å), survives trimming, stays in
the fit — so misorienting it **tilts the whole molecule in the common frame**, displacing the α1/α2/
N-term contacts that carry the signal. **The true-positive signature is diffuse loss.**

**Replacement discriminator: score 2LE4 twice**, superposed on (a) the whole HMG (47–103) and (b) the
**major wing only, α1+α2 (47–81)**.

| pattern | reading |
|---|---|
| fnat poor on whole-HMG fit, **recovers** on major-wing fit | **a rigid-body hinge at α3** — the paper's finding, cleanly isolated. This is the result. |
| poor under both | distributed structural error; quality artifact not excluded |
| good under both | the read-out is insensitive to this rearrangement (see §3.5) |

One extra superposition, and a better test than per-residue localization.

### 3.5 — Pre-registered expectations

**Magnitudes from the paper:**

| quantity | value | use |
|---|---|---|
| 9QPF vs 1GT0, backbone RMSD over the HMG | **0.8 Å** | 9QPF should land **within noise of the holo-redocked control**, not merely above the gate |
| 2LE4 vs 9QPF, backbone RMSD (47–103) | **3.7 Å** | the deviation the dock must be sensitive to |
| α1–α3 interhelical angle | **107°** refined vs **77°** in 2LE4 | the ~30° hinge |
| RDC Q-factor | **0.08** refined, 0.35 NOESY-only, **0.76** for 2LE4 | independent confirmation 2LE4 is wrong |

**Pre-registered modal outcome: pass/pass.** The trim removes 50% of the contact set, and the surviving
half is where the two apo structures agree. α3 contributes 3% of contacts directly. The hinge can still
register through the §3.4 superposition mechanism, but **the direct interface signal is largely gone.**

**If pass/pass, that is the result, not a null.** Report the 50% figure and the scope statement in §7 as
the finding: the read-out scores half the native contact set, and the half it scores is the half that
doesn't change.

### 3.6 — Dock, score, and the BioEmu arm

Unmodified Stage-2 dock, `fnat_gate` at `FNAT_FLOOR` = 0.5.

**Two documented traps (Chapter 19½ README):** scoring the deposited holo *with* its own DNA is
near-circular — strip the DNA and re-dock. And an all-`NA` fnat file means the **reference chain
selection is wrong**, not a zero score; with a ternary complex this is the likeliest first failure.

**BioEmu:** one run on the Sox2 construct, ~13 min wall at ~80 residues / 100 samples.

The free-state ensemble genuinely contains partially unfolded members, so a generator sampling them is
behaving **correctly** while scoring low fnat for the right reason. "Below both apo references" is
ambiguous between *the generator is bad* and *the generator is right and fnat is measuring foldedness*.

**Discriminator:** split by whether the HMG core (47–103) is folded — the paper's MD shows it stable
across 15 µs — and score the folded subpopulation separately.

| result | reading |
|---|---|
| Folded frames reach 9QPF-like fnat; median dragged down by unfolded ones | **The generator reproduces a real equilibrium.** Reframes Finding 3. |
| Folded frames still well below 9QPF | The generator never reaches the bound-competent fold. **Finding 3 is substantially stronger.** |
| Lands nearer 2LE4 than 9QPF | **A named failure mode:** the generator reproduces a conformation the field has shown to be wrong. |

> **This exposure is not confined to Sox2.** Deposited apo bundles are biased toward folded states — NMR
> determination selects for a determinable fold — so across all five pilots, "best frame below the apo
> median" may partly measure that bias. **Run the folded-subpopulation split on the pilots too.**

---

## 4 · Reading the whole result

| 9QPF | 2LE4 whole-HMG | 2LE4 major-wing | reading |
|---|---|---|---|
| near holo control | poor | **recovers** | **The strongest validation available.** The hinge is isolated as a rigid-body rearrangement, recovered from interface geometry alone. |
| near holo control | poor | poor | Distributed error; quality artifact not excluded. Report with caveat. |
| pass | **pass** | pass | **The modal expectation.** A quantified blind spot: 50% of contacts trimmed, α3 only 3% of the rest. Report §7's scope statement as the result. |
| fail | — | — | Check §3.2's four assertions and the holo canary, then §3.1's chain handling. Do not interpret. |

---

## 5 · Cost

| step | time |
|---|---|
| Fetch; title + UniProt assertion; verify two DNA chains | 30 min |
| Alignment + UniProt bridge assertions (4 checks) | 30 min |
| Shared core across four sources | 45 min |
| Dock + `fnat_gate`, three deposited sources | 1 h |
| **Major-wing vs whole-HMG superposition** | 30 min |
| BioEmu run | ~13 min, unattended |
| BioEmu scored, folded/unfolded split | 45 min |
| *(optional)* 8BX1 replicate | 30 min |

**~4.5 hours hands-on.**

---

## 6 · The ternary confound — smaller than feared

Only **three** Sox2 residues contact Oct-1: 113, 114 and 116. The two that also touch DNA (113, 114) are
both inside the trimmed C-terminal segment. **After trimming, the scored interface contains no Sox2–POU
contact residues at all.**

The confound is real for the *conformation* — Sox2's fold in 1GT0 may be partly shaped by the partner
interface — but essentially absent from the *contact set*. State it that way.

**Replicate:** if the primary result is clean, **8BX1 chain A** is nearly free — different partner (Oct4
rather than Oct1), different crystal form, one more dock.

---

## 7 · What the test sees — the sharpened scope statement

The C-terminal segment 104–117 carries **33% of the DNA contacts** and is part of the tail the paper says
rigidifies on specific binding. The trim is not discarding decoration; it is discarding a third of the
interface, and specifically the part whose ordering the paper identifies as the binding-specific event.

> The read-out scores the **major-wing interface** — α1, α2 and the ordered half of the N-terminal arm,
> about half the native contact set. It is **blind to the C-terminal arm**, which contributes a third of
> the contacts and is where this protein's binding-coupled ordering happens.

**This sharpens the ets1/Sox2 pairing rather than weakening it.** Both are now the same statement:

> The read-out sees the **preformed core** and misses the **coupling region**. Ets-1's autoinhibitory
> module unfolds on DNA binding but is flanking and was trimmed; Sox2's C-terminal arm orders on binding
> and is trimmed. In both cases "selection-like" is correct *for the fold* and silent on the mechanism.

---

## 8 · What this buys

- **Fills a flagged gap.** Chapter 19½ notes the HMG-box split — NHP6A induced-fit, Sox2 conformational
  selection — and **lef1, the HMG-box pilot, has no apo structure.**
- **Upgrades Finding 3 from a statistic to a named case**, and via §3.6 may force a sharper version.
- **Validates a method rather than adding a datapoint** — and in the modal outcome, quantifies its blind
  spot, which is equally useful.

---

## 9 · Risk table

| risk | check | status |
|---|---|---|
| `label_seq_id` read as `auth_seq_id` | Parse via mmCIF loop header; identity-fraction assertion | **blocking, bit v3** |
| UniProt bridge wrong — α3 over wrong residues | §3.2 assertion 4: **+0 / +37 / +38** | **blocking** |
| Numbering frames mismatch | §3.2 assertions 1–3 | **blocking** |
| Wrong protein (6T78 = Sox11) | Assert title **and** UniProt ID | **blocking** |
| Ternary complex breaks chain extraction | Two DNA chains before scoring | **mandatory** |
| Native contact set includes POU contacts | Restrict to Sox2 chains | mandatory — but see §6 |
| All-`NA` fnat output | Wrong reference chain selection, not zero | — |
| Alignment silently broken | Holo-redocked arm must score ~0.8+ | canary |
| Diffuse loss misread as quality artifact | §3.4 major-wing superposition | **inverted in v4** |
| BioEmu low fnat from real partial unfolding | §3.6 folded-subpopulation split | — |
| Shared core trims α3 | — | closed, negative |

---

## Reference

Orsetti, A., Slejfer, J., Ha, S., Kevelam, D.I., Tekkelenburg, J., Van Duijn, T., Leppäkoski, A.,
Sedrakyan, A., Szilagyi, Á., Schellevis, R.D., Soufi, A., Cojocaru, V. & van Ingen, H.
**Solution structure of the Sox2 DNA-binding domain reveals conformational selection in DNA binding.**
*Nucleic Acids Research* **53**(20) (2025). doi:10.1093/nar/gkaf1121 (PMC12585912).

Structures: **2LE4** (CESG, 2011, unpublished), **9QPF** (Orsetti *et al.*), **1GT0** chain D (Oct-1 POU
+ Sox2 on DNA, 2.60 Å), **8BX1** chain A (Oct4 + Sox2, 2.50 Å — replicate).
