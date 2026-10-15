# Addendum to `sox2_validation_guide_v3.md` — a correction I owe you, and a power problem

**2026-10-05.** Two things. First, **the UniProt bridge I gave you last turn was partly wrong and v3 has
it as a blocking assertion** — implemented as written, the assertion fails on correct data. Second, I
computed the Sox2–DNA interface from 1GT0 to check whether the test can see what it is looking for, and
it mostly cannot. §3.5's discriminator is inverted.

---

## 1. Correction: 9QPF is **+0**, not +30

My v2 addendum reported `9QPF auth 1–97, UniProt = auth +30`. That was a parsing error on my side — I
indexed `atom_site` by a hard-coded column position and read **`label_seq_id`** (which does run 1–97)
instead of `auth_seq_id`. Re-parsed from the mmCIF loop header:

| structure | auth span | **UniProt = auth +** | identity | covers (UniProt) |
|---|---|---|---|---|
| 9QPF chain A | **31–127** | **+0** | 96/97 | 31–127 |
| 2LE4 chain A | 1–81 | +37 | 80/81 | 38–118 |
| 1GT0 chain D | 1–79 | +38 | 79/79 | 39–117 |

**9QPF is already deposited in UniProt numbering.** Consequences for v3 §3.2:

- assertion 4 — "the UniProt bridge is **+30 / +37 / +38**" → **+0 / +37 / +38**
- assertion 2 — "pairwise offsets **+7** (2LE4→9QPF) and **+8** (1GT0.D→9QPF)" → auth-to-auth these are
  **+37** and **+38**. The +7/+8 figures are 2LE4/1GT0 auth against 9QPF's *label* numbering, which is
  not what an aligner working on deposited files will see.
- the helix table's 9QPF column is wrong throughout: **α3 is 9QPF auth 84–103**, not 54–73 (and
  identically for α1 47–62, α2 68–81 — 9QPF auth *is* UniProt).

**Unchanged and re-verified:** 2LE4 +37, 1GT0.D +38, the shared core at UniProt **39–117** (79
contiguous residues), and all four secondary-structure elements complete in all three structures. The
2LE4 column and the 1GT0.D column of v3's helix table are correct as printed.

The irony is that v3's own warning applies to my error: an offset check that passes on a
plausible-but-wrong frame. The identity-fraction assertion (#3) is what would have caught it.

## 2. The power problem: α3 carries 3% of the DNA interface

I computed the Sox2–DNA contact set from 1GT0 chain D directly (4.5 Å heavy-atom, waters excluded):
**25 Sox2 residues, 72 residue-level contact pairs.** By segment:

| segment | interface residues | contact pairs | share |
|---|---|---|---|
| N-term 39–46 | 6 | 20 | 28% |
| α1 47–62 | 4 | 11 | 15% |
| α2 68–81 | 7 | 15 | 21% |
| **α3 84–103** | **1** (R98) | **2** | **3%** |
| C-term 104–117 | 7 | 24 | 33% |

**The helix whose orientation distinguishes 2LE4 from the refined structure makes two DNA contacts.**
This inverts §3.5's discriminator:

> v3: *"Loss concentrated in α3 → the read-out measured the real difference. Loss diffuse → quality
> artifact, inconclusive."*

A genuine 30° α3 hinge **cannot** produce localized α3 loss — there are only 2 contact pairs there to
lose, 3% of fnat. What it produces is a **distorted global superposition**: α3 is well ordered in both
bundles (spread ≤ 1.5 Å in 2LE4, ≤ 2.0 Å in 9QPF) so it survives trimming and stays in the fit, and
misorienting it tilts the whole molecule in the common frame, displacing the α1/α2/N-term contacts that
carry the signal. **The true-positive signature is diffuse loss — which v3 reads as inconclusive.**

**Replacement discriminator.** Score 2LE4 twice: superposed on the whole HMG (47–103) and on the **major
wing only (α1+α2, 47–81)**.

| pattern | reading |
|---|---|
| fnat poor on whole-HMG fit, **recovers** on major-wing fit | a rigid-body hinge at α3 — the paper's finding, cleanly isolated |
| poor under both | distributed structural error; quality artifact not excluded |
| good under both | the read-out is insensitive to this rearrangement |

This is a better test than the per-residue decomposition and costs one extra superposition.

## 3. Half the interface is trimmed — and it is the half that differs

Per-residue CA spread across each 20-model bundle, superposed on the HMG core:

| segment | 2LE4 median | 9QPF median | above 3 Å |
|---|---|---|---|
| α1 47–62 | 0.65 Å | **0.27 Å** | 0 / 0 |
| α2 68–81 | 0.64 Å | **0.28 Å** | 0 / 0 |
| α3 84–103 | 0.94 Å | **0.35 Å** | 0 / 0 |
| N-term 39–46 | 2.80 Å | 3.37 Å | 4 / 4 |
| C-term 104–117 | **13.26 Å** | **14.81 Å** | 13 / 13 |

Applying the 3 Å cutoff, **both apo structures lose exactly the same 10 interface residues** — 40, 41,
42, 108, 109, 110, 112, 113, 114, 115 — which is **36 of the 72 contact pairs, exactly 50%**: the entire
C-terminal interface segment and half the N-terminal one.

What remains is dominated by α1 + α2 — **the region where 2LE4 and 9QPF agree**, both being tight,
well-ordered helices. The test trims the regions that differ and scores the region that doesn't.

**This raises the prior on v3 §4's pass/pass row from "a live possibility" to the modal expectation.**
The hinge can still register through the superposition mechanism above, but the direct interface signal
is largely gone. Worth knowing before committing 4.5 hours — and worth pre-registering as the expected
outcome rather than discovering.

## 4. Two things this makes easier

**The ternary confound (§6) is smaller than feared.** Only **three** Sox2 residues contact Oct-1 — 113,
114 and 116 — and the two that also touch DNA (113, 114) are both inside the trimmed C-terminal segment.
After trimming, the scored interface contains no Sox2–POU contact residues at all. The confound is real
for the *conformation* but essentially absent from the *contact set*.

**§7's scope statement should be sharper.** v3 says the tails are "deliberately out of scope" and the
test is "about the fold." But the C-terminal segment 104–117 — part of the tail the paper says rigidifies
on specific binding — carries **33% of the DNA contacts**. The trim is not discarding decoration; it is
discarding a third of the interface, and specifically the part whose ordering the paper identifies as
the binding-specific event. The honest statement is narrower:

> The read-out scores the **major-wing interface** (α1, α2 and the ordered half of the N-terminal arm) —
> about half of the native contact set. It is blind to the C-terminal arm, which contributes a third of
> the contacts and is where this protein's binding-coupled ordering happens.

That is still a defensible scope, and it sharpens the ets1/Sox2 pairing rather than weakening it.

## 5. Recommendation

The test is still worth running — it is half a day and the structures are verified — but with:

1. the corrected offsets (**+0 / +37 / +38**) in §3.2's assertions;
2. the **major-wing vs whole-HMG superposition** as the discriminator, replacing per-residue α3 localization;
3. **pass/pass pre-registered as the expected outcome**, with the scope statement in §4 above as the
   result rather than as a disappointment.

Read that way the test cannot fail to be informative: either the read-out detects an α3 hinge through a
distorted fit, or it demonstrates a quantified blind spot — half the interface, including the segment
that orders on binding.

## Reproduce

```bash
# all three numbers are from deposited files; files.rcsb.org is reachable from endeavour
# offsets: parse atom_site via the loop header and use auth_seq_id, NOT a fixed column index
# interface: 1GT0 chain D vs chains A+B, 4.5 A heavy atom, waters excluded -> 25 res / 72 pairs
# spread:   per-residue CA RMSF across the 20-model bundle, superposed on UniProt 47-103
```
