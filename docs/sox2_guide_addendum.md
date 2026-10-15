# Addendum to `sox2_validation_guide.md` — the three gating checks, resolved

**2026-10-05.** The guide flags §3.0 as the step that "gates everything" and leaves it open, and §7
lists two further risks that must be checked before the dock. All three are resolved below from the
deposited files. Two of them change the plan.

---

## 1. §3.0 — there is no binary Sox2–DNA structure. The ternary path is the only path.

RCSB returns 18 entries containing Sox2 with at least one DNA entity. Classified by protein-chain
composition:

| class | entries | usable |
|---|---|---|
| ternary, Sox2 + a POU partner | 1GT0, 1O4X, 6HT5, 8BX1, 8BX2 | yes, with chain surgery |
| nucleosomal (4–15 protein chains) | 6T78*, 6T79, 6T7A, 6T7C, 6T7D, 6T90, 6T93, 6YOV, 9RL4, 9RMC, 9RN1, 9RN2 | no |

**The guide's preferred option — "search for a binary Sox2–DNA structure, prefer it if one exists" —
does not exist.** The contingency in §3.0 step 2 is the only route, which means the irf failure mode is
not a risk to be avoided here but a condition to be handled. §3.0 step 3's verification is mandatory,
not precautionary.

**One trap worth naming.** `6T78` is X-ray, 2.50 Å, **one** protein entity plus two DNA chains — it
looks exactly like the binary hit you want, and a text search for "SOX2" returns it. It is
**Sox11**: *"Structure of human Sox11 transcription factor in complex with a short DNA fragment."* I
nearly recommended it. Check the title, not the composition.

**Recommended reference: 1GT0 chain D** (2.60 Å). It is the structure the paper itself aligns against,
so a disagreement is attributable to your pipeline rather than to a different choice of holo. 8BX1
chain A is marginally higher resolution (2.50 Å) but pairs Sox2 with Oct4 rather than Oct1.

**Sox2 is a different chain in every entry** — D in 1GT0 and 1O4X, C in 6HT5, A in 8BX1 and 8BX2.
Specify it explicitly; do not let the prep script pick "the first protein chain."

**A confound to state, not to fix.** In every available reference Sox2 is bound to DNA *and* to a POU
partner, and the POU/HMG interface is the point of those structures. So the comparison is strictly
"free Sox2 vs partner-and-DNA-bound Sox2." The paper has the same limitation — it also aligns to
1GT0 — so this does not invalidate the test, but it belongs in the caption.

## 2. §3.1 — the three structures use three different numbering frames

This is the check most likely to produce a silently wrong answer.

| structure | auth numbering | residues | offset to 9QPF |
|---|---|---|---|
| 9QPF (correct apo) | 1–97 | 97 | — |
| 2LE4 (misleading apo) | 1–81 | 81 | **+7** |
| 1GT0 chain D (holo) | 1–79 | 79 | **+8** |

At the offsets above the sequences are essentially perfect matches: 2LE4 → 9QPF gives **80 of 81**
identical, 1GT0.D → 9QPF gives **79 of 79**.

**By raw residue number, 2LE4 and 9QPF agree at only 7 of 81 positions.** An intersection taken on auth
numbering would superpose residue *i* of one structure onto a different residue of the other, and
return a large RMSD and a confident wrong answer with no error raised. The three-helix rearrangement
the test is designed to detect would be indistinguishable from the misalignment.

**Before the dock, assert that the core builder is aligning by sequence and not by residue number**,
and check the recovered offsets are +7 and +8.

## 3. §7 — the shared core is not a problem. It covers the whole HMG box.

The guide's largest listed risk is that trimming removes the α3 interface and the test becomes
impossible. It does not:

> **79 contiguous residues, span 9–87 in 9QPF numbering, no gaps.**
>
> `DRVKRPMNAFMVWSRGQRRKMAQENPKMHNSEISKRLGAEWKLLSETEKRPFIDEAKRLRALHMKEHPDYKYRPRRKTK`

That is the complete HMG box — the `RPMNAFMVW` signature at the N-terminus through the basic
`RPRRKTK` tail — so all three helices including α3 are inside the shared set, and the per-residue
decomposition in §3.3 has the positions it needs. Equivalent ranges: 2LE4 2–80, 1GT0 chain D 1–79.

The intersection is limited by the *holo* chain (79 residues), not by 2LE4. 9QPF's extra N- and
C-terminal residues — the tails carrying the paper's rigidity result — fall outside it, which is
consistent with §6: this test is about the fold, and the tails are deliberately out of scope.

## 4. Net effect on the plan

The guide's logic is sound and the test is executable. Three changes:

1. **Drop the binary-holo search**; go straight to 1GT0 chain D with explicit chain selection and the
   mandatory two-DNA-chain assertion. Saves the 30 min in §8 and removes a dead end.
2. **Add a sequence-alignment assertion** before the core build, checking the +7 / +8 offsets. This is
   new — it is not in the risk table and it is the failure most likely to produce a plausible wrong
   result.
3. **Close the §7 "core trims too much" risk** as checked and negative.

Cost is unchanged at roughly half a day. The BioEmu arm is cheaper than the guide assumes: on today's
batch, an ~80-residue construct at 100 samples runs in about 13 minutes wall, not an hour.

## Reproduce

```python
# chain composition of every Sox2-DNA entry, and the numbering offsets
# (files.rcsb.org is reachable from endeavour; data.rcsb.org is not allowlisted in the sandbox)
import gzip, urllib.request
raw = gzip.decompress(urllib.request.urlopen(
    "https://files.rcsb.org/download/1GT0.cif.gz").read()).decode()
# CA residues of chain D, first model -> compare with 2LE4 chain A and 9QPF chain A
```
