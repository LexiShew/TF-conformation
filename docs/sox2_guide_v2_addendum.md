# Addendum to `sox2_validation_guide_v2.md` — paper read, one new blocking frame

**2026-10-05.** I fetched and read Orsetti *et al.* (doi:10.1093/nar/gkaf1121, PMC12585912) rather than
working from the abstract. Three results: one of my own corrections is **withdrawn**, the guide gains
**quantitative expectations** it was missing, and there is **a fourth numbering frame** that §3.2 does
not cover and §3.4 cannot run without.

---

## 1. The 1GT0 attribution is correct — my retraction was wrong

Last turn I removed "the structure the paper itself aligns against" from both of my documents because I
had no evidence for it. The paper confirms it:

> *"Structural alignment of the RDC-refined Sox2 31–127 structure (HMG residues 47–104) with an
> unpublished free state structure (PDB: 2LE4) and the DNA-bound structure (PDB: 1GT0)."* — Fig. 2C
> caption

and the quantitative result:

> *"the DNA-bound structure of Sox2 HMG domain (PDB 1GT0) resembles the free state Sox2 31–137 closely,
> with heavy atom backbone RMSD of 0.8 Å."*

**v2 §3.1 and §5 are right as written; restore them.** My caution was procedurally correct at the time —
I had asserted it without a source — but the claim itself holds. One caveat the paper adds: 1GT0 is
**missing the N- and C-terminal tails**, which the authors had to build with MODELLER. That is
consistent with chain D covering only UniProt 39–117.

## 2. A fourth numbering frame — the one the paper speaks in

§3.2 tabulates three auth frames and asserts offsets +7 and +8. Correct, but incomplete: **the paper
reports everything in UniProt P48431 numbering**, and that is the frame in which α3 is defined. Without
the bridge, §3.4's "is the loss localized to α3?" cannot be evaluated — you would be decomposing over
the wrong residues.

Verified against the paper's own hydrophobic-core landmarks (A47, F48, V50, W51, I71, L75, W79, F90,
A94, L97 — all ten match):

| structure | auth span | **UniProt = auth +** | identity | covers (UniProt) |
|---|---|---|---|---|
| 9QPF chain A | 1–97 | **+30** | 96/97 | 31–127 |
| 2LE4 chain A | 1–81 | **+37** | 80/81 | 38–118 |
| 1GT0 chain D | 1–79 | **+38** | 79/79 | 39–117 |

Consistent with §3.2's pairwise offsets (37−30 = +7, 38−30 = +8), so that table is not wrong — it is
just not sufficient to place the paper's landmarks.

**Secondary-structure windows, from the paper** (*"The HMG domain spans from residue 47 to 103, with
helix α1 spanning between residue 47 and 62, α2 from residue 68 to 81 and α3 from residue 84 to 103"*),
translated into each file's own numbering. **All are fully present in all three structures** — checked
residue by residue:

| element | UniProt | 9QPF auth | 2LE4 auth | 1GT0.D auth |
|---|---|---|---|---|
| HMG domain | 47–103 | 17–73 | 10–66 | 9–65 |
| α1 | 47–62 | 17–32 | 10–25 | 9–24 |
| α2 | 68–81 | 38–51 | 31–44 | 30–43 |
| **α3** | **84–103** | **54–73** | **47–66** | **46–65** |

**Add to §3.2's blocking assertions:** a fourth check that the UniProt bridge is +30 / +37 / +38, and
that α3 maps into the shared core. The per-residue decomposition should be reported in **UniProt
numbering** so it can be compared to the paper directly.

This also strengthens §3.3: the shared core is UniProt 39–117, which contains the complete HMG domain
(47–103) with 8 residues of margin at the N-terminus and 14 at the C-terminus. The paper's own
comparison region for 2LE4 — residues 47–103 — sits entirely inside it.

## 3. Quantitative expectations the guide was missing

v2's reading table (§4) is qualitative — pass/fail. The paper supplies numbers that let you
pre-register magnitudes instead:

| quantity | value | use |
|---|---|---|
| 9QPF vs 1GT0, backbone RMSD over the HMG | **0.8 Å** | 9QPF should score *close to the holo-redocked control*, not merely above the gate |
| 2LE4 vs 9QPF, backbone RMSD, residues 47–103 | **3.7 Å** | the deviation the dock must be sensitive to |
| α1–α3 interhelical angle | **107°** (refined) vs **77°** (2LE4) | a 30° hinge — the specific defect §3.4 is testing for |
| RDC Q-factor | **0.08** (refined), **0.35** (NOESY-only), **0.76** (2LE4) | independent confirmation that 2LE4 is the wrong one |

The 0.8 Å figure is the most useful: it converts §4's top row from "9QPF passes" to "**9QPF scores
within noise of the holo control**", which is a much harder test to pass by accident and a much
stronger statement if it does.

Note also that a 30° hinge in one helix may or may not drop fnat below 0.5 after rigid-body docking.
The "pass / pass" row of §4 — the read-out being insensitive to this rearrangement — is a live
possibility, not a remote one. v2 already handles it correctly as a scope limitation.

## 4. One correction to the BioEmu arm's reading

§3.5 reads "lands below both apo references" as *"consistent with the pilots, extends Finding 3 to a
15th system."* The paper complicates that:

> the folded free protein *"coexists in dynamic equilibrium with partially unfolded states, which are
> quenched upon specific DNA binding."*

So the free-state ensemble genuinely contains partially unfolded members. A generator that samples them
is behaving **correctly**, and would score low fnat for the right reason. "Below both" is therefore
ambiguous between *the generator is bad* and *the generator is right and the metric is measuring
foldedness*.

**Discriminator:** split the BioEmu ensemble by whether the HMG core (UniProt 47–103) is folded — the
paper's own MD shows this region stable across 15 µs — and score the folded subpopulation separately.
If the folded frames reach 9QPF-like fnat and the overall median is dragged down by unfolded ones,
that is a very different result from a generator that never reaches the bound-competent fold. It also
makes the Sox2 arm a sharper test of Finding 3 than the pilots were.

## 5. Unverified in v2

The provenance line *"9QPF … 200 conformers calculated"* — the paper states the final ensemble
*"comprises the 20 lowest-energy conformers"*, and I did not find a 200-conformer statement in the
text. The deposited entry has 20 models. Either soften the claim or source it from the PDB entry's
deposition record.

## Reproduce

```bash
# full text
python -c "print(open('articles/10.1093_nar_gkaf1121.txt').read().count('1GT0'))"   # 7 mentions
# numbering bridge: align each chain's CA sequence to UniProt P48431 residues 31-127
# verified offsets +30 (9QPF.A), +37 (2LE4.A), +38 (1GT0.D)
```
