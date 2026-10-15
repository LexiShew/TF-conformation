# Sox2 validation of the apo-fnat read-out — a guide

**Written 2026-09-29.** A half-day, one-GPU-hour test that uses a published NMR result to validate the
Chapter 19½ mechanism read-out against a case with both a **known right answer and a known wrong one**.

**Why it is worth doing.** Most validation asks "does the method agree with the literature." This asks
something harder: **can the read-out tell a good apo structure from a misleading one?** Sox2 is the rare
system where two deposited free-state structures of the same domain disagree, and an orthogonal
measurement has settled which is right.

**When.** After the talk. Ahead of the 2euv NDT80 pilot in priority — it is cheaper and it validates a
method rather than adding a datapoint.

---

## 1 · The setup

Orsetti *et al.*, *Nucleic Acids Research* **53**(20), 2025 (doi:10.1093/nar/gkaf1121) determined the
solution structure of the free human Sox2 DNA-binding domain and overturned the prevailing induced-fit
model.

**The two free-state structures:**

| PDB | provenance | conformation | what it implied |
|---|---|---|---|
| **2LE4** | unpublished CESG structural-genomics deposit (2011); 100 conformers calculated, 20 submitted; selection by target function; 81 residues | "the relative orientation of the three helices is significantly different compared to the DNA-bound state" | **induced fit** — and the paper says this model was *"seemingly supported"* by it |
| **9QPF** | Orsetti *et al.* 2025; solution NMR, 200 conformers calculated, lowest-energy selection | "a three-helical bundle conformation identical to its DNA-bound state, indicating that folding is already encoded in the free state" | **conformational selection** |

**The paper's conclusion:** *"Binding of Sox2 to DNA thus involves conformational selection, rather than
exclusively induced fit, as was previously proposed."*

Established by NOESY and RDC analysis, ¹⁵N backbone relaxation across the full domain, thermal stability
assays, and MD — none of which is interface geometry. So if the apo-fnat read-out separates 9QPF from
2LE4, it recovered the same conclusion by an independent route.

**The additional dynamics result, which matters for §6:** the HMG core has very similar rigidity in free
and bound states, while the N- and C-terminal tails become much more rigid upon binding.

---

## 2 · What the test is

Exactly the Chapter 19½ Layer-2 procedure, with one new pair and one extra arm:

```
dock a free-state structure onto the holo DNA  →  score fnat against the holo native contact set
```

Four sources through the **same unmodified dock**, all trimmed to one shared residue set:

| source | expectation |
|---|---|
| holo Sox2, DNA stripped, re-docked | passes (the positive control — proves the dock is lossless) |
| **9QPF** (correct apo) | **passes** |
| **2LE4** (misleading apo) | **fails** |
| BioEmu ensemble of the Sox2 construct | tests Finding 3 against two apo references |

---

## 3 · Procedure

### 3.0 — Find the holo structure (do this first; it gates everything)

Sox2's canonical co-crystals are **ternary** — Sox2 plus Oct-1 POU on a single duplex.

**This is the irf failure mode by construction**: two proteins on one duplex is exactly what caused the
chain-splitting step to keep one protein and one DNA strand, leaving no base-paired helix. It went
undetected for months.

1. Search for a **binary Sox2–DNA** structure. Prefer it if one exists.
2. If only ternary is available, proceed but restrict the native contact set to **Sox2-chain contacts
   only**, and treat the Oct-1 chain as absent for scoring.
3. **Verify both DNA chains survive into `structures/source_chains/<pdb>_chains/`** before scoring
   anything. RCSB metadata will not predict this — 1IF1 deposits both strands and local prep still
   dropped one.

### 3.1 — Shared core

Standard Chapter 19½ treatment: sequence intersection across 2LE4, 9QPF, the holo chain and the BioEmu
ensemble; strip expression tags; trim disordered termini inward while intra-bundle spread exceeds 3.0 Å.

**Two Sox2-specific notes:**

- 2LE4 models 81 residues; 9QPF's construct is an "extended" HMG domain. The intersection will be
  smaller than either. Record exactly which residues survive — §6 depends on it.
- The terminal tails are where the paper's dynamics result lives, and trimming removes them. **That is
  correct for this test** (see §6) but must be stated.

### 3.2 — Dock and score

Run all four sources through the unmodified Stage-2 dock and `fnat_gate` at `FNAT_FLOOR` = 0.5. Report
median fnat and pass rate per source, as in the three-way table.

**Two traps already documented in the Chapter 19½ README:**

- Scoring the deposited holo *with* its own DNA is near-circular — fnat's native contact set comes from
  that same structure. **Strip the DNA and re-dock.**
- The gate writes `NA` when the reference has no native contacts (`fnat = rec/len(nat) if nat else nan`).
  An all-`NA` file means the **reference chain selection is wrong**, not that the ensemble scored zero.
  With a ternary complex this is the likely first failure.

### 3.3 — Per-residue decomposition ⭐ **the step that makes this a test rather than a coincidence**

2LE4 could fail for boring reasons. It is a 2011 structural-genomics entry selected by target function,
and a lower-quality structure docks worse regardless of conformation. That would be a **quality
artifact**, not a conformational finding.

**The discriminator:** decompose fnat per residue and ask whether 2LE4's loss is **localized to helix
α3** — the specific rearrangement the paper identifies.

| pattern | reading |
|---|---|
| Loss concentrated in α3, rest of the interface scoring like 9QPF | The read-out measured the actual structural difference. **This is the result.** |
| Loss diffuse across the interface | Measuring structure quality, not conformation. Report with the caveat; the test is inconclusive. |

Run the same decomposition on 9QPF as the contrast.

### 3.4 — BioEmu arm

One BioEmu run on the Sox2 construct (~1 GPU-hour), scored through the same dock.

The Finding 3 question: **does BioEmu's best frame beat the deposited apo median** — now with two apo
references at different quality?

| result | reading |
|---|---|
| BioEmu lands nearer 9QPF | Better than the 5-of-5 pilot result would predict. Worth explaining. |
| BioEmu lands nearer 2LE4 | **A concrete, named failure mode** rather than an aggregate statistic: the generator reproduces a conformation the field has established is wrong. Strong, quotable. |
| BioEmu below both | Consistent with the pilots. Extends Finding 3 to a 15th system. |

---

## 4 · Reading the whole result

| 9QPF | 2LE4 | α3 localized | reading |
|---|---|---|---|
| pass | fail | yes | **The strongest validation available.** The read-out independently recovers a published NMR conclusion *and* correctly rejects the misleading deposit, from interface geometry alone. Goes in the paper and the next talk. |
| pass | fail | no | Suggestive, confounded with structure quality. Report with the caveat. |
| pass | pass | — | The read-out is insensitive to the α3 rearrangement. **A real scope limitation** — worth knowing and worth reporting, since it bounds what the mechanism split can resolve. |
| fail | fail | — | Ternary-complex handling is broken. Debug §3.0 before concluding anything. |
| fail | pass | — | Would require explanation before any other conclusion. Check the dock and the contact set first. |

---

## 5 · What it buys beyond validation

**It fills a gap already flagged in Chapter 19½.** The limits section notes the HMG-box split — NHP6A
induced-fit, Sox2 reported as conformational selection — and **lef1, the project's HMG-box pilot, has no
apo structure**. Sox2 gives the family an anchor, and the paper itself notes that *"analysis of Sox2
orthologue structures in free and bound state presents conflicting evidence"* for the structural-change
model, which is exactly the kind of disagreement a geometric read-out can arbitrate.

**It upgrades Finding 3 from a statistic to a named case.** "BioEmu's best frame is below the apo median
in 5 of 5 pilots" is strong; "and on Sox2 it reproduces the deposited structure the field has since
shown to be wrong" is memorable.

**It is a validation of a method, not another datapoint** — which is why it ranks above the 2euv pilot.

---

## 6 · The scope point worth stating explicitly

The paper reports that the HMG **core** has very similar rigidity in free and bound states, while the
**N- and C-terminal tails** become much more rigid upon binding.

The read-out trims disordered termini by design. Chapter 19½ already flags this as its central
limitation, with **ets1** as the case where trimming hid genuine induced fit: Ets-1's autoinhibitory
module unfolds on DNA binding, but it is flanking and was trimmed, so the 103 residues scored were the
preformed ETS domain.

**Sox2 is the mirror image:** core preformed, tails ordered on binding. So "selection-like for the fold"
is the *correct* answer here, and the literature agrees.

Together the two cases define what the measurement sees:

> The read-out answers **"is the DBD fold preformed?"** — reliably. It does **not** answer "does this
> protein undergo induced fit," because the regions where induced fit most often lives are the ones it
> trims. Ets-1 is where that distinction bites; Sox2 is where it doesn't.

**Report them as a pair.** It converts a limitation paragraph into a characterization of the method's
scope, which is a stronger thing to have.

---

## 7 · Risks and checks

| risk | check |
|---|---|
| Ternary complex breaks chain extraction (the irf failure) | Verify two DNA chains in `source_chains/` **before** scoring. Assert at extraction, not at scoring. |
| Native contact set includes Oct-1 contacts | Restrict to Sox2-chain contacts explicitly; inspect the contact list. |
| All-`NA` fnat output | Means wrong reference chain selection, not a zero score. Both scripts should raise `BrokenReference`. |
| 2LE4 fails on quality, not conformation | §3.3 per-residue decomposition. Do not skip it. |
| Shared core trims too much | Record surviving residue count; if the α3 interface residues are trimmed out, the test cannot work. **Check this before running the dock.** |
| Construct differences between 2LE4 and 9QPF | 2LE4 models 81 residues; 9QPF is an extended domain. Confirm the intersection covers the DNA-binding interface. |

---

## 8 · Cost

| step | time |
|---|---|
| Locate holo structure; verify chains | 30 min |
| Shared core across four sources | 1 h |
| Dock + `fnat_gate`, three deposited sources | 1 h |
| Per-residue decomposition | 1 h |
| BioEmu run | ~1 GPU-h, unattended |
| BioEmu scored | 30 min |

**~4 hours hands-on, one GPU-hour.** It is `fnat_on_pairs.py` with a new row and one extra arm.

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

Structures: **2LE4** (Sahu, Markley, Tonelli, Bahrami, Eghbalnia; CESG; deposited 2011, unpublished) and
**9QPF** (Orsetti *et al.*, solution NMR, 200 conformers calculated).
