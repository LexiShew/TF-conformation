# Pilot-candidate selection & multi-DBD policy
**Scope:** `/project2/rohs_102/shewchuk/TF-conformation` · test set = `data/folds/id.txt` (130 entries, 84 PDBs)

## 0 · Two findings that precede the candidate list

### 0.1 The TBP pilot carries the wrong PWM label (verified)
`config/pilots/tbp.sh` sets `PWM_LABEL="MA0343.1.jaspar"`. Audited against the master dataset
(`stage4_preprocess/jaspar_h11mo_cluster_wise_dna_containing_dataset.npy`), **12/13 pilots match
their structure's assigned PWM; TBP is the sole mismatch**:

| | |
|---|---|
| `1tgh_A` dataset label | `TBP_HUMAN.H11MO.0.A` |
| `tbp.sh` config label | `MA0343.1.jaspar` |
| JASPAR MA0343.1 identity | **NDT80**, *S. cerevisiae*, "NDT80 domain factors" (UniProt P38830) |
| PDB 1TGH | human TATA-binding protein |

MA0343.1 is not a TBP matrix. The five `MA0343.1` test entries (`2etw/2euz/2evf/2evi/2evj`) are all
**NDT80** structures (Pfam PF05224, NDT80/PhoG-like) — a different clan from TBP's PF00352/TBP-like.
The real TBP-family test anchor is a single entry, `7oh9_K_MA0386.1` (SPT15, yeast TBP) — and 7OH9 is
a nucleosome+TFIIA complex, not a clean TBP–DNA monomer.

Why this matters beyond bookkeeping: `stage5_build_aug` uses `PWM_LABEL` as a **substring filter**
(`build_augmented_fold.py`, `build_combined_assembly.py`) to decide which entries the ~95 augmented
TBP conformers join. With `MA0343.1`, every TBP conformer was injected under NDT80's label, so the
"TBP" arm trained TBP backbones against an NDT80 motif and was scored on NDT80 test entries. TBP is
the project's headline negative result ("uniformly worse, all significant" — `pilot_assessment_memo`)
and one of the two protagonists of the book outline. **That result should be re-checked before it is
written up**; the label error is a live alternative explanation for it that is independent of the
DNA-deformation mechanism.

I have not changed anything — flagging only. Fix is one line in `tbp.sh` plus a stage-5-onward rerun.

### 0.2 The IRF failure is not a missing strand in the PDB
Worth correcting because it changes how you screen candidates. PDB **1IF1 deposits both strands**
(entity `polydeoxyribonucleotide`, auth chains C and D, 26 nt, self-complementary). The loss happens
locally: `structures/source_chains/1if1_chains/1if1.cif` contains only **2 chains** — one DNA, one
protein — and `1if1_dna.pdb` holds a single strand. The ASU has **two IRF copies on one duplex**, and
the chain-splitting step kept one protein + one strand, so stage-4 `processDNA` gets no base-paired
helix (the "helix count problem 0" class already documented for 5D5U/HSF).

So there are two distinct hazards to screen for, and they need different fixes:
- **Type A — deposition has one strand** (true palindromic ASU; duplex is a symmetry mate).
  Fix: regenerate the biological assembly from `_pdbx_struct_oper_list` (the documented 5D5U recipe).
- **Type B — deposition has both strands, local prep dropped one** (IRF/1IF1). Fix: re-extract the
  reference CIF keeping both DNA chains; no symmetry operator needed.

**Screen used for every candidate below:** reject if the entry has only one DNA chain, and reject if
the protein:duplex stoichiometry is >1 (two protein copies per duplex is what broke IRF).

---

## 1 · What has already run (13 pilots, 8 Pfam clans)

| pilot | PDB | DBD family (Pfam) | clan | on-target test entries |
|---|---|---|---|---|
| csl | 3brg | PF09270 + PF09271 (+PF20144 TIG) | Trefoil / P53-like / E-set | 5 |
| dux4 | 5z6z | PF00046 Homeodomain ×2 | HTH | 9 |
| egr1 | 1aay | PF00096 C2H2 ×3 | C2H2-zf | 19 |
| engrailed | 3hdd | PF00046 Homeodomain | HTH | 9 |
| err | 1lo1 | PF00105 C4 zinc finger | Tbcl_zf | 10 |
| ets1 | 1k79 | PF00178 Ets | HTH | 10 |
| foxa | 1vtn | PF00250 Forkhead | HTH | 5 |
| hsf | 5d5u | PF00447 HSF-type | HTH | 3 |
| irf | 1if1 | PF00605 IRF | HTH | 4 | 
| lef1 | 2lef | PF00505 HMG-box | HMG-box | **0** |
| nfat | 1a66 | PF00554 Rel homology | P53-like | 15 |
| runx | 1hjc | PF00853 Runt | P53-like | 5 |
| tbp | 1tgh | PF00352 TBP | TBP-like | 1 (mislabelled, §0.1) |

Note **lef1 has zero same-family test entries** — criterion 4 is not actually satisfied for it, which
is consistent with it only ever being evaluated on the general-130.

### Clan-level saturation
Of the 13 clans present in the test set, 7 are already covered. The uncovered ones are
**bZIP (8 entries), MADS/SRF (12), bHLH (4), Zn2Cys6 (2)** — and all four are *obligate dimers*:
every non-test structure in those clans has ≥2 protein chains on the duplex, i.e. the IRF failure
mode by construction. They are closed to BioEmu as whole clans, with exactly one exception (SKN-1,
below). The remaining headroom is therefore mostly **family-level novelty inside already-covered
clans**, which still satisfies "new Pfam family" and keeps criterion 4 satisfiable.

---

## 2 · Recommended candidates

All verified against RCSB: single protein chain, two DNA chains (real duplex), not in `id.txt`.

| # | TF | PDB | chain | new family | clan | DBDs | res (Å) | PWM | test anchors |
|---|---|---|---|---|---|---|---|---|---|
| 1 | **NDT80** | 2euv | A | PF05224 NDT80/PhoG | P53-like | 1 | 1.94 | MA0343.1 | 5 (2etw,2euz,2evf,2evi,2evj) |
| 2 | **p53 core** | 7eds | A | PF00870 p53 DBD | P53-like | 1 | 1.77 | MA0106.3 / P53_HUMAN | 5 (2ady,3q05,6fj5×3) |
| 3 | **Oct-1 POU** | 1oct | C | PF00157 POU-spec + PF00046 | HTH | **2** | 3.00 | MA0785.1 | 1 POU + 9 homeodomain |
| 4 | **Brn-1/POU3F1** | 2xsd | C | PF00157 POU-spec + PF00046 | HTH | **2** | 2.05 | PO3F1_MOUSE | 1 POU + 9 homeodomain |
| 5 | **PBX1 (TALE)** | 1lfu | P | PF05920 Homeobox-KN + PF00046 | HTH | **2** | NMR | PBX1_MOUSE | 3 (4xrm×2, 5bng) |
| 6 | **SKN-1** | 1skn | P | PF03131 bZIP Maf | bZIP | 1 | 2.50 | MA0547.1 | 2 (2wty A/B) |

**Why these, briefly.**

1. **NDT80 (2euv)** — strongest pick. New family, 1.94 Å, clean monomer+duplex, and it is the
   *only* candidate that is in **neither** `train*.txt` nor `valid*.txt` nor `id.txt` (verified),
   so augmenting with it cannot double-count an entry the baseline already trains on. Its five test
   anchors are a tight homogeneous family. It also directly repairs the MA0343.1 confusion: run
   under MA0343.1, NDT80 conformers finally reach the test entries that label actually belongs to.
   Sibling structures 2euw/2eux/2evg/2evh/1mnn are in the training folds — prefer **2euv**.
2. **p53 (7eds)** — new family in a covered clan, best resolution on the list (1.77 Å), 30-bp duplex,
   5 anchors. Deliberately **not** 7b4d/7b4f/7b4e: those deposit a single DNA strand (two are
   self-complementary) = exactly the IRF/5D5U trap. 7EDS is the one clean monomeric p53–DNA entry.
3–5. **The multi-DBD trio** — see §3; these are the cases your question is about, and they are the
   scientifically interesting ones. 2xsd is the better-resolved twin of 1oct (2.05 vs 3.00 Å); run
   **2xsd first** and keep 1oct as the replicate. 1lfu is NMR, so it has no crystal-packing artifact
   but also no resolution — a useful contrast given the project's crystal-is-one-frame thesis.
6. **SKN-1 (1skn)** — the single exception to the dimer wall. SKN-1 binds DNA as a *monomer* using a
   Maf-like basic region with no leucine zipper, so it opens the bZIP clan (8 test entries) to BioEmu
   when no other bZIP can be modelled. 2.50 Å, 2 direct Maf anchors. Lower priority only because
   the basic region is largely disordered when free — expect low fnat yield, which is itself a
   datapoint for the rigidity/fnat relationship in the assessment memo.

**Suggested order:** 2euv → 7eds → 2xsd → 1lfu → 1skn → 1oct. The first two are low-risk family
additions; 2xsd/1lfu test the multi-DBD policy; 1skn is the high-variance probe.

---

## 3 · How to handle TFs with multiple DNA-binding domains

This is the part of the request that needs a policy rather than a list. First, the scale of it —
measured, not assumed:

| set | chains | >1 Pfam clan | >1 Pfam family |
|---|---|---|---|
| test set (130) | 124 | 25 (20%) | 31 (25%) |
| monomer pool | 153 | 2 (1%) | 12 (8%) |
| all DeepPBS | 3245 | 993 (31%) | 1141 (35%)|

A quarter of the test set is multi-domain, so this is not an edge case; and the existing tooling does
not handle it. `structures/pfam/filter_pfam.py` groups by `(pdb_id, chain, clan_acc)` and keeps the
dominant family **within each clan** — so a chain with two DBDs from *different* clans keeps both
rows and is silently counted in two clan bars, while a chain with two *different* families in the
*same* clan (POU-specific + homeodomain, both HTH) has one of its two real DBDs deleted. Neither
behaviour is what you want for pilot selection.

### Recommended policy — classify by the *recognition unit*, not by the chain

**Rule 1 — Distinguish three cases before anything else.**
- (a) **Tandem repeats of one family** (EGR1 PF00096 ×3, DUX4 PF00046 ×2, TBP PF00352 ×2).
  Already handled correctly: one family, N copies. Record `n_domains` but classify as that family.
  This is the common case and needs no new policy.
- (b) **Accessory/non-DBD domains alongside the DBD** (NF-κB Rel PF00554 + Ig-like PF16179;
  CSL PF20144 TIG; p53 PF07710 tetramerisation; Ets PF19525 flanking region).
  These inflate the clan count and should be **excluded from the classification key** via an
  explicit non-DBD denylist. They are dimerisation/scaffold modules, not readout modules.
  This alone removes most of the test set's 20% multi-clan rate — E-set/P53-like pairs dominate it.
- (c) **Genuinely multiple distinct DBDs that both touch DNA** (POU: PF00157 + PF00046;
  TALE/PBX: PF05920 + PF00046). This is the real case.

**Rule 2 — For case (c), make the classification key the sorted *set* of contacting DBD families,
not a single label.** Classify Oct-1 as `{PF00046, PF00157}`, not "POU" or "homeodomain". Then:
- criterion 3 (*covers a new family/clan*) is satisfied if the **set** is new — `{PF00046, PF00157}`
  is new even though PF00046 alone is covered by dux4/engrailed;
- criterion 4 (*test set has a within-family entry*) should be evaluated **per component family**,
  and you should record *which* component the anchor matches. For 1oct you get 1 anchor via PF00157
  (`1e3o_C`) and 9 via PF00046. That asymmetry is informative, not a problem.

**Rule 3 — Decide "does it contact DNA" from the structure, not from Pfam.** Pfam tells you a domain
exists; it does not tell you whether it is at the interface in *this* complex. You already compute
exactly what is needed — `fnat_gate/` / `interface_residues` — so define a domain as DNA-contacting
if ≥1 residue within its Pfam span is in the crystal interface set (≤4.5 Å to DNA). Run this once over
all candidate chains and cache it; it converts the denylist in Rule 1b from a hand-curated judgement
into a measured property, and it will catch cases the denylist misses.

**Rule 4 — Report ensemble drift per domain, not per chain.** This is where multi-DBD pilots pay off
scientifically. The memo's central mechanism is that fnat yield and augmentation benefit track Cα
drift from the bound pose. For a two-DBD TF, the dominant free-state motion is usually the **relative
rigid-body placement of the two domains about their linker**, which a single whole-chain RMSD
averages into noise. Compute (i) per-domain drift after superposing on that domain alone, and
(ii) inter-domain drift (superpose on domain 1, measure domain 2). POU is the textbook case: two
small rigid subdomains on a flexible 15–25 residue linker. The prediction is per-domain drift low,
inter-domain drift high — i.e. BioEmu gets each domain right and the *assembly* wrong. If
augmentation hurts there while per-domain drift stays low, that is a cleaner demonstration of the
project's thesis than any single-DBD pilot can give, because it separates "wrong fold" from
"wrong arrangement."

**Rule 5 — Keep the fnat gate per-domain too.** A whole-interface fnat can pass at 0.5 while one of
two DBDs has completely left its site. Gate on the minimum across domains, not the mean.

### Concretely, for the tooling
- add a `non_dbd_pfam.txt` denylist (seed: PF07710, PF19525, PF16179, PF20144, plus PAS/Myc-LZ as encountered);
- replace `filter_pfam.py`'s `(pdb,chain,clan)` grouping with `(pdb,chain)` → sorted DBD-family set,
  collapsing tandem repeats to `family×N`;
- add the interface-overlap test from Rule 3 as the gate on set membership;
- emit `n_dbd_families` and `dbd_set` columns into the pilot config so stages 3–7 can key off them.

---

## 4 · Caveats
- Pfam annotations here come from the cached SIFTS pulls in `structures/pfam/` (June 2026), not a
  fresh InterPro query; family assignments for the 6 candidates were spot-checked against
  RCSB/JASPAR but the full tables were not re-derived.
- Structural screening used RCSB deposition metadata (chain counts, entity types, self-complementarity).
  It catches Type-A single-strand entries but **cannot** predict Type-B local-prep losses like IRF —
  verify both DNA chains survive into `structures/source_chains/<pdb>_chains/` before launching.
- `1lfu` is an NMR ensemble; stage-1/stage-2 assume a single reference conformer — check which model
  the prep picks.
- Resolution 3.00 Å for 1oct is marginal for interface work; 2xsd is preferred for that reason.
