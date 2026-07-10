# TF-conformation — Analysis

Comprehensive analysis of the six-pilot DeepPBS conformational-augmentation study.
All figures use a cool-pastel palette; all numbers are computed from pipeline outputs on the cluster
(`/project2/rohs_102/shewchuk/TF-conformation`). Single-seed results (5-seed re-run in flight).

Pilots span six DNA-binding-domain families, ordered throughout by fnat gate pass-rate:
**ETS1** (ETS, 100%) · **TBP** (β-saddle, 100%) · **EGR1** (C2H2 zinc-finger, 91%) ·
**engrailed** (homeodomain, 91%) · **FOXA** (forkhead, 75%) · **LEF1** (HMG-box, 19%).

---

## 1 · fnat & gate behaviour

### F1 — fnat distributions per pilot
![F1]({{artifact:4f2502d8-f7ef-4a02-adab-3e57564d2088}})
Post-minimization fraction-of-native-contacts per state, with the 0.5 gate floor. Points colored pass/fail.
ETS1 & TBP sit entirely above the floor; LEF1 has most frames below it (19% pass).

### F2 — gate pass-rate
![F2]({{artifact:12f3eb4b-1015-4f15-ba9a-a21d3689cbd9}})
Fraction of BioEmu frames surviving the fnat gate, with frame counts. Spans 100% (ETS1/TBP) to 19% (LEF1).

### F3 — fnat vs interface-RMSD
![F3]({{artifact:598471bf-7500-41fe-aa9a-e59ac55f6140}})
The two gate criteria co-vary tightly (Spearman ρ = −0.84, n=561): higher iRMSD ⟶ lower fnat.
Pilots cluster along the curve by rigidity.

### F4 — does interface size predict fidelity?
![F4]({{artifact:2b5dc400-dcf2-4199-8ac9-c1aa643408ac}})
No (ρ = −0.06, n.s.). TBP (40 interface residues, 100% pass) and LEF1 (39, 19% pass) are near-identical in
size but opposite in fidelity — faithfulness is about module rigidity, not contact count.

---

## 2 · Interface RMSD (iRMSD)

### I1 — iRMSD distributions
![I1]({{artifact:bf83f8ae-bfa9-4c15-9c70-be7b6da37c32}})
Per-pilot interface RMSD from crystal (post-minimization), same rigidity ordering as fnat.

### I2 — localized vs distributed distortion
![I2]({{artifact:f788b595-0b7e-40cf-a224-958e9783a697}})
Per-segment MAX vs MEAN iRMSD. All points sit above y=x ⟶ interface distortion is **localized**
(one segment always moves more than the rest). LEF1 sprawls to ~6 Å; rigid pilots cluster near origin.

### I3 — Stage 2 → Stage 3 change
![I3]({{artifact:3f8760d0-d252-4326-96a6-f0ea98094ce8}})
Minimization is a **local relaxation**: global backbone barely moves (left), and the per-state fnat change
(right, egr1/foxa — the pilots with Stage-2 interface metrics) is centered near zero with a slight negative skew.

### I4 — interface geometry
![I4]({{artifact:02286a24-3b6e-4ad5-8abd-c8a280099634}})
Interface residue count & segmentation per pilot (crystal-defined constants, not frame-varying).

---

## 3 · Cα-RMSD stage evolution

### R1 — backbone deviation, Stage 2 vs Stage 3
![R1]({{artifact:822c87aa-b1df-4927-8742-2ba39a2ffa50}})
Whole-protein Cα-RMSD from the crystal bound pose, docked (Stage 2) vs minimized (Stage 3).
Rigid recognition modules stay near crystal; the mobile HMG-box (LEF1) drifts far. Stage 2 ≈ Stage 3.

### R2 — per-residue mobility profiles
![R2]({{artifact:b1e514c6-1505-4a76-9bf6-370facad9e31}})
Mean per-residue Cα-RMSD. Universal signature: **rigid core recognition regions, floppy termini**.
LEF1/FOXA have large C-terminal excursions (17–20 Å); TBP's β-saddle is flat (~1 Å) throughout.

### R3 — effect of minimization
![R3]({{artifact:7b834bed-ff67-4536-b735-ea873ebcc645}})
Per-state Cα-RMSD change (Stage 3 − Stage 2). All medians slightly negative (−0.00 to −0.03 Å):
minimization nudges frames marginally **toward** the crystal, never away.

---

## 4 · Structure renders (PyMOL, edu license)

### S1 — best vs worst fnat frame
![S1 ETS1]({{artifact:39991425-335f-48aa-a3c9-f001d81ce85a}})
ETS1 (rigid): best (teal, fnat 0.91) and worst (rose, fnat 0.61) frames both hug the crystal (grey) on the DNA (orange).
![S1 LEF1]({{artifact:e92996b7-6e16-445a-af9e-0820652e6281}})
LEF1 (mobile): worst frame (rose, fnat 0.12) has its N-terminal helix flung off the DNA vs best (teal, 0.71) — why 81% of LEF1 frames fail.

### S2 — crystal vs docked vs minimized
![S2]({{artifact:90100efd-e6c7-4964-80c6-21ae346bc4d6}})
ETS1 state 026: crystal reference (firebrick) vs Stage-2 docked (marine) vs Stage-3 minimized (forest).
Three structures overlap tightly — for a rigid module, docking and minimization preserve the bound pose.

---

## 5 · Pfam family & accuracy

### P1 — pilot → family map
![P1]({{artifact:b0da844a-6f76-4b27-b068-e513ed35ce24}})
The six pilots span six DBD families (Pfam assignments from RCSB PDB).

### P2 — baseline accuracy by family
![P2]({{artifact:fd8be735-32ef-4294-a58c-cb09d85b94b3}})
DeepPBS baseline accuracy over the general-130 benchmark, grouped by family (mean over pilot models).
Family assigned at the **motif level** (each entry inherits the family of the TF whose PWM it evaluates — important because several ETS1 entries come from ETS1–RUNX1 co-crystals that a per-crystal assignment would misfile as Runt). Forkhead and TBP read best; ETS and bZIP are among the harder families.

### P3 — augmentation effect by family  ⭐
![P3]({{artifact:a1708536-68a1-4057-8d59-62de5652f7f1}})
**Key finding (motif-level family assignment):** ETS (n=10, median ΔPearson **+0.042**, 60% of entries improve) and IRF (n=4, +0.039) are the only net-positive families; **every other family is net-negative** (TBP and bHLH worst). ETS is the strongest and best-represented positive family — the ETS1 pilot's gain is **echoed by its ERG/FLI1 paralogs**, not confined to one structure (biggest gains on the low-baseline ERG entries, 0.35→0.51). Both positive families are small-n and single-seed, so this is a lead to test, not a settled result.

> **Note:** an earlier version of this figure assigned family per-PDB from the crystal's first protein chain, which mislabeled ETS1–RUNX1 co-crystal entries as Runt and inflated the ETS median to +0.078. The motif-level assignment above is correct.

---

## Data tables (`analysis/data/`)
- [perstate_metrics.csv]({{artifact:0c000826-9248-4768-8801-bf7fed3a200d}}) — per-state fnat/iRMSD/n_iface_res/seq_ident, Stage 2 + Stage 3, 6 pilots (751 rows). Stage-2 interface metrics exist only for egr1/foxa.
- [ca_rmsd_perstate.csv]({{artifact:fd9f58ab-12c0-429f-9a1c-192f13b09a19}}) — whole-protein Cα-RMSD from crystal (superposed on all protein Cα), Stage 2 + Stage 3 (1134 rows).
- [ca_rmsd_perresidue.csv]({{artifact:f9edd9db-318b-4c76-bb9d-2b89315c79af}}) — per-residue mean Cα-RMSD profile, Stage 3 (612 rows).
- [perentry_accuracy.csv]({{artifact:472230c7-677b-45ec-afca-bf1dd1f3846a}}) — per-entry baseline+augmented metrics over general-130, 6 pilots (780 rows), 15 held-out subset flags.
- [family_annotation.csv]({{artifact:46e5ce42-caf5-44c2-ac9f-68b00952a566}}) — per-entry Pfam family / TF name / motif (RCSB + curated).

## Method notes
- **Cα-RMSD** is whole-protein (superposed on all protein Cα atoms vs the crystal reference), not interface-restricted.
- **fnat / iRMSD** come from the pipeline's fnat-gate output (fraction of native interface contacts; interface RMSD global & per-segment).
- **Family** assigned at the **motif level**: each benchmark entry inherits the family of the TF whose PWM it evaluates (via RCSB Pfam on the motif's PDBs, curated for co-crystals and structures lacking a deposited Pfam mapping). This correctly handles multi-protein complexes (e.g. ETS1–RUNX1 co-crystals) where a per-crystal first-chain assignment would misfile entries.
- Single-seed models; the 5-seed paired re-run will refine the per-family effect estimates.
