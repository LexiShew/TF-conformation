#!/usr/bin/env python3
"""Build a base-pair-step -> stiffness lookup table from the hexABC ensemble.

Reads sum_stiffness.json (summed diagonal elastic constant per dinucleotide step)
across ALL hexABC sequences, groups by step label (AA/TT, AC/GT, TA/TA, ...) and reports mean/std/n per step.
hexABC keeps both reading directions, so this yields 16 distinct dinucleotide-step
labels (NOT collapsed to the 10 unique WC pairs). This step-stiffness table is the
sequence-dependent prior for mapping onto per-position k_dna.

Also emits the per-step breakdown for the 6 helical DOF (twist/roll/tilt/shift/
slide/rise) so a richer per-DOF restraint is possible later.

Output: stiffness_table.json  { "sum": {step: {mean,std,n}}, "per_dof": {...} }
"""
import os, json, glob
from collections import defaultdict
import statistics as st

HEXABC="/project2/rohs_102/share/HexABC_data"
DOF=["sum","twist","roll","tilt","shift","slide","rise"]

def norm_step(label):
    # labels are "XY/WZ" (step on strand1 / its WC complement); use as-is.
    # hexABC keeps both reading directions -> 16 distinct step labels, not 10.
    return label.strip()

acc={dof:defaultdict(list) for dof in DOF}
seqdirs=sorted(glob.glob(f"{HEXABC}/seq*"))
n_seq=0
for sd in seqdirs:
    adir=f"{sd}/analyses/average"
    if not os.path.isdir(adir): continue
    n_seq+=1
    for dof in DOF:
        fn=f"{adir}/{dof}_stiffness.json"
        if not os.path.exists(fn): continue
        try:
            js=json.load(open(fn))
            key=f"{dof}_stiffness"
            for entry in js.get(key,[]):
                lab=norm_step(entry.get("label",""))
                val=entry.get("value")
                if lab and val is not None:
                    acc[dof][lab].append(float(val))
        except (json.JSONDecodeError, OSError, ValueError):
            continue

out={}
for dof in DOF:
    out[dof]={}
    for step,vals in acc[dof].items():
        if len(vals)==0: continue
        out[dof][step]={"mean":round(st.mean(vals),5),
                        "std":round(st.pstdev(vals),5) if len(vals)>1 else 0.0,
                        "n":len(vals)}

meta={"n_sequences":n_seq,"source":HEXABC,"note":"per-bp-step elastic constants averaged across hexABC MD ensemble; 'sum'=summed diagonal stiffness"}
result={"meta":meta,"stiffness":out}
json.dump(result, open("stiffness_table.json","w"), indent=2)
# console summary: sum_stiffness per canonical step, sorted
print(f"n_sequences={n_seq}")
print("=== sum_stiffness per step (mean, sorted flexible->stiff) ===")
sm=sorted(out["sum"].items(), key=lambda kv: kv[1]["mean"])
for step,d in sm:
    print(f"  {step:8s} mean={d['mean']:7.3f}  std={d['std']:6.3f}  n={d['n']}")
