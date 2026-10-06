#!/usr/bin/env python
"""Quantify the near-duplicate claim: how close are a state's features to the
crystal's, in the DNA-shape channels vs the protein channels?

If the frozen arm restrains DNA at crystal geometry, the DNA-shape channels of
every state npz should be near-identical to the crystal's, making each state a
near-copy of the crystal in that subspace. Measured, not assumed.
"""
import os, re, glob, json
import numpy as np

R = "/project2/rohs_102/shewchuk/TF-conformation"
PIL = ["csl","egr1","engrailed","err","ets1","foxa","hsf","irf","lef1","nfat","runx","tbp"]
leak = json.load(open(f"{R}/analysis/analyses/conf_vs_crystal_pwm/crystal_leakage.json"))

def meta(tf):
    t = open(f"{R}/config/pilots/{tf}.sh").read()
    g = lambda k: re.search(k+r'="([^"]+)"', t).group(1)
    return g("PDB_ID"), g("PWM_LABEL"), g("BINDING_CHAIN")

z0 = np.load(f"{R}/output/stage5_aug/combined_assembly_ets1/1k79_state_001_MA0098.2.jaspar.npz", allow_pickle=True)
print("FEATURE KEYS:")
for k in z0.files:
    a = z0[k]
    print(f"   {k:22s} {str(getattr(a,'shape',None)):18s} {getattr(a,'dtype',None)}")

def rel_diff(a, b):
    """Median relative deviation over matched-shape arrays."""
    if a.shape != b.shape: return np.nan
    a, b = a.astype(float).ravel(), b.astype(float).ravel()
    sc = np.abs(b).mean()
    if not np.isfinite(sc) or sc == 0: return np.nan
    return float(np.abs(a-b).mean() / sc)

rows = []
for tf in PIL:
    pid, pwm, bc = meta(tf)
    d = f"{R}/output/stage5_aug/combined_assembly_{tf}"
    cry = leak[tf]["crystal_entry"]
    cp = os.path.join(d, cry)
    if not os.path.exists(cp): continue
    zc = np.load(cp, allow_pickle=True)
    states = sorted(glob.glob(os.path.join(d, f"{pid}_state_*_{pwm}.npz")))[:40]
    per = {}
    for sp in states:
        zs = np.load(sp, allow_pickle=True)
        for k in zc.files:
            if k not in zs.files: continue
            try:
                v = rel_diff(zs[k], zc[k])
            except Exception:
                continue
            if v is not None and np.isfinite(v):
                per.setdefault(k, []).append(v)
    for k, vs in per.items():
        rows.append(dict(tf=tf, key=k, n=len(vs), med_rel_dev=float(np.median(vs))))

import csv
with open(f"{R}/analysis/analyses/conf_vs_crystal_pwm/feature_similarity.csv","w",newline="") as f:
    w = csv.DictWriter(f, fieldnames=["tf","key","n","med_rel_dev"]); w.writeheader(); w.writerows(rows)

print("\nMEDIAN RELATIVE DEVIATION of state features from the crystal's (per key, across pilots):")
keys = sorted(set(r["key"] for r in rows))
for k in keys:
    v = [r["med_rel_dev"] for r in rows if r["key"] == k]
    print(f"   {k:22s} n_pilots={len(v):2d}  median={np.median(v):.4f}  range={np.min(v):.4f}-{np.max(v):.4f}")
print("\nper-pilot, shape-only keys:")
for tf in PIL:
    sel = {r["key"]: r["med_rel_dev"] for r in rows if r["tf"] == tf}
    if sel: print(f"   {tf:10s}", {k: round(v,4) for k,v in sorted(sel.items())})
