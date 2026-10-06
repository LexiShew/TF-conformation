#!/usr/bin/env python
"""Is a conformation a NEAR-DUPLICATE of its crystal, in the model's own input space?

The right yardstick is not absolute deviation but deviation relative to the
spread of the training set. For each pilot:
   d_own   = distance from each of its states to its own crystal
   d_other = distance from unrelated training entries to that same crystal
A near-duplicate regime means d_own << d_other. If d_own is comparable to
d_other, the states are not near-copies in the space the model actually sees.

Distances are computed per DNA position on the DNA node-feature block (which
carries the 3DNA shape parameters), aggregated as RMS over positions after
matching length, and reported as the ratio d_own / d_other.
"""
import os, re, glob, json
import numpy as np

R = "/project2/rohs_102/shewchuk/TF-conformation"
PIL = ["csl","egr1","engrailed","err","ets1","foxa","hsf","irf","lef1","nfat","runx","tbp"]
leak = json.load(open(f"{R}/analysis/analyses/conf_vs_crystal_pwm/crystal_leakage.json"))
ASSEMBLY = f"{R}/data/assembly2024"

def meta(tf):
    t = open(f"{R}/config/pilots/{tf}.sh").read()
    g = lambda k: re.search(k+r'="([^"]+)"', t).group(1)
    return g("PDB_ID"), g("PWM_LABEL"), g("BINDING_CHAIN")

def dna_block(path):
    z = np.load(path, allow_pickle=True)
    if "X_dna" not in z.files: return None
    a = np.asarray(z["X_dna"], dtype=float)
    return a if a.ndim == 2 else None

def rms_dist(a, b):
    """RMS difference over the overlapping leading positions, per channel-mean."""
    if a is None or b is None: return np.nan
    if a.shape[1] != b.shape[1]: return np.nan
    n = min(a.shape[0], b.shape[0])
    if n < 4: return np.nan
    d = a[:n] - b[:n]
    return float(np.sqrt(np.mean(d**2)))

# unrelated reference entries: a fixed sample of the upstream training assembly
others = sorted(glob.glob(os.path.join(ASSEMBLY, "*.npz")))
rng = np.random.default_rng(0)
others = list(rng.choice(others, size=min(60, len(others)), replace=False))
other_blocks = []
for p in others:
    b = dna_block(p)
    if b is not None: other_blocks.append((os.path.basename(p), b))

print(f"reference pool: {len(other_blocks)} unrelated training entries")
print(f"{'tf':10s} {'n_st':>5s} {'d_own':>8s} {'d_other':>8s} {'ratio':>7s} {'r_chan':>7s}")
rows = []
for tf in PIL:
    pid, pwm, bc = meta(tf)
    d = f"{R}/output/stage5_aug/combined_assembly_{tf}"
    cp = os.path.join(d, leak[tf]["crystal_entry"])
    C = dna_block(cp)
    if C is None: continue
    sts = sorted(glob.glob(os.path.join(d, f"{pid}_state_*_{pwm}.npz")))[:60]
    dn = [rms_dist(dna_block(s), C) for s in sts]
    dn = [v for v in dn if np.isfinite(v)]
    do = [rms_dist(B, C) for _, B in other_blocks]
    do = [v for v in do if np.isfinite(v)]
    if not dn or not do: continue
    # channel-wise correlation of a median state against the crystal
    mid = dna_block(sts[len(sts)//2])
    n = min(mid.shape[0], C.shape[0])
    rc = float(np.corrcoef(mid[:n].ravel(), C[:n].ravel())[0,1])
    ratio = np.median(dn)/np.median(do)
    rows.append(dict(tf=tf, n_states=len(dn), d_own=np.median(dn), d_other=np.median(do),
                     ratio=ratio, r_channels=rc))
    print(f"{tf:10s} {len(dn):5d} {np.median(dn):8.4f} {np.median(do):8.4f} {ratio:7.3f} {rc:7.3f}")

import csv
with open(f"{R}/analysis/analyses/conf_vs_crystal_pwm/nearduplicate_ratio.csv","w",newline="") as f:
    w = csv.DictWriter(f, fieldnames=["tf","n_states","d_own","d_other","ratio","r_channels"])
    w.writeheader(); w.writerows(rows)
rr = [r["ratio"] for r in rows]
print(f"\nmedian ratio d_own/d_other = {np.median(rr):.3f}  (<<1 => near-duplicate regime)")
print(f"mean channel-correlation state-vs-crystal = {np.mean([r['r_channels'] for r in rows]):.3f}")
