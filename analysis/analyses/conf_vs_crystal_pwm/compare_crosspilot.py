#!/usr/bin/env python
"""Cross-pilot generalization: does augmentation help on a complex the model
has NEVER seen, at the structure level?

Design. For target pilot X and donor pilot Y (Y != X), score X's conformations
and X's crystal under Y's baseline and Y's augmented checkpoints. Y's augmented
arm was trained on Y's conformations only, so X's complex is unseen at the
structure level -- not merely absent from the entry list. The paired quantity

    delta(X, Y) = r_augmented_Y(X) - r_baseline_Y(X)

isolates what augmentation buys on a genuinely novel complex, holding the
scored structures fixed. The within-pilot value delta(X, X) from the earlier
report is the leaky upper bound on the same quantity.

Two readouts, matching the earlier report:
  state -> crystal  : self-consistency of the conformer ensemble
  crystal -> experiment : does the model read an unseen crystal better
"""
import os, re, glob, json
import numpy as np
import pandas as pd

R = "/project2/rohs_102/shewchuk/TF-conformation"
TRAIN = f"{R}/output/stage6_train"
IDD = f"{R}/stage7_eval/confcmp_idfiles"
OUT = f"{R}/analysis/analyses/conf_vs_crystal_pwm"
MIN_OVERLAP = 4
PIL = ["csl","egr1","engrailed","err","ets1","foxa","hsf","irf","lef1","nfat","runx","tbp"]
leak = json.load(open(f"{OUT}/crystal_leakage.json"))

def norm(p):
    s = p.sum(-1, keepdims=True)
    return np.divide(p, s, out=np.full_like(p, 0.25), where=s > 0)

def to_frame(z, key="P"):
    n = z["Y"].shape[0]
    f = np.full((n, 4), np.nan)
    cols = np.flatnonzero(z["Y_mask"].astype(bool))
    vals = z[key][z[key+"_mask"].astype(bool)]
    k = min(len(cols), len(vals))
    f[cols[:k]] = norm(vals[:k])
    return f

def exp_frame(z):
    n = z["Y"].shape[0]
    f = np.full((n, 4), np.nan)
    m = z["Y_mask"].astype(bool)
    f[m] = norm(z["Y"][m])
    return f

def pear(a, b):
    m = ~(np.isnan(a).any(-1) | np.isnan(b).any(-1))
    if m.sum() < MIN_OVERLAP: return np.nan, int(m.sum())
    x, y = a[m].ravel(), b[m].ravel()
    if x.std() == 0 or y.std() == 0: return np.nan, int(m.sum())
    return float(np.corrcoef(x, y)[0,1]), int(m.sum())

def ic(p):
    q = p[~np.isnan(p).any(-1)]
    if q.size == 0: return np.nan
    return float(np.sum(2.0 + np.sum(np.where(q>0, q*np.log2(np.where(q>0,q,1)), 0.0), -1)))

rows = []
for X in PIL:
    cry = leak[X]["crystal_entry"]
    idf = f"{IDD}/id_{X}_confcmp.txt"
    if not os.path.exists(idf): continue
    entries = [l.strip() for l in open(idf) if l.strip()]
    states = [e for e in entries if "_state_" in e]
    for Y in PIL:
        same = (X == Y)
        for kind in ["baseline","augmented"]:
            sub = "predictions_confcmp" if same else "predictions_crosspilot"
            for outer in sorted(glob.glob(f"{TRAIN}/{kind}_{Y}_fold0_s*")):
                b = os.path.basename(outer)
                if "dnarelax" in b or ".incomplete" in b: continue
                m = re.search(r"_s(\d+)$", b)
                if not m: continue
                seed = int(m.group(1))
                pdir = f"{outer}/{b}/{sub}"
                cp = f"{pdir}/{cry}_predict.npz"
                if not os.path.exists(cp): continue
                zc = np.load(cp)
                Pc, Ye = to_frame(zc), exp_frame(zc)
                r_cry_exp, n_ce = pear(Pc, Ye)
                srs, sxs = [], []
                for e in states:
                    sp = f"{pdir}/{e}_predict.npz"
                    if not os.path.exists(sp): continue
                    Ps = to_frame(np.load(sp))
                    if Ps.shape != Pc.shape: continue
                    r1, n1 = pear(Ps, Pc)
                    r2, _ = pear(Ps, Ye)
                    if np.isfinite(r1): srs.append(r1)
                    if np.isfinite(r2): sxs.append(r2)
                if not srs: continue
                rows.append(dict(target=X, donor=Y, same_pilot=same, kind=kind, seed=seed,
                                 n_states=len(srs),
                                 r_state_crystal=float(np.mean(srs)),
                                 r_state_exp=float(np.mean(sxs)) if sxs else np.nan,
                                 r_crystal_exp=r_cry_exp, crystal_ic=ic(Pc),
                                 target_leakage=("exact" if leak[X]["exact_in_baseline"]
                                     else "structure_only" if leak[X]["same_structure_other_label_baseline"]
                                     else "held_out")))

df = pd.DataFrame(rows)
df.to_csv(f"{OUT}/crosspilot_per_pair.csv", index=False)

# paired delta per (target, donor, seed)
w = df.pivot_table(index=["target","donor","same_pilot","seed","target_leakage"],
                   columns="kind", values=["r_state_crystal","r_crystal_exp","r_state_exp"]).reset_index()
w.columns = ["__".join([c for c in col if c]) if isinstance(col, tuple) else col for col in w.columns]
for met in ["r_state_crystal","r_crystal_exp","r_state_exp"]:
    a, b = f"{met}__augmented", f"{met}__baseline"
    if a in w.columns and b in w.columns:
        w[f"d_{met}"] = w[a] - w[b]
w.to_csv(f"{OUT}/crosspilot_deltas.csv", index=False)

pd.set_option("display.width", 220); pd.set_option("display.max_columns", 30)
print("pairs scored:", len(df), "| unique target-donor:", df.groupby(['target','donor']).ngroups)
print("\n=== augmentation delta: WITHIN pilot (leaky) vs CROSS pilot (clean) ===")
g = w.groupby("same_pilot")[["d_r_state_crystal","d_r_crystal_exp","d_r_state_exp"]].agg(["mean","median","size"])
print(g.round(3).to_string())
print("\n=== cross-pilot delta, by target leakage class ===")
print(w[~w.same_pilot].groupby("target_leakage")[["d_r_state_crystal","d_r_crystal_exp"]].agg(["mean","size"]).round(3).to_string())
print("\n=== cross-pilot delta per target (averaged over donors and seeds) ===")
pt = w[~w.same_pilot].groupby("target")[["d_r_state_crystal","d_r_crystal_exp"]].mean()
wi = w[w.same_pilot].groupby("target")[["d_r_state_crystal","d_r_crystal_exp"]].mean()
cmpt = pt.join(wi, lsuffix="_cross", rsuffix="_within")
print(cmpt.round(3).sort_values("d_r_state_crystal_cross").to_string())
