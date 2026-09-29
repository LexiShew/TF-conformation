#!/usr/bin/env python
"""Baseline-DeepPBS PWM predictions: BioEmu/minimized bound conformations vs
the crystal complex prediction, compared in shared PWM coordinates.

DeepPBS emits, per entry, a per-DNA-position prediction P (length = duplex
length) plus P_mask marking the positions that were aligned to the reference
PWM, and Y/Y_mask giving the experimental PWM and which of ITS columns were
covered. The masked P positions correspond IN ORDER to the masked Y columns,
so every entry can be lifted into the PWM's own coordinate frame regardless of
how many duplex positions its alignment covered. Conformations whose alignment
window differs in length from the crystal's are therefore still comparable, on
the intersection of covered PWM columns.
"""
import os, sys, glob, json, re
import numpy as np
import pandas as pd

R = "/project2/rohs_102/shewchuk/TF-conformation"
TRAIN = os.path.join(R, "output/stage6_train")
IDD = os.path.join(R, "stage7_eval/confcmp_idfiles")
OUT = os.path.join(R, "analysis/analyses/conf_vs_crystal_pwm")
os.makedirs(OUT, exist_ok=True)
MIN_OVERLAP = 4   # PWM columns required for a comparison to be reported

PILOTS = ["csl","egr1","engrailed","err","ets1","foxa","hsf","irf","lef1","nfat","runx","tbp"]

def pilot_meta(tf):
    txt = open(os.path.join(R, "config/pilots", tf + ".sh")).read()
    g = lambda k: re.search(k + r'="([^"]+)"', txt).group(1)
    return g("PDB_ID"), g("PWM_LABEL"), g("BINDING_CHAIN")

def crystal_entry(tf):
    pid, pwm, bc = pilot_meta(tf)
    return "5d5u_duplex_MA0486.2.jaspar.npz" if tf == "hsf" else f"{pid}_{bc}_{pwm}.npz"

def norm(p):
    s = p.sum(axis=-1, keepdims=True)
    return np.divide(p, s, out=np.full_like(p, 0.25), where=s > 0)

def to_pwm_frame(z, key="P"):
    """Lift a prediction into PWM coordinates: (n_pwm_cols, 4) with NaN where
    this entry's alignment did not cover the column."""
    n = z["Y"].shape[0]
    frame = np.full((n, 4), np.nan)
    cols = np.flatnonzero(z["Y_mask"].astype(bool))
    vals = z[key][z[key + "_mask"].astype(bool)]
    k = min(len(cols), len(vals))
    frame[cols[:k]] = norm(vals[:k])
    return frame

def exp_pwm(z):
    n = z["Y"].shape[0]
    frame = np.full((n, 4), np.nan)
    m = z["Y_mask"].astype(bool)
    frame[m] = norm(z["Y"][m])
    return frame

def shared(a, b):
    m = ~(np.isnan(a).any(axis=-1) | np.isnan(b).any(axis=-1))
    return a[m], b[m], int(m.sum())

def pearson(a, b):
    a, b = a.ravel(), b.ravel()
    if a.size < 2 or a.std() == 0 or b.std() == 0: return np.nan
    return float(np.corrcoef(a, b)[0, 1])

def mae(a, b): return float(np.abs(a - b).mean())

def jsd(a, b):
    m = 0.5 * (a + b)
    def kl(x, y):
        t = np.where(x > 0, x * np.log2(np.divide(x, y, out=np.ones_like(x), where=y > 0)), 0.0)
        return np.sum(t, axis=-1)
    return float(np.mean(0.5 * kl(a, m) + 0.5 * kl(b, m)))

def ic(p):
    p = p[~np.isnan(p).any(axis=-1)]
    if p.size == 0: return np.nan
    return float(np.sum(2.0 + np.sum(np.where(p > 0, p * np.log2(np.where(p > 0, p, 1)), 0.0), axis=-1)))

def cons(p):
    return "".join("." if np.isnan(row).any() else "ACGT"[int(np.argmax(row))] for row in p)

def score(a, b, tag):
    A, B, n = shared(a, b)
    if n < MIN_OVERLAP:
        return {f"r_vs_{tag}": np.nan, f"mae_vs_{tag}": np.nan,
                f"jsd_vs_{tag}": np.nan, f"n_overlap_{tag}": n}
    return {f"r_vs_{tag}": pearson(A, B), f"mae_vs_{tag}": mae(A, B),
            f"jsd_vs_{tag}": jsd(A, B), f"n_overlap_{tag}": n}

rows, skipped = [], []
bundles = {}

for tf in PILOTS:
    pid, pwmlab, bc = pilot_meta(tf)
    cry_name = crystal_entry(tf)
    entries = [l.strip() for l in open(os.path.join(IDD, f"id_{tf}_confcmp.txt")) if l.strip()]
    states = [e for e in entries if "_state_" in e]

    for outer in sorted(glob.glob(os.path.join(TRAIN, f"baseline_{tf}_fold0_s*"))):
        if "dnarelax" in outer or ".incomplete" in outer: continue
        m = re.search(r"_s(\d+)$", os.path.basename(outer))
        if not m: continue
        seed = int(m.group(1))
        pdir = os.path.join(outer, os.path.basename(outer), "predictions_confcmp")
        cpath = os.path.join(pdir, cry_name + "_predict.npz")
        if not os.path.exists(cpath):
            skipped.append((tf, seed, "ALL", "no crystal prediction")); continue

        zc = np.load(cpath)
        Pc, Yexp = to_pwm_frame(zc, "P"), exp_pwm(zc)
        base = dict(tf=tf, pdb=pid, pwm_label=pwmlab, seed=seed)

        rows.append({**base, "entry": cry_name, "kind": "crystal", "state": np.nan,
                     "n_cols_covered": int((~np.isnan(Pc).any(axis=-1)).sum()),
                     **score(Pc, Pc, "crystal"), **score(Pc, Yexp, "exp"),
                     "ic_bits": ic(Pc), "consensus": cons(Pc)})

        stack, sids = [], []
        for e in states:
            sp = os.path.join(pdir, e + "_predict.npz")
            if not os.path.exists(sp):
                skipped.append((tf, seed, e, "prediction file absent")); continue
            zs = np.load(sp)
            Ps = to_pwm_frame(zs, "P")
            if Ps.shape != Pc.shape:
                skipped.append((tf, seed, e, f"different PWM length {Ps.shape} vs {Pc.shape}")); continue
            sc = score(Ps, Pc, "crystal")
            if np.isnan(sc["r_vs_crystal"]):
                skipped.append((tf, seed, e, f"overlap {sc['n_overlap_crystal']} < {MIN_OVERLAP}")); continue
            stack.append(Ps); sids.append(int(re.search(r"_state_(\d+)_", e).group(1)))
            A, B, _ = shared(Ps, Pc)
            rows.append({**base, "entry": e, "kind": "state",
                         "state": sids[-1],
                         "n_cols_covered": int((~np.isnan(Ps).any(axis=-1)).sum()),
                         **sc, **score(Ps, Yexp, "exp"),
                         "ic_bits": ic(Ps), "consensus": cons(Ps),
                         "cons_match_crystal": int(sum(x == y for x, y in
                                                       zip(cons(Ps), cons(Pc)) if x != "." and y != "."))})

        if stack:
            S = np.stack(stack)
            Pens = norm(np.nanmean(S, axis=0))
            Pens[np.isnan(S).all(axis=0).any(axis=-1)] = np.nan
            rows.append({**base, "entry": "ENSEMBLE_MEAN", "kind": "ensemble_mean", "state": np.nan,
                         "n_cols_covered": int((~np.isnan(Pens).any(axis=-1)).sum()),
                         **score(Pens, Pc, "crystal"), **score(Pens, Yexp, "exp"),
                         "ic_bits": ic(Pens), "consensus": cons(Pens),
                         "cons_match_crystal": int(sum(x == y for x, y in
                                                       zip(cons(Pens), cons(Pc)) if x != "." and y != "."))})
            np.savez_compressed(os.path.join(OUT, f"pwms_{tf}_s{seed}.npz"),
                                states=S, state_ids=np.array(sids),
                                crystal=Pc, experiment=Yexp, ensemble_mean=Pens)
            if seed == 1:
                bundles[tf] = dict(states=S, state_ids=np.array(sids), crystal=Pc,
                                   experiment=Yexp, ensemble_mean=Pens)

df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUT, "per_state_pwm_vs_crystal.csv"), index=False)

st = df[df.kind == "state"]
summ = (st.groupby(["tf", "pdb", "pwm_label"])
          .agg(n_states=("state", "nunique"), n_seeds=("seed", "nunique"),
               n_cols=("n_overlap_crystal", "median"),
               r_mean=("r_vs_crystal", "mean"), r_sd=("r_vs_crystal", "std"),
               r_min=("r_vs_crystal", "min"), r_max=("r_vs_crystal", "max"),
               mae_mean=("mae_vs_crystal", "mean"), jsd_mean=("jsd_vs_crystal", "mean"),
               ic_mean=("ic_bits", "mean"),
               cons_match_mean=("cons_match_crystal", "mean"),
               r_exp_mean=("r_vs_exp", "mean"), mae_exp_mean=("mae_vs_exp", "mean"))
          .reset_index())
cry = (df[df.kind == "crystal"].groupby("tf")
       .agg(crystal_ic=("ic_bits", "mean"), crystal_r_exp=("r_vs_exp", "mean"),
            crystal_mae_exp=("mae_vs_exp", "mean"), crystal_consensus=("consensus", "first")).reset_index())
ens = (df[df.kind == "ensemble_mean"].groupby("tf")
       .agg(ens_r_crystal=("r_vs_crystal", "mean"), ens_mae_crystal=("mae_vs_crystal", "mean"),
            ens_r_exp=("r_vs_exp", "mean"), ens_mae_exp=("mae_vs_exp", "mean"),
            ens_ic=("ic_bits", "mean"), ens_consensus=("consensus", "first")).reset_index())
summ = summ.merge(cry, on="tf", how="left").merge(ens, on="tf", how="left")
summ["ic_deficit"] = summ.crystal_ic - summ.ic_mean
summ["exp_gain_ensemble"] = summ.ens_r_exp - summ.crystal_r_exp
summ.to_csv(os.path.join(OUT, "per_tf_summary.csv"), index=False)

np.savez_compressed(os.path.join(OUT, "pwm_bundle_seed1.npz"),
                    **{f"{tf}__{k}": v for tf, d in bundles.items() for k, v in d.items()})
with open(os.path.join(OUT, "skipped.json"), "w") as f:
    json.dump([list(map(str, s)) for s in skipped], f, indent=1)

pd.set_option("display.width", 250); pd.set_option("display.max_columns", 60)
print("rows", len(df), "| unique states scored", st.entry.nunique(), "| skipped", len(skipped))
print(summ[["tf","pdb","n_states","n_seeds","n_cols","r_mean","r_sd","r_min","r_max","mae_mean","jsd_mean",
            "ic_mean","crystal_ic","ic_deficit","ens_r_crystal","r_exp_mean","crystal_r_exp","ens_r_exp"]].round(3).to_string(index=False))
print()
print(summ[["tf","crystal_consensus","ens_consensus","cons_match_mean"]].to_string(index=False))
