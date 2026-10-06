#!/usr/bin/env python
"""Score conformation PWMs against the crystal prediction for all three
DeepPBS conditions: baseline (no augmentation), augmented/frozen-DNA, and
augmented/relaxed-DNA.

Each condition is scored SELF-CONSISTENTLY: a conformation's PWM is compared
to the crystal PWM produced by the SAME checkpoint, so differences reflect the
model, not a shifted reference. Relaxed-arm entries come from the relaxed
feature directory (same entry names, relaxed DNA geometry).

Leakage is recorded per pilot, not assumed: the augmented arms were trained on
these very conformations, so their state-level agreement is IN-SAMPLE by
construction and is reported as such.
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

# condition -> (run-dir glob fold token, feature dir template, id-file template)
COND = {
 "baseline":      ("0",          f"{R}/output/stage5_aug/combined_assembly_{{tf}}",          f"{IDD}/id_{{tf}}_confcmp.txt"),
 "aug_frozen":    ("0",          f"{R}/output/stage5_aug/combined_assembly_{{tf}}",          f"{IDD}/id_{{tf}}_confcmp.txt"),
 "aug_relax":     ("0_dnarelax", f"{R}/output/stage5_aug_dnarelax/combined_assembly_{{tf}}", f"{IDD}/id_{{tf}}_confcmp_relax.txt"),
}
KIND = {"baseline":"baseline", "aug_frozen":"augmented", "aug_relax":"augmented"}

leak = json.load(open(f"{OUT}/crystal_leakage.json"))

def meta(tf):
    t = open(f"{R}/config/pilots/{tf}.sh").read()
    g = lambda k: re.search(k+r'="([^"]+)"', t).group(1)
    return g("PDB_ID"), g("PWM_LABEL"), g("BINDING_CHAIN")

def norm(p):
    s = p.sum(axis=-1, keepdims=True)
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

def shared(a, b):
    m = ~(np.isnan(a).any(-1) | np.isnan(b).any(-1))
    return a[m], b[m], int(m.sum())

def pear(a, b):
    a, b = a.ravel(), b.ravel()
    if a.size < 2 or a.std() == 0 or b.std() == 0: return np.nan
    return float(np.corrcoef(a, b)[0, 1])

def jsd(a, b):
    m = 0.5*(a+b)
    kl = lambda x, y: np.sum(np.where(x > 0, x*np.log2(np.divide(x, y, out=np.ones_like(x), where=y > 0)), 0.0), -1)
    return float(np.mean(0.5*kl(a, m) + 0.5*kl(b, m)))

def ic(p):
    q = p[~np.isnan(p).any(-1)]
    if q.size == 0: return np.nan
    return float(np.sum(2.0 + np.sum(np.where(q > 0, q*np.log2(np.where(q > 0, q, 1)), 0.0), -1)))

def cons(p):
    return "".join("." if np.isnan(r).any() else "ACGT"[int(np.argmax(r))] for r in p)

def sc(a, b, tag):
    A, B, n = shared(a, b)
    if n < MIN_OVERLAP:
        return {f"r_vs_{tag}": np.nan, f"mae_vs_{tag}": np.nan, f"jsd_vs_{tag}": np.nan, f"n_overlap_{tag}": n}
    return {f"r_vs_{tag}": pear(A, B), f"mae_vs_{tag}": float(np.abs(A-B).mean()),
            f"jsd_vs_{tag}": jsd(A, B), f"n_overlap_{tag}": n}

rows, skipped = [], []
for cond, (foldtok, featT, idT) in COND.items():
    kind = KIND[cond]
    for tf in PIL:
        pid, pwm, bc = meta(tf)
        cry = leak[tf]["crystal_entry"]
        idf = idT.format(tf=tf)
        if not os.path.exists(idf):
            skipped.append((cond, tf, "ALL", "no id file")); continue
        entries = [l.strip() for l in open(idf) if l.strip()]
        states = [e for e in entries if "_state_" in e]
        lk = ("exact" if leak[tf]["exact_in_baseline"]
              else "structure_only" if leak[tf]["same_structure_other_label_baseline"] else "held_out")
        for outer in sorted(glob.glob(f"{TRAIN}/{kind}_{tf}_fold{foldtok}_s*")):
            b = os.path.basename(outer)
            if ".incomplete" in b: continue
            if foldtok == "0" and "dnarelax" in b: continue
            m = re.search(r"_s(\d+)$", b)
            if not m: continue
            seed = int(m.group(1))
            pdir = f"{outer}/{b}/predictions_confcmp"
            cp = f"{pdir}/{cry}_predict.npz"
            if not os.path.exists(cp):
                skipped.append((cond, tf, seed, "no crystal prediction")); continue
            zc = np.load(cp)
            Pc, Ye = to_frame(zc), exp_frame(zc)
            base = dict(condition=cond, tf=tf, pdb=pid, seed=seed, crystal_leakage=lk,
                        states_in_own_training=(cond != "baseline"))
            rows.append({**base, "entry": cry, "kind": "crystal", "state": np.nan,
                         **sc(Pc, Pc, "crystal"), **sc(Pc, Ye, "exp"),
                         "ic_bits": ic(Pc), "consensus": cons(Pc)})
            stack, sids = [], []
            for e in states:
                sp = f"{pdir}/{e}_predict.npz"
                if not os.path.exists(sp):
                    skipped.append((cond, tf, seed, f"missing {e}")); continue
                zs = np.load(sp)
                Ps = to_frame(zs)
                if Ps.shape != Pc.shape:
                    skipped.append((cond, tf, seed, f"len {Ps.shape} vs {Pc.shape}")); continue
                s1 = sc(Ps, Pc, "crystal")
                if np.isnan(s1["r_vs_crystal"]):
                    skipped.append((cond, tf, seed, f"overlap {s1['n_overlap_crystal']}")); continue
                stack.append(Ps); sids.append(int(re.search(r"_state_(\d+)_", e).group(1)))
                rows.append({**base, "entry": e, "kind": "state", "state": sids[-1],
                             **s1, **sc(Ps, Ye, "exp"), "ic_bits": ic(Ps), "consensus": cons(Ps)})
            if stack:
                S = np.stack(stack)
                Pe = norm(np.nanmean(S, 0))
                Pe[np.isnan(S).all(0).any(-1)] = np.nan
                rows.append({**base, "entry": "ENSEMBLE_MEAN", "kind": "ensemble_mean", "state": np.nan,
                             **sc(Pe, Pc, "crystal"), **sc(Pe, Ye, "exp"),
                             "ic_bits": ic(Pe), "consensus": cons(Pe)})

df = pd.DataFrame(rows)
df.to_csv(f"{OUT}/per_state_all_conditions.csv", index=False)
json.dump([list(map(str, s)) for s in skipped], open(f"{OUT}/skipped_all_conditions.json","w"), indent=1)

st = df[df.kind == "state"]
summ = (st.groupby(["condition","tf","crystal_leakage"])
        .agg(n_states=("state","nunique"), n_seeds=("seed","nunique"),
             r_mean=("r_vs_crystal","mean"), r_sd=("r_vs_crystal","std"),
             mae_mean=("mae_vs_crystal","mean"), jsd_mean=("jsd_vs_crystal","mean"),
             ic_mean=("ic_bits","mean"), r_exp_mean=("r_vs_exp","mean")).reset_index())
cry = df[df.kind=="crystal"].groupby(["condition","tf"]).agg(
        crystal_ic=("ic_bits","mean"), crystal_r_exp=("r_vs_exp","mean")).reset_index()
ens = df[df.kind=="ensemble_mean"].groupby(["condition","tf"]).agg(
        ens_r_crystal=("r_vs_crystal","mean"), ens_r_exp=("r_vs_exp","mean"),
        ens_ic=("ic_bits","mean")).reset_index()
summ = summ.merge(cry, on=["condition","tf"]).merge(ens, on=["condition","tf"])
summ["ic_deficit"] = summ.crystal_ic - summ.ic_mean
summ.to_csv(f"{OUT}/per_tf_all_conditions.csv", index=False)

pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
print("rows", len(df), "| skipped", len(skipped))
print("\n== state-level agreement with the same model's crystal PWM ==")
piv = summ.pivot_table(index=["tf","crystal_leakage"], columns="condition",
                       values=["r_mean","ic_mean","ens_r_crystal"]).round(3)
print(piv.to_string())
print("\n== condition means (states) ==")
print(st.groupby("condition").agg(n=("entry","size"), r_mean=("r_vs_crystal","mean"),
      r_med=("r_vs_crystal","median"), ic=("ic_bits","mean"), r_exp=("r_vs_exp","mean")).round(3).to_string())
print("\n== crystal PWM agreement with experiment, per condition ==")
print(df[df.kind=="crystal"].groupby("condition").agg(
      r_exp=("r_vs_exp","mean"), ic=("ic_bits","mean")).round(3).to_string())
