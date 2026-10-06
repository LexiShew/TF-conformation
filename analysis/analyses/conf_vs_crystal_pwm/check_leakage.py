#!/usr/bin/env python
"""Two-level training membership for each pilot's crystal entry:
  exact  = same structure AND same PWM label present in the training file
  struct = same PDB+chain present under a DIFFERENT PWM label (structure seen,
           label not) -- still structural leakage, weaker than exact.
"""
import os, re, json
R = "/project2/rohs_102/shewchuk/TF-conformation"
PIL = ["csl","egr1","engrailed","err","ets1","foxa","hsf","irf","lef1","nfat","runx","tbp"]
def meta(tf):
    t = open(f"{R}/config/pilots/{tf}.sh").read()
    g = lambda k: re.search(k+r'="([^"]+)"', t).group(1)
    return g("PDB_ID"), g("PWM_LABEL"), g("BINDING_CHAIN")
def rd(p): return set(l.strip() for l in open(p)) if os.path.exists(p) else set()
base_tr = rd(f"{R}/data/folds/train0.txt")

rows = {}
for tf in PIL:
    pid, pwm, bc = meta(tf)
    cry = "5d5u_duplex_MA0486.2.jaspar.npz" if tf=="hsf" else f"{pid}_{bc}_{pwm}.npz"
    af = rd(f"{R}/output/stage5_aug/folds_aug/train0_aug_{tf}.txt")
    ar = rd(f"{R}/output/stage5_aug/folds_aug/train0_aug_dnarelax_{tf}.txt")
    pref = f"{pid}_{bc}_"
    same_struct_base = sorted(e for e in base_tr if e.startswith(pref))
    same_struct_aug  = sorted(e for e in af      if e.startswith(pref))
    rows[tf] = dict(crystal_entry=cry,
        exact_in_baseline = cry in base_tr,
        exact_in_aug_frozen = cry in af,
        exact_in_aug_relax  = cry in ar,
        same_structure_other_label_baseline = [e for e in same_struct_base if e != cry],
        same_structure_other_label_aug      = [e for e in same_struct_aug  if e != cry])
json.dump(rows, open(f"{R}/analysis/analyses/conf_vs_crystal_pwm/crystal_leakage.json","w"), indent=1)
print(f"{'tf':10s} {'exact_base':11s} {'exact_augF':11s} {'same_struct_other_label_in_base'}")
for tf,d in rows.items():
    o = ",".join(d["same_structure_other_label_baseline"]) or "-"
    print(f'{tf:10s} {str(d["exact_in_baseline"]):11s} {str(d["exact_in_aug_frozen"]):11s} {o}')
ex = [tf for tf,d in rows.items() if d["exact_in_baseline"]]
st = [tf for tf,d in rows.items() if not d["exact_in_baseline"] and d["same_structure_other_label_baseline"]]
cl = [tf for tf,d in rows.items() if not d["exact_in_baseline"] and not d["same_structure_other_label_baseline"]]
print("\nexact-entry leakage :", ex)
print("structure-only leak :", st)
print("fully held out      :", cl)
