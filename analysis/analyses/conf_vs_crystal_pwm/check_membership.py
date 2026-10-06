#!/usr/bin/env python
"""Exact per-entry training membership for every condition, by entry NAME."""
import os, re, json
R = "/project2/rohs_102/shewchuk/TF-conformation"
PIL = ["csl","egr1","engrailed","err","ets1","foxa","hsf","irf","lef1","nfat","runx","tbp"]
def meta(tf):
    t = open(f"{R}/config/pilots/{tf}.sh").read()
    g = lambda k: re.search(k+r'="([^"]+)"', t).group(1)
    return g("PDB_ID"), g("PWM_LABEL"), g("BINDING_CHAIN")
def rd(p):
    return set(l.strip() for l in open(p)) if os.path.exists(p) else set()

base_tr, base_va = rd(f"{R}/data/folds/train0.txt"), rd(f"{R}/data/folds/valid0.txt")
out = {}
for tf in PIL:
    pid, pwm, bc = meta(tf)
    cry = "5d5u_duplex_MA0486.2.jaspar.npz" if tf=="hsf" else f"{pid}_{bc}_{pwm}.npz"
    af = rd(f"{R}/output/stage5_aug/folds_aug/train0_aug_{tf}.txt")
    ar = rd(f"{R}/output/stage5_aug/folds_aug/train0_aug_dnarelax_{tf}.txt")
    av = rd(f"{R}/output/stage5_aug/folds_aug/valid0_{tf}.txt")
    sf = rd(f"{R}/stage7_eval/confcmp_idfiles/id_{tf}_confcmp.txt") - {cry}
    sr = rd(f"{R}/stage7_eval/confcmp_idfiles/id_{tf}_confcmp_relax.txt") - {cry}
    out[tf] = dict(
      crystal_entry=cry,
      crystal_in={"baseline_train":cry in base_tr, "baseline_valid":cry in base_va,
                  "aug_frozen_train":cry in af, "aug_relax_train":cry in ar,
                  "aug_valid":cry in av},
      n_states_frozen=len(sf), n_states_relax=len(sr),
      states_in_aug_frozen_train=len(sf & af), states_in_aug_relax_train=len(sr & ar),
      states_in_baseline_train=len(sf & base_tr),
      frozen_states_heldout_from_frozen=sorted(sf - af)[:6],
      relax_states_heldout_from_relax=sorted(sr - ar)[:6])
json.dump(out, open(f"{R}/analysis/analyses/conf_vs_crystal_pwm/train_membership.json","w"), indent=1)
print(f"{'tf':10s} {'cry_base':9s} {'cry_augF':9s} {'cry_augR':9s} {'nF':>4s} {'inF':>4s} {'nR':>4s} {'inR':>4s}")
for tf,d in out.items():
    c = d["crystal_in"]
    print(f'{tf:10s} {str(c["baseline_train"]):9s} {str(c["aug_frozen_train"]):9s} {str(c["aug_relax_train"]):9s} '
          f'{d["n_states_frozen"]:4d} {d["states_in_aug_frozen_train"]:4d} {d["n_states_relax"]:4d} {d["states_in_aug_relax_train"]:4d}')
print("\nstates NOT in their own arm's training set (frozen):",
      {tf:d["frozen_states_heldout_from_frozen"] for tf,d in out.items() if d["frozen_states_heldout_from_frozen"]})
print("\nany crystal in a validation set?:", {tf:d["crystal_in"] for tf,d in out.items() if d["crystal_in"]["baseline_valid"] or d["crystal_in"]["aug_valid"]})
