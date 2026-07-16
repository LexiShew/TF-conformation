#!/usr/bin/env python3
"""Extract train/val overfitting diagnostics from DeepPBS Model_metrics.json across all stage6_train runs.
Reads output/stage6_train/*/*/Model_metrics.json; writes analysis/dna_relax/data/training_overfit_metrics.csv.
best_epoch = epoch of min validation MAE (the checkpoint used for benchmarking, Model.best.tar).
Key signals: val_rise_after_min (val loss increase from its min to final epoch),
gap_at_best (val-train loss gap at the selected checkpoint)."""
import json, os, re, glob, numpy as np, pandas as pd
ROOT=os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
rows=[]
for f in glob.glob(os.path.join(ROOT,"output/stage6_train/*/*/Model_metrics.json")):
    run=os.path.basename(os.path.dirname(f))
    m=re.match(r"(baseline|augmented)_([a-z0-9]+)_fold0(_dnarelax)?(?:_s(\d+))?$",run)
    if not m: continue
    d=json.load(open(f)); va=np.array(d["validation"]["loss"],dtype=float); tr=np.array(d["train"]["loss"],dtype=float)
    if np.all(np.isnan(va)): continue
    be=d.get("best_epoch"); n=len(tr); bi=min(max((be-1) if be else int(np.nanargmin(va)),0),n-1)
    rows.append(dict(arm=m.group(1),tf=m.group(2),cond="relaxed" if m.group(3) else "frozen",
        seed=int(m.group(4)) if m.group(4) else 0,run=run,n_epochs=n,best_epoch=be,
        val_loss_min=float(np.nanmin(va)),val_loss_final=float(va[-1]),
        val_rise_after_min=float(va[-1]-np.nanmin(va)),gap_at_best=float(va[bi]-tr[bi])))
pd.DataFrame(rows).to_csv(os.path.join(ROOT,"analysis/dna_relax/data/training_overfit_metrics.csv"),index=False)
print(f"wrote {len(rows)} runs")
