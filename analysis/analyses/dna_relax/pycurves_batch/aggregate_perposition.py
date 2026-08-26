#!/usr/bin/env python3
"""Per-position DNA shape profiles from the pyCurves ensemble.

For each TF x condition, reads every state JSON and builds a per-helical-position
profile (median +/- IQR across states) for a set of shape features. This localizes
WHERE along the duplex the frozen->relaxed deformation concentrates.

Features:
  - minor_width   (legacy groove: data[level].sub_levels[0].minor_width)
  - inclin, tip, xdisp, ydisp  (curvesplus_base_pair_axis per base pair)

Output: analysis/dna_relax/data/pycurves_perposition.csv
  columns: tf, cond, feature, position, n_states, median, q25, q75
Position is the pyCurves level/bp index (1-based along the duplex).
"""
import os, json, glob, csv, sys
from collections import defaultdict
import statistics as st

BASE = "/project2/rohs_102/shewchuk/TF-conformation"
PYC  = f"{BASE}/analysis/dna_relax/pycurves"
OUT  = f"{BASE}/analysis/dna_relax/data/pycurves_perposition.csv"
TFS  = ["tbp","lef1","engrailed","egr1","ets1","foxa","dux4"]

def legacy_minorwidth_by_level(js):
    """Return {level_index(int): minor_width} from legacy groove data (sub_level 0)."""
    out={}
    try:
        g=js["dataframes"]["groove"]["data"]
    except (KeyError,TypeError):
        return out
    for lvl,entry in g.items():
        try:
            sub=entry["sub_levels"]
            # sub_levels keyed "0","1","2"; take level-0 (canonical)
            s0=sub.get("0") or sub.get(0)
            if s0 and s0.get("minor_width") is not None:
                out[int(lvl)]=float(s0["minor_width"])
        except (KeyError,TypeError,ValueError):
            continue
    return out

def cp_localparams_by_bp(js):
    """Return {bp_index(int): {inclin,tip,xdisp,ydisp}} from curvesplus_base_pair_axis."""
    out={}
    try:
        rows=js["dataframes"]["curvesplus_base_pair_axis"]
    except (KeyError,TypeError):
        return out
    for i,r in enumerate(rows,start=1):
        d={}
        for k in ("inclin","tip","xdisp","ydisp"):
            v=r.get(k)
            if v is not None:
                try: d[k]=float(v)
                except (ValueError,TypeError): pass
        if d: out[i]=d
    return out

# accumulate: values[(tf,cond,feature,position)] = [vals across states]
values=defaultdict(list)

def add_state(tf,cond,legacy_path,cp_path):
    # minor width from legacy
    if os.path.exists(legacy_path):
        try:
            js=json.load(open(legacy_path))
            for lvl,mw in legacy_minorwidth_by_level(js).items():
                values[(tf,cond,"minor_width",lvl)].append(mw)
        except (json.JSONDecodeError,OSError): pass
    # local params from curvesplus
    if os.path.exists(cp_path):
        try:
            js=json.load(open(cp_path))
            for bp,d in cp_localparams_by_bp(js).items():
                for k,v in d.items():
                    values[(tf,cond,k,bp)].append(v)
        except (json.JSONDecodeError,OSError): pass

for tf in TFS:
    d=f"{PYC}/{tf}"
    if not os.path.isdir(d): 
        print(f"  no dir for {tf}", file=sys.stderr); continue
    # crystal
    add_state(tf,"crystal",f"{d}/crystal_legacy.json",f"{d}/crystal_curvesplus.json")
    # frozen + relaxed states
    for cond in ("frozen","relaxed"):
        for lp in sorted(glob.glob(f"{d}/{cond}_state_*_legacy.json")):
            state=lp.split(f"{cond}_state_")[1].split("_legacy")[0]
            cp=f"{d}/{cond}_state_{state}_curvesplus.json"
            add_state(tf,cond,lp,cp)
    print(f"  {tf}: accumulated", file=sys.stderr)

# write medians
rows=[]
for (tf,cond,feat,pos),vals in values.items():
    if len(vals)==0: continue
    sv=sorted(vals)
    n=len(sv)
    med=st.median(sv)
    q25=sv[int(0.25*(n-1))]; q75=sv[int(0.75*(n-1))]
    rows.append((tf,cond,feat,pos,n,round(med,4),round(q25,4),round(q75,4)))
rows.sort(key=lambda r:(r[0],r[2],r[1],r[3]))
os.makedirs(os.path.dirname(OUT),exist_ok=True)
with open(OUT,"w",newline="") as fh:
    w=csv.writer(fh)
    w.writerow(["tf","cond","feature","position","n_states","median","q25","q75"])
    w.writerows(rows)
print(f"wrote {OUT}: {len(rows)} rows")
