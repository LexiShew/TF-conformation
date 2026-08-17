#!/usr/bin/env python3
"""Parse Stage-3 minimization slurm logs into tidy CSVs for the Ch.9 figure suite.

Reads slurm_output/min_*.out (+ min_recover*.out), classifies each by variant
(frozen / dnarelax), dedups per (variant,pilot,state) preferring the successful
record, and emits:
  summary.csv    one row per state (init/final clashes, collapse stage, runtime,
                 metals, clusters, cage drift, final PE, min distances)
  traj_sample.csv  sampled per-stage trajectories (<=60 states/pilot) for faint-line plots
  quantiles.csv    per (variant,pilot,point) clash/PE quantiles over ALL states
  recovery.csv     per pilot: states recovered / still failed by the gentle ramp
"""
import glob, re, os, sys, csv, statistics as st
SO   = "/project2/rohs_102/shewchuk/TF-conformation/slurm_output"
OUT  = "/project2/rohs_102/shewchuk/TF-conformation/analysis/stage3_figs/data"
os.makedirs(OUT, exist_ok=True)

re_pilot = re.compile(r"Loaded pilot config: (\w+)")
re_tag   = re.compile(r"\[([0-9a-zA-Z]+)_state_(\d+)\]")
re_init  = re.compile(r"Initial heavy-atom clashes: (\d+), min_dist: ([\d.]+)")
re_stage = re.compile(r"Stage (\d+)/(\d+): [σs].*?=([\d.]+), PE=([-\d.eE+]+), clashes=(\d+), min_d=([\d.]+)")
re_final = re.compile(r"Final: PE=([-\d.eE+]+), clashes=(\d+), min_d=([\d.]+)")
re_done  = re.compile(r"DONE in ([\d.]+)s")
re_metal = re.compile(r"Found (\d+) structural metal ion\(s\); (\d+) coordination cluster")
re_drift = re.compile(r"Cage drift.*?mean=([\d.]+), max=([\d.]+)")
re_hatoms= re.compile(r"After H addition: (\d+) atoms")
re_recout= re.compile(r"state (\d+) (RECOVERED|STILL FAILED)")

def classify(txt):
    if "soft tether" in txt or "k_dna" in txt or "DNA released" in txt or "FULLY RELAXED" in txt:
        return "dnarelax"
    if "ignore-metals set" in txt or "legacy mode" in txt:
        return "cage_off"
    return "frozen"

records = {}   # key (variant,pilot,pdb,state) -> dict
pdb2pilot = {}
rec_outcomes = {}  # (pilot,state)->outcome from recovery logs

files = sorted(glob.glob(f"{SO}/min_*.out"))
n_parsed = 0
for f in files:
    try:
        txt = open(f, errors="replace").read()
    except Exception:
        continue
    mt = re_tag.search(txt)
    if not mt:
        continue
    pdb, state = mt.group(1).lower(), int(mt.group(2))
    mp = re_pilot.search(txt)
    pilot = mp.group(1) if mp else pdb2pilot.get(pdb)
    if mp: pdb2pilot.setdefault(pdb, pilot)
    if pilot is None:
        continue
    variant = classify(txt)
    # recovery outcome
    for m in re_recout.finditer(txt):
        rec_outcomes[(pilot, int(m.group(1)))] = (variant, m.group(2))
    # jobid from filename for tie-break
    mj = re.search(r"min(?:_recover(?:_arr)?)?_(\d+)_", os.path.basename(f))
    jobid = int(mj.group(1)) if mj else 0
    stages = [(int(a), float(sg), float(pe), int(cl), float(md))
              for (a, _tot, sg, pe, cl, md) in re_stage.findall(txt)]
    mi = re_init.search(txt); mf = re_final.search(txt)
    md_ = re_done.search(txt); mm = re_metal.search(txt); mdr = re_drift.search(txt)
    success = mf is not None
    rec = dict(variant=variant, pilot=pilot, pdb=pdb, state=state, jobid=jobid,
               success=success,
               init_clash=int(mi.group(1)) if mi else None,
               init_min_d=float(mi.group(2)) if mi else None,
               final_pe=float(mf.group(1)) if mf else None,
               final_clash=int(mf.group(2)) if mf else None,
               final_min_d=float(mf.group(3)) if mf else None,
               runtime_s=float(md_.group(1)) if md_ else None,
               n_metals=int(mm.group(1)) if mm else 0,
               n_clusters=int(mm.group(2)) if mm else 0,
               drift_mean=float(mdr.group(1)) if mdr else None,
               drift_max=float(mdr.group(2)) if mdr else None,
               stages=stages)
    key = (variant, pilot, state)
    old = records.get(key)
    # prefer successful; then higher jobid
    if old is None or (rec["success"] and not old["success"]) or \
       (rec["success"] == old["success"] and jobid >= old["jobid"]):
        records[key] = rec
    n_parsed += 1

# collapse stage = first ramp stage index where clashes drop below max(2, 5% of init)
def collapse_stage(rec):
    if not rec["stages"] or not rec["init_clash"]:
        return None
    thr = max(2, 0.05 * rec["init_clash"])
    for (idx, sg, pe, cl, md) in rec["stages"]:
        if cl <= thr:
            return idx
    return None

# ---- summary.csv ----
with open(f"{OUT}/summary.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["variant","pilot","pdb","state","success","init_clash","final_clash",
                "collapse_stage","runtime_s","n_metals","n_clusters","drift_mean",
                "drift_max","final_pe","init_min_d","final_min_d"])
    for rec in records.values():
        w.writerow([rec["variant"],rec["pilot"],rec["pdb"],rec["state"],int(rec["success"]),
                    rec["init_clash"],rec["final_clash"],collapse_stage(rec),rec["runtime_s"],
                    rec["n_metals"],rec["n_clusters"],rec["drift_mean"],rec["drift_max"],
                    rec["final_pe"],rec["init_min_d"],rec["final_min_d"]])

# ---- trajectory: build point series (init, H-min, stage1..5, final) ----
def series(rec):
    """list of (point_idx,label,sigma,clashes,pe,min_d)"""
    if rec["init_clash"] is None or not rec["stages"]:
        return None
    pts = [(0,"init",None,rec["init_clash"],None,rec["init_min_d"]),
           (1,"H-min",None,rec["init_clash"],None,rec["init_min_d"])]
    for (idx,sg,pe,cl,md) in rec["stages"]:
        pts.append((1+idx, f"s{sg}", sg, cl, pe, md))
    if rec["success"]:
        pts.append((2+len(rec["stages"]),"final",None,rec["final_clash"],rec["final_pe"],rec["final_min_d"]))
    return pts

# sample <=60 states/pilot for faint lines (frozen + dnarelax)
from collections import defaultdict
bypilot = defaultdict(list)
for rec in records.values():
    if rec["success"] and series(rec):
        bypilot[(rec["variant"],rec["pilot"])].append(rec)
with open(f"{OUT}/traj_sample.csv","w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["variant","pilot","state","point","label","sigma","clashes","pe","min_d"])
    for (variant,pilot),recs in bypilot.items():
        recs=sorted(recs,key=lambda r:r["state"])[:60]
        for rec in recs:
            for (pi,lab,sg,cl,pe,md) in series(rec):
                w.writerow([variant,pilot,rec["state"],pi,lab,sg,cl,pe,md])

# ---- quantiles over ALL states, per (variant,pilot,point) ----
agg=defaultdict(lambda: defaultdict(lambda: {"clash":[], "pe":[], "md":[]}))
labels_order={}
for rec in records.values():
    s=series(rec)
    if not s: continue
    for (pi,lab,sg,cl,pe,md) in s:
        d=agg[(rec["variant"],rec["pilot"])][pi]
        if cl is not None: d["clash"].append(cl)
        if pe is not None: d["pe"].append(pe)
        if md is not None: d["md"].append(md)
        labels_order[pi]=lab
def q(xs,p):
    if not xs: return ""
    xs=sorted(xs); import math
    k=(len(xs)-1)*p; f=math.floor(k); ceil=math.ceil(k)
    if f==ceil: return xs[int(k)]
    return xs[f]*(ceil-k)+xs[ceil]*(k-f)
with open(f"{OUT}/quantiles.csv","w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["variant","pilot","point","label","n","clash_q25","clash_med","clash_q75","pe_med","md_med"])
    for (variant,pilot),pts in agg.items():
        for pi in sorted(pts):
            d=pts[pi]
            w.writerow([variant,pilot,pi,labels_order.get(pi,""),len(d["clash"]),
                        q(d["clash"],.25),q(d["clash"],.5),q(d["clash"],.75),
                        q(d["pe"],.5),q(d["md"],.5)])

# ---- recovery.csv ----
recpilot=defaultdict(lambda:{"recovered":0,"failed":0})
for (pilot,state),(variant,outcome) in rec_outcomes.items():
    if outcome=="RECOVERED": recpilot[pilot]["recovered"]+=1
    else: recpilot[pilot]["failed"]+=1
with open(f"{OUT}/recovery.csv","w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["pilot","n_recovered","n_still_failed","n_needed_recovery"])
    for pilot,d in sorted(recpilot.items()):
        w.writerow([pilot,d["recovered"],d["failed"],d["recovered"]+d["failed"]])

# diagnostics
from collections import Counter
vc=Counter((r["variant"] for r in records.values()))
pc=Counter((r["pilot"] for r in records.values() if r["variant"]=="frozen"))
print("parsed_files",n_parsed,"unique_states",len(records),"variants",dict(vc))
print("frozen_states_per_pilot",dict(sorted(pc.items())))
print("recovery_pilots",dict(recpilot))
print("wrote", os.listdir(OUT))
