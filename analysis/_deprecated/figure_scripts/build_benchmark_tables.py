#!/usr/bin/env python3
"""build_benchmark_tables.py - regenerate the figure_scripts input CSVs from
the Stage-7 eval JSONs. Reconstructed from figure_scripts/README.md spec.

Outputs (into --outdir, default figure_scripts/):
  perseed_perentry.csv    per (tf,arm,dna,seed,entry) metrics + family/selffam (KEEPS seed 0)
  perentry_condition.csv  seed-mean per (tf,arm,dna,entry), seeds s1-s5 only (n_seeds)
  perseed_summary.csv     per (tf,arm,dna,seed) benchmark-mean over 130 entries

Reads output/stage7_eval/id_benchmark_<tf>[_dnarelax].json. Checkpoint names:
  {baseline,augmented}_<tf>_fold0[_dnarelax][_sN]  (bare fold0 = seed 0).
Family per entry from analysis/data/family_annotation.csv (entry->family,motif,tf_name,
global 130-entry map). Each pilot's OWN family (for selffam) from PILOT_FAMILY below.
"""
import os, sys, json, glob, csv
from collections import defaultdict
BASE="/project2/rohs_102/shewchuk/TF-conformation"
EVAL=os.path.join(BASE,"output","stage7_eval")
FAM =os.path.join(BASE,"analysis","data","family_annotation.csv")
METRICS=["pearsonr","spearmanr","auroc","ic_weighted_pcc","mae"]

# pilot -> its own Pfam family label (for selffam flag). Matches _common.pilot_entryfam.
PILOT_FAMILY={'ets1':'ETS','tbp':'TBP / β-saddle','egr1':'C2H2 zinc finger',
  'engrailed':'Homeodomain','foxa':'Forkhead','lef1':'HMG-box','csl':'CSL/RBPJ',
  'err':'Nuclear receptor (Zn)','nfat':'Rel/NF-κB (RHD)','runx':'Runt',
  'hsf':'HSF','irf':'IRF'}

def load_family_map():
    m={}
    for r in csv.DictReader(open(FAM)):
        e=r["entry"]
        if e not in m:
            m[e]=dict(family=r["family"],motif=r["motif"],tf_name=r["tf_name"])
    return m

def parse_ckpt(name):
    # {arm}_{tf}_fold0[_dnarelax][_sN]
    arm="baseline" if name.startswith("baseline") else "augmented"
    dna="relax" if "dnarelax" in name else "frozen"
    seed=0
    if "_s" in name.split("fold0")[-1]:
        tail=name.split("fold0")[-1]
        for tok in tail.split("_"):
            if tok.startswith("s") and tok[1:].isdigit(): seed=int(tok[1:])
    return arm,dna,seed

def main(pilots, outdir):
    fammap=load_family_map()
    perseed_rows=[]   # keeps seed 0
    for tf in pilots:
        own=PILOT_FAMILY.get(tf)
        for suf,dna in [("","frozen"),("_dnarelax","relax")]:
            f=os.path.join(EVAL,f"id_benchmark_{tf}{suf}.json")
            if not os.path.isfile(f): continue
            res=json.load(open(f)).get("results",{})
            for ckpt,entries in res.items():
                arm,dna2,seed=parse_ckpt(ckpt)
                for ekey,mets in entries.items():
                    e=ekey if ekey.endswith(".npz") else ekey+".npz"
                    fam=fammap.get(e,{})
                    row=dict(tf=tf,arm=arm,dna=dna,seed=seed,entry=ekey.replace(".npz",""),
                             family=fam.get("family",""),motif=fam.get("motif",""),
                             tf_name=fam.get("tf_name",""))
                    for m in METRICS: row["m_"+m]=mets.get(m)
                    row["entry_key"]=e
                    row["selffam"]=(fam.get("family")==own)
                    perseed_rows.append(row)
    os.makedirs(outdir,exist_ok=True)
    # 1) perseed_perentry.csv
    cols=["tf","arm","dna","seed","entry","family","motif","tf_name",
          "m_pearsonr","m_spearmanr","m_auroc","m_ic_weighted_pcc","m_mae","entry_key","selffam"]
    with open(os.path.join(outdir,"perseed_perentry.csv"),"w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=cols); w.writeheader()
        for r in perseed_rows: w.writerow({k:r[k] for k in cols})
    # 2) perentry_condition.csv  (seeds s1-s5 only, mean over seeds)
    agg=defaultdict(list)
    for r in perseed_rows:
        if r["seed"]==0: continue
        agg[(r["tf"],r["arm"],r["dna"],r["entry_key"])].append(r)
    with open(os.path.join(outdir,"perentry_condition.csv"),"w",newline="") as fh:
        cols2=["tf","arm","dna","entry","family","motif","tf_name","n_seeds",
               "mean_pearsonr","mean_spearmanr","mean_auroc","mean_ic_weighted_pcc","mean_mae"]
        w=csv.DictWriter(fh,fieldnames=cols2); w.writeheader()
        for (tf,arm,dna,ekey),rs in agg.items():
            def mean(m):
                vals=[x["m_"+m] for x in rs if x["m_"+m] is not None]
                return sum(vals)/len(vals) if vals else ""
            r0=rs[0]
            w.writerow(dict(tf=tf,arm=arm,dna=dna,entry=ekey,family=r0["family"],
                motif=r0["motif"],tf_name=r0["tf_name"],n_seeds=len(rs),
                mean_pearsonr=mean("pearsonr"),mean_spearmanr=mean("spearmanr"),
                mean_auroc=mean("auroc"),mean_ic_weighted_pcc=mean("ic_weighted_pcc"),
                mean_mae=mean("mae")))
    # 3) perseed_summary.csv (per tf/arm/dna/seed benchmark mean; seeds s1-s5)
    agg2=defaultdict(list)
    for r in perseed_rows:
        if r["seed"]==0: continue
        agg2[(r["tf"],r["arm"],r["dna"],r["seed"])].append(r)
    with open(os.path.join(outdir,"perseed_summary.csv"),"w",newline="") as fh:
        cols3=["tf","arm","dna","seed","n_entries","mean_pearsonr","mean_spearmanr",
               "mean_auroc","mean_ic_weighted_pcc","mean_mae"]
        w=csv.DictWriter(fh,fieldnames=cols3); w.writeheader()
        for (tf,arm,dna,seed),rs in agg2.items():
            def mean(m):
                vals=[x["m_"+m] for x in rs if x["m_"+m] is not None]
                return sum(vals)/len(vals) if vals else ""
            w.writerow(dict(tf=tf,arm=arm,dna=dna,seed=seed,n_entries=len(rs),
                mean_pearsonr=mean("pearsonr"),mean_spearmanr=mean("spearmanr"),
                mean_auroc=mean("auroc"),mean_ic_weighted_pcc=mean("ic_weighted_pcc"),
                mean_mae=mean("mae")))
    print(f"wrote 3 CSVs for {len(pilots)} pilots -> {outdir}")

if __name__=="__main__":
    pilots=["ets1","tbp","egr1","engrailed","foxa","lef1","csl","err","nfat","runx","hsf","irf"]
    outdir=sys.argv[1] if len(sys.argv)>1 else os.path.join(BASE,"analysis","figure_scripts")
    main(pilots,outdir)
