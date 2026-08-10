#!/usr/bin/env python3
"""interface_mgwfl.py — see interface_mgwfl.README.py for full method docs.
Run under `deeppbs` conda env from the TF-conformation repo root:
    python analysis/dna_relax/interface_mgwfl.py
Reproduces interface-restricted MGW-FL and its correlation with DeepPBS
augmentation accuracy (Jiang et al. Biophys J 2026 metric on our pilots).
"""
import json, glob, os, sys, subprocess, re
import statistics as st
sys.path.insert(0, "fnat_gate")
import interface_rmsd as ir
from Bio.PDB import MMCIFParser

ROOT="/project2/rohs_102/shewchuk/TF-conformation"
os.chdir(ROOT)
DNA={'DA':'A','DC':'C','DG':'G','DT':'T'}
PILOTS={"egr1":"1aay","engrailed":"3hdd","ets1":"1k79","foxa":"1vtn","lef1":"2lef",
        "tbp":"1tgh","csl":"3brg","err":"1lo1","nfat":"1a66","runx":"1hjc","hsf":"5d5u","irf":"1if1"}
CONTACT_CUT, IFACE_CUT, GAP = 4.5, 5.0, 4
p=MMCIFParser(QUIET=True)

def cfg(tf,key):
    r=subprocess.run(f"grep -E 'export {key}=' config/pilots/{tf}.sh",shell=True,capture_output=True,text=True).stdout
    m=re.search(rf'{key}="?([0-9,]+)"?',r); return m.group(1) if m else None

def contacting_dna(tf,pid):
    ref=f"structures/source_chains/{pid}_chains/{pid}.cif"
    pc=cfg(tf,"PROTEIN_CHAIN"); dc=cfg(tf,"DNA_CHAINS")
    pc=int(pc) if pc is not None else None
    dc=[int(x) for x in dc.split(",")] if dc else None
    R=ir.ref_side(ir.load_models(ref)[0],IFACE_CUT,CONTACT_CUT,GAP,protein_chain=pc,dna_chains=dc)
    rd=R["rd"]; cj=sorted(set(j for i,j in R["native"]))
    contacts=set((rd[j].get_parent().id,rd[j].id[1]) for j in cj)
    chains=sorted(set(r.get_parent().id for r in rd))
    return contacts, chains

def crystal_strand(pid,chain):
    st=p.get_structure("x",f"structures/source_chains/{pid}_chains/{pid}.cif")
    for m in st:
        for ch in m:
            if ch.id!=chain: continue
            seq=[];nums=[]
            for r in ch:
                rn=r.get_resname().strip()
                if rn in DNA: seq.append(DNA[rn]); nums.append(r.id[1])
            return "".join(seq),nums
    return "",[]

def mgw_vec(path, want="fl"):
    d=json.load(open(path)); g=d.get("dataframes",{}).get("groove")
    if not g or "data" not in g: return None
    lv=sorted(g["data"].keys(),key=lambda x:int(x)); prof=[]
    for l in lv:
        subs=g["data"][l].get("sub_levels",{})
        vals=[subs[s].get("minor_width") for s in subs if isinstance(subs[s].get("minor_width"),(int,float))]
        prof.append(sum(vals)/len(vals) if vals else None)
    return prof

def ensemble_fl(paths):
    profs=[q for q in (mgw_vec(f) for f in paths) if q]
    if not profs: return None
    L=min(len(q) for q in profs); profs=[q[:L] for q in profs]
    return [round(st.pstdev([q[i] for q in profs if q[i] is not None]),3)
            if len([q[i] for q in profs if q[i] is not None])>1 else None for i in range(L)]

def bp_seq(path):
    d=json.load(open(path)); g=d["dataframes"]["groove"]["data"]
    lv=sorted(g.keys(),key=lambda x:int(x))
    return "".join(g[l]["base_pair"].split("-")[0] for l in lv)

def iface_levels(tf,pid,contacts,chains):
    cf=glob.glob(f"analysis/dna_relax/pycurves/{tf}/crystal_legacy.json")
    if not cf: return None
    pyc=bp_seq(cf[0]); s1=chains[0]; cseq,cnums=crystal_strand(pid,s1)
    idx=cseq.find(pyc)
    if idx<0:  # best-offset fallback
        best=(-1,-1)
        for off in range(max(1,len(cseq)-len(pyc)+1)):
            mt=sum(1 for a,b in zip(pyc,cseq[off:off+len(pyc)]) if a==b)
            if mt>best[1]: best=(off,mt)
        idx=best[0]
    resmap=[cnums[idx+k] if 0<=idx+k<len(cnums) else None for k in range(len(pyc))]
    return [i for i,rn in enumerate(resmap) if rn is not None and (s1,rn) in contacts]

def bench(tf):
    f=f"output/stage7_eval/id_benchmark_{tf}.json"
    if not os.path.exists(f): return None
    res=json.load(open(f)).get("results",{})
    def arm(w):
        P=[];M=[]
        for model,structs in res.items():
            m=model.lower()
            if 'dnarelax' in m or 'legacy' in m or not m.startswith(w): continue
            ps=[v.get("pearsonr") for v in structs.values() if isinstance(v.get("pearsonr"),(int,float))]
            ms=[v.get("mae") for v in structs.values() if isinstance(v.get("mae"),(int,float))]
            if ps: P.append(sum(ps)/len(ps))
            if ms: M.append(sum(ms)/len(ms))
        return (sum(P)/len(P) if P else None, sum(M)/len(M) if M else None)
    bP,bM=arm("baseline"); aP,aM=arm("augmented")
    return dict(base_P=bP,aug_P=aP,base_MAE=bM,aug_MAE=aM)

if __name__=="__main__":
    rows=[]
    for tf,pid in PILOTS.items():
        try:
            contacts,chains=contacting_dna(tf,pid)
            idx=iface_levels(tf,pid,contacts,chains) or []
            froz=ensemble_fl(sorted(glob.glob(f"analysis/dna_relax/pycurves/{tf}/frozen_state_*_legacy.json")))
            relax=ensemble_fl(sorted(glob.glob(f"analysis/dna_relax/pycurves/{tf}/relaxed_state_*_legacy.json")))
            af3=ensemble_fl(sorted(glob.glob(f"af3/af3_dna/{tf}_{pid}/*_legacy.json")))
            def im(v): 
                if not v: return None
                x=[v[i] for i in idx if i<len(v) and v[i] is not None]; return sum(x)/len(x) if x else None
            b=bench(tf) or {}
            rows.append(dict(pilot=tf,n_iface=len(idx),
                iface_mgwfl_froz=im(froz),iface_mgwfl_relax=im(relax),iface_mgwfl_af3=im(af3),
                **{k:b.get(k) for k in ("base_P","aug_P","base_MAE","aug_MAE")}))
        except Exception as e:
            rows.append(dict(pilot=tf,error=str(e)[:100]))
    import csv
    with open("analysis/dna_relax/data/iface_mgwfl_vs_accuracy.csv","w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=sorted({k for r in rows for k in r})); w.writeheader(); w.writerows(rows)
    print("wrote analysis/dna_relax/data/iface_mgwfl_vs_accuracy.csv", len(rows),"pilots")
