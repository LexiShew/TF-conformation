#!/usr/bin/env python3
"""
Batch DNA-geometry analysis for the DNA-relaxation study.

For every state in a frozen-DNA Stage-3 dir and a relaxed (_dnarelax) Stage-3 dir,
compute, relative to the docked Stage-2 input:
  - DNA backbone (P, C1') RMSD vs docked (raw, no alignment: protein bb is pinned
    to absolute coords by the Stage-3 restraint, so DNA motion is measured directly)
  - global bend angle (principal-axis split of the base-pair-center helix, two halves)
  - max adjacent P-P gap (fraying proxy; canonical B-DNA ~6-7 A, >~9 A = unwinding)
  - per-residue P displacement (localization of motion)
  - per-base-pair C1'-C1' width (groove/pairing proxy)

Register is FIXED antiparallel: chain B (101+k) <-> chain C (124-k), k=0..11
(the 1tgh TBP TATA duplex). Adjust REG for other constructs.

Usage:
  python batch_dna_shape.py <docked_dir> <frozen_dir> <relaxed_dir> <out_prefix>
Runs in the deeppbs env (numpy only). CPU. No alignment/superposition.
"""
import numpy as np, sys, glob, os, csv
DNA={"DA","DG","DC","DT","DA5","DG5","DC5","DT5","DA3","DG3","DC3","DT3"}
def load(path):
    d={}
    for ln in open(path):
        if not ln.startswith(("ATOM","HETATM")): continue
        if ln[17:20].strip() not in DNA: continue
        a=ln[12:16].strip()
        if a not in ("P","C1'"): continue
        d[(ln[21],int(ln[22:26]),a)]=np.array([float(ln[30:38]),float(ln[38:46]),float(ln[46:54])])
    return d
REG=[(("B",101+k),("C",124-k)) for k in range(12)]
def rmsd_vs(dref,d,atoms):
    v=[np.linalg.norm(d[a]-dref[a]) for a in atoms if a in d and a in dref]
    return float(np.sqrt(np.mean(np.square(v)))) if v else float("nan")
def per_res_P_disp(dref,d):
    out={}
    for (ch,rs,a) in dref:
        if a!="P": continue
        k=(ch,rs,a)
        if k in d: out[(ch,rs)]=float(np.linalg.norm(d[k]-dref[k]))
    return out
def bp_center(d,bp):
    (cA,rA),(cB,rB)=bp
    a=d.get((cA,rA,"C1'")); b=d.get((cB,rB,"C1'"))
    if a is None or b is None: return None
    return (a+b)/2.0, np.linalg.norm(a-b)
def bend_angle(d):
    cens=[]
    for bp in REG:
        r=bp_center(d,bp)
        if r is None: return float("nan")
        cens.append(r[0])
    cens=np.array(cens)
    def paxis(pts):
        pts=pts-pts.mean(0); u,s,vt=np.linalg.svd(pts); v=vt[0]
        if np.dot(pts[-1]-pts[0], v)<0: v=-v
        return v
    v1=paxis(cens[:6]); v2=paxis(cens[6:])
    return float(np.degrees(np.arccos(np.clip(np.dot(v1,v2),-1,1))))
def maxppgap(d):
    gaps=[]
    for ch in ("B","C"):
        rs=sorted(r for (c,r,a) in d if c==ch and a=="P")
        for i in range(len(rs)-1):
            k1=(ch,rs[i],"P"); k2=(ch,rs[i+1],"P")
            if k1 in d and k2 in d: gaps.append(np.linalg.norm(d[k1]-d[k2]))
    return float(max(gaps)) if gaps else float("nan")
docked_dir,frozen_dir,relaxed_dir,outpref=sys.argv[1:5]
rows=[]; perres={"frozen":{}, "relaxed":{}}; perbp={"frozen":{}, "relaxed":{}}
for cond,ddir in [("frozen",frozen_dir),("relaxed",relaxed_dir)]:
    for f in sorted(glob.glob(f"{ddir}/*_state_*.pdb")):
        st=os.path.basename(f).split("_state_")[1].split(".")[0]
        pdbid=os.path.basename(f).split("_state_")[0]
        drefp=f"{docked_dir}/{pdbid}_state_{st}.pdb"
        if not os.path.exists(drefp): continue
        dref=load(drefp); d=load(f)
        atoms=[k for k in sorted(dref) if k[2] in ("P","C1'") and k in d]
        r_bb=rmsd_vs(dref,d,atoms)
        bd=bend_angle(dref); bn=bend_angle(d)
        rows.append(dict(cond=cond,state=st,dna_bb_rmsd=round(r_bb,4),
                         bend_docked=round(bd,2),bend=round(bn,2),dbend=round(bn-bd,2),
                         maxppgap=round(maxppgap(d),3)))
        for k,v in per_res_P_disp(dref,d).items(): perres[cond].setdefault(k,[]).append(v)
        for bp in REG:
            r=bp_center(d,bp)
            if r: perbp[cond].setdefault(bp[0][1],[]).append(r[1])
with open(f"{outpref}_perstate.csv","w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
with open(f"{outpref}_perres.csv","w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["chain","resseq","cond","mean_Pdisp","sd_Pdisp","n"])
    for cond in ("frozen","relaxed"):
        for (ch,rs),vals in sorted(perres[cond].items()):
            w.writerow([ch,rs,cond,round(np.mean(vals),4),round(np.std(vals),4),len(vals)])
with open(f"{outpref}_perbp.csv","w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["bp_Bres","cond","mean_c1c1","sd_c1c1","n"])
    for cond in ("frozen","relaxed"):
        for br,vals in sorted(perbp[cond].items()):
            w.writerow([br,cond,round(np.mean(vals),4),round(np.std(vals),4),len(vals)])
for cond in ("frozen","relaxed"):
    r=[x["dna_bb_rmsd"] for x in rows if x["cond"]==cond]
    b=[x["dbend"] for x in rows if x["cond"]==cond]
    g=[x["maxppgap"] for x in rows if x["cond"]==cond]
    print(f"{cond}: n={len(r)} | DNAbbRMSD med={np.median(r):.3f} mean={np.mean(r):.3f} max={np.max(r):.3f} "
          f"| dBend med={np.median(b):+.2f} | maxPPgap med={np.median(g):.2f} max={np.max(g):.2f} "
          f"| >8A tail={sum(1 for x in g if x>8.0)} >9A={sum(1 for x in g if x>9.0)}")
