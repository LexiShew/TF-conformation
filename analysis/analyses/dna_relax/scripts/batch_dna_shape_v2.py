#!/usr/bin/env python3
"""
Batch DNA-geometry analysis for the DNA-relaxation study (generalized, per-TF).

Generalizes scripts/batch_dna_shape.py (which hardcoded the 1tgh TBP register)
so it runs on ANY of the TF pilots. The antiparallel base-pair register is
auto-detected geometrically from the docked reference: strand-B DNA residues
are paired to their nearest strand-C partner by C1'-C1' distance, keeping only
mutual nearest neighbours within a Watson-Crick distance window. Everything
downstream (RMSD, bend, P-P fraying, per-residue P displacement, per-bp
C1'-C1' width) is identical in definition to the original TBP script.

For every state in the frozen-DNA Stage-3 dir and the relaxed (_dnarelax)
Stage-3 dir, computes -- relative to the docked Stage-2 input:
  - DNA backbone (P, C1') RMSD vs docked (raw, no alignment; protein bb pinned
    by the Stage-3 restraint so DNA motion is measured directly)
  - global bend angle (principal-axis split of the bp-center helix, two halves)
  - max adjacent P-P gap (fraying proxy; ~6-7 A canonical, >~9 A = unwinding)
  - per-residue P displacement (localization of motion)
  - per-base-pair C1'-C1' width (groove/pairing proxy)

Usage:
  python batch_dna_shape_v2.py <docked_dir> <frozen_dir> <relaxed_dir> <out_prefix>
Runs in the deeppbs env (numpy only). CPU. No alignment/superposition.
Set OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 (login-node thread cap segfaults BLAS).
"""
import numpy as np, sys, glob, os, csv, json
DNA={"DA","DG","DC","DT","DA5","DG5","DC5","DT5","DA3","DG3","DC3","DT3"}

def load(path):
    """Return {(chain,resseq,atom): xyz} for DNA P and C1' atoms."""
    d={}
    for ln in open(path):
        if not ln.startswith(("ATOM","HETATM")): continue
        if ln[17:20].strip() not in DNA: continue
        a=ln[12:16].strip()
        if a not in ("P","C1'"): continue
        d[(ln[21],int(ln[22:26]),a)]=np.array([float(ln[30:38]),float(ln[38:46]),float(ln[46:54])])
    return d

def detect_register(dref):
    """Auto-detect the antiparallel base-pair register positionally.

    The two DNA strands are the first two chains (sorted); strand B is paired to
    strand C position-by-position. Duplex orientation (whether strand C runs
    5'->3' with strand B, i.e. ascending, or antiparallel by numbering, i.e.
    descending) is decided geometrically from the strand endpoints: pair B's
    lowest-numbered C1' with whichever C endpoint is closer. This reproduces the
    original hand-set TBP register (B101+k <-> C124-k) exactly and generalizes to
    every equal-length pilot duplex without a distance cutoff that would drop
    widened-groove or terminal pairs.
    Returns [((chB,resB),(chC,resC)), ...] ordered along strand B ascending.
    """
    c1={}
    for (ch,rs,a),xyz in dref.items():
        if a=="C1'": c1[(ch,rs)]=xyz
    chains=sorted({ch for ch,_ in c1})
    if len(chains)<2: return []
    chB,chC=chains[0],chains[1]
    Bres=sorted(rs for (ch,rs) in c1 if ch==chB)
    Cres=sorted(rs for (ch,rs) in c1 if ch==chC)
    if not Bres or not Cres: return []
    xB0=c1[(chB,Bres[0])]
    d_first=np.linalg.norm(xB0-c1[(chC,Cres[0])])
    d_last =np.linalg.norm(xB0-c1[(chC,Cres[-1])])
    Cseq = Cres if d_first<=d_last else list(reversed(Cres))  # antiparallel by numbering => reversed
    n=min(len(Bres),len(Cseq))
    return [((chB,Bres[i]),(chC,Cseq[i])) for i in range(n)]

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

def bend_angle(d,REG):
    cens=[]
    for bp in REG:
        r=bp_center(d,bp)
        if r is None: return float("nan")
        cens.append(r[0])
    cens=np.array(cens)
    if len(cens)<4: return float("nan")
    def paxis(pts):
        pts=pts-pts.mean(0); u,s,vt=np.linalg.svd(pts); v=vt[0]
        if np.dot(pts[-1]-pts[0], v)<0: v=-v
        return v
    h=len(cens)//2
    v1=paxis(cens[:h]); v2=paxis(cens[h:])
    return float(np.degrees(np.arccos(np.clip(np.dot(v1,v2),-1,1))))

def maxppgap(d):
    gaps=[]
    for ch in sorted({c for (c,r,a) in d if a=="P"}):
        rs=sorted(r for (c,r,a) in d if c==ch and a=="P")
        for i in range(len(rs)-1):
            k1=(ch,rs[i],"P"); k2=(ch,rs[i+1],"P")
            if k1 in d and k2 in d: gaps.append(np.linalg.norm(d[k1]-d[k2]))
    return float(max(gaps)) if gaps else float("nan")

def main():
    docked_dir,frozen_dir,relaxed_dir,outpref=sys.argv[1:5]
    ref_states=sorted(glob.glob(f"{docked_dir}/*_state_*.pdb"))
    if not ref_states:
        print("ERROR: no docked states in",docked_dir,file=sys.stderr); sys.exit(2)
    REG=detect_register(load(ref_states[0]))
    reg_bres=[bp[0][1] for bp in REG]
    print(f"register: n_bp={len(REG)} chB={REG[0][0][0] if REG else '?'} "
          f"chC={REG[0][1][0] if REG else '?'} Bres={reg_bres}")
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
            bd=bend_angle(dref,REG); bn=bend_angle(d,REG)
            rows.append(dict(cond=cond,state=st,dna_bb_rmsd=round(r_bb,4),
                             bend_docked=round(bd,2),bend=round(bn,2),dbend=round(bn-bd,2),
                             maxppgap=round(maxppgap(d),3)))
            for k,v in per_res_P_disp(dref,d).items(): perres[cond].setdefault(k,[]).append(v)
            for bp in REG:
                r=bp_center(d,bp)
                if r: perbp[cond].setdefault(bp[0][1],[]).append(r[1])
    if not rows:
        print("ERROR: no matched states",file=sys.stderr); sys.exit(3)
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
    with open(f"{outpref}_register.json","w") as fh:
        json.dump({"n_bp":len(REG),
                   "chB":REG[0][0][0] if REG else None,
                   "chC":REG[0][1][0] if REG else None,
                   "Bres":reg_bres}, fh)
    for cond in ("frozen","relaxed"):
        r=[x["dna_bb_rmsd"] for x in rows if x["cond"]==cond]
        b=[x["dbend"] for x in rows if x["cond"]==cond]
        g=[x["maxppgap"] for x in rows if x["cond"]==cond]
        if not r: continue
        print(f"{cond}: n={len(r)} | DNAbbRMSD med={np.median(r):.3f} mean={np.mean(r):.3f} max={np.max(r):.3f} "
              f"| dBend med={np.median(b):+.2f} | maxPPgap med={np.median(g):.2f} max={np.max(g):.2f} "
              f"| >8A tail={sum(1 for x in g if x>8.0)} >9A={sum(1 for x in g if x>9.0)}")

if __name__=="__main__":
    main()
