#!/usr/bin/env python3
"""
Render the 4-panel DNA-shape figure for one TF pilot, matching the TBP figure
(figures/tbp_dna_shape.png). Reads the CSVs produced by batch_dna_shape_v2.py.

Panels (same layout/definitions as the TBP figure):
  a  DNA backbone RMSD vs docked        (jitter strip, frozen grey / relaxed blue, median bars)
  b  per-residue mean P displacement    (line + SD band, strand split marker)
  c  max adjacent P-P gap               (jitter strip, 9 A helix-unwind line, >8 A tail note)
  d  delta global bend angle vs docked  (jitter strip, zero line, median bars)

Panel titles are data-driven (computed per TF), so each panel states a claim
tested against that TF's own data rather than the TBP-specific text.

Usage:
  python plot_dna_shape.py <data_prefix> <tf_label> <out_png>
  e.g. python plot_dna_shape.py data/egr1_dna EGR1 figures/egr1_dna_shape.png
"""
import sys, csv, json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import sys as _sys, os as _os
_sys.path.insert(0, _os.path.abspath(_os.path.join(_os.path.dirname(__file__), "..", "..", "..")))
from palette import TEAL as _FROZ, GREEN as _RELAX, ALARM as RED
GREY=_FROZ; BLUE=_RELAX  # frozen=TEAL, relaxed=GREEN (names kept for downstream use)
CONDS=[("frozen",GREY,"frozen"),("relaxed",BLUE,"DNA-relaxed")]

def read_perstate(p):
    d={"frozen":{}, "relaxed":{}}
    for row in csv.DictReader(open(p)):
        c=row["cond"]
        for k in ("dna_bb_rmsd","dbend","maxppgap"):
            v=row[k]
            d[c].setdefault(k,[]).append(float(v) if v not in ("","nan") else np.nan)
    return {c:{k:np.array(v,dtype=float) for k,v in dd.items()} for c,dd in d.items()}

def read_perres(p):
    d={"frozen":[], "relaxed":[]}
    for row in csv.DictReader(open(p)):
        d[row["cond"]].append((row["chain"],int(row["resseq"]),
                               float(row["mean_Pdisp"]),float(row["sd_Pdisp"])))
    return d

def jitter(ax,x0,vals,color,seed=0):
    vals=vals[np.isfinite(vals)]
    rng=np.random.default_rng(seed)
    xs=x0+rng.uniform(-0.13,0.13,size=len(vals))
    ax.scatter(xs,vals,s=16,color=color,alpha=0.55,linewidths=0,zorder=2)

def median_bar(ax,x0,vals,color,fmt="med {:.2f}",side="right",unit=""):
    vals=vals[np.isfinite(vals)]
    if not len(vals): return np.nan
    m=float(np.median(vals))
    ax.plot([x0-0.22,x0+0.22],[m,m],color=color,lw=3,zorder=3,solid_capstyle="round")
    if fmt:
        dx=0.28 if side=="right" else -0.28
        ha="left" if side=="right" else "right"
        ax.annotate(fmt.format(m)+unit,(x0+dx,m),va="center",ha=ha,color=color,
                    fontsize=10,fontweight="bold")
    return m

def main():
    pref,label,out=sys.argv[1:4]
    ps=read_perstate(f"{pref}_perstate.csv")
    pr=read_perres(f"{pref}_perres.csv")
    reg=json.load(open(f"{pref}_register.json")) if os.path.exists(f"{pref}_register.json") else {}
    def finite(cond,key): 
        v=ps[cond].get(key,np.array([])); return v[np.isfinite(v)]
    nf=len(finite("frozen","dna_bb_rmsd")); nr=len(finite("relaxed","dna_bb_rmsd"))

    fig,axes=plt.subplots(2,2,figsize=(15.5,11))
    (a,b),(cc,dd)=axes
    fig.suptitle(f"{label} DNA minimization: soft-tether relaxation vs frozen-DNA baseline "
                 f"({nr} vs {nf} states)",fontsize=15,y=0.985)

    # ---- panel a: DNA backbone RMSD ----
    mf=np.median(finite("frozen","dna_bb_rmsd")) if nf else np.nan
    mr=np.median(finite("relaxed","dna_bb_rmsd")) if nr else np.nan
    xt=[]
    for i,(cond,col,lab) in enumerate(CONDS):
        v=ps[cond].get("dna_bb_rmsd",np.array([]))
        jitter(a,i,v,col,seed=1+i); median_bar(a,i,v,col,side="right")
        xt.append(f"{lab}\n(n={len(finite(cond,'dna_bb_rmsd'))})")
    a.set_xticks([0,1]); a.set_xticklabels(xt)
    a.set_xlim(-0.5,1.6); a.set_ylabel("DNA backbone RMSD vs docked (\u00c5)")
    ratio=(mr/mf) if (np.isfinite(mf) and mf>0) else np.nan
    ta=(f"DNA moves ~{ratio:.1f}\u00d7 more under relaxation" if np.isfinite(ratio) and ratio>=1.05
        else ("Relaxed and frozen DNA move comparably" if np.isfinite(ratio) else "DNA backbone motion vs docked"))
    a.set_title(ta,loc="center",fontsize=12)
    a.set_title("a",loc="left",fontweight="bold",fontsize=15)

    # ---- panel b: per-residue P displacement ----
    base=pr["frozen"] if pr["frozen"] else pr["relaxed"]
    keyseq=[(ch,rs) for (ch,rs,_,_) in base]
    xidx={k:i for i,k in enumerate(keyseq)}
    labels=[f"{k[0]}{k[1]}" for k in keyseq]
    chains=[k[0] for k in keyseq]; splitx=None
    for i in range(1,len(chains)):
        if chains[i]!=chains[i-1]: splitx=i-0.5; break
    for cond,col,lab in CONDS:
        pts=pr[cond]
        if not pts: continue
        xs=[xidx[(ch,rs)] for (ch,rs,_,_) in pts if (ch,rs) in xidx]
        mu=np.array([m for (ch,rs,m,s) in pts if (ch,rs) in xidx])
        sd=np.array([s for (ch,rs,m,s) in pts if (ch,rs) in xidx])
        b.plot(xs,mu,"-o",color=col,ms=4,lw=1.6,label=lab,zorder=3)
        b.fill_between(xs,mu-sd,mu+sd,color=col,alpha=0.15,zorder=1)
    b.set_ylim(bottom=0)
    if splitx is not None:
        b.axvline(splitx,ls=":",color="0.6",lw=1)
        chB=reg.get("chB","B"); chC=reg.get("chC","C")
        ytop=b.get_ylim()[1]
        b.text(splitx*0.5,ytop*0.96,f"strand {chB}",ha="center",color="0.5",fontsize=10)
        b.text((splitx+len(keyseq))*0.5,ytop*0.96,f"strand {chC}",ha="center",color="0.5",fontsize=10)
    b.set_xticks(range(len(labels))); b.set_xticklabels(labels,rotation=90,fontsize=7)
    b.set_xlabel("DNA residue"); b.set_ylabel("mean P displacement vs docked (\u00c5)")
    b.legend(frameon=False,loc="upper right")
    b.set_title("Per-residue displacement localization",loc="center",fontsize=12)
    b.set_title("b",loc="left",fontweight="bold",fontsize=15)

    # ---- panel c: max adjacent P-P gap ----
    xt=[]; tail8={}; tail9={}
    for i,(cond,col,lab) in enumerate(CONDS):
        v=ps[cond].get("maxppgap",np.array([]))
        jitter(cc,i,v,col,seed=5+i); median_bar(cc,i,v,col,fmt="",side="right")
        vf=finite(cond,"maxppgap")
        tail8[cond]=int(np.sum(vf>8.0)); tail9[cond]=int(np.sum(vf>=9.0))
        xt.append(lab)
    # ensure headroom above the 9 A line so the label never collides with the title
    lo,hi=cc.get_ylim()
    cc.set_ylim(lo,max(hi,9.45))
    cc.axhline(9.0,ls="--",color=RED,lw=1.5)
    # red label just ABOVE the line, tail note just BELOW it -> never overlap each other
    cc.text(0.5,9.05,"helix-unwind regime (\u22659 \u00c5)",ha="center",va="bottom",color=RED,fontsize=10)
    cc.text(0.5,8.95,
            f">8 \u00c5 tail: {tail8.get('frozen',0)}/{nf} frozen, {tail8.get('relaxed',0)}/{nr} relaxed",
            ha="center",va="top",color="0.4",fontsize=10)
    cc.set_xticks([0,1]); cc.set_xticklabels(xt); cc.set_xlim(-0.5,1.6)
    cc.set_ylabel("max adjacent P-P gap (\u00c5)")
    n9=tail9.get('frozen',0)+tail9.get('relaxed',0)
    tc=("No excess fraying vs frozen: none \u22659\u00c5" if n9==0 and tail8.get('relaxed',0)<=tail8.get('frozen',0)+1
        else (f"{n9} state(s) reach the \u22659\u00c5 unwind regime" if n9>0 else "Fraying / helix-unwind check"))
    cc.set_title(tc,loc="center",fontsize=12)
    cc.set_title("c",loc="left",fontweight="bold",fontsize=15)

    # ---- panel d: delta global bend angle ----
    md=np.median(finite("relaxed","dbend")) if nr else np.nan
    xt=[]
    for i,(cond,col,lab) in enumerate(CONDS):
        v=ps[cond].get("dbend",np.array([]))
        jitter(dd,i,v,col,seed=9+i); median_bar(dd,i,v,col,fmt="med {:.2f}",side="right",unit="\u00b0")
        xt.append(lab)
    dd.axhline(0.0,ls="--",color="0.6",lw=1)
    dd.set_xticks([0,1]); dd.set_xticklabels(xt); dd.set_xlim(-0.5,1.6)
    dd.set_ylabel("\u0394 global bend angle vs docked (\u00b0)")
    if np.isfinite(md):
        td=(f"Relaxation increases bend (median {md:+.2f}\u00b0)" if md>0.1
            else (f"Relaxation reduces bend (median {md:+.2f}\u00b0)" if md<-0.1
                  else "Relaxation leaves global bend unchanged"))
    else:
        td="Global bend change vs docked"
    dd.set_title(td,loc="center",fontsize=12)
    dd.set_title("d",loc="left",fontweight="bold",fontsize=15)

    for ax in (a,b,cc,dd):
        ax.spines[["top","right"]].set_visible(False)
    fig.tight_layout(rect=[0,0,1,0.965])
    os.makedirs(os.path.dirname(out) or ".",exist_ok=True)
    fig.savefig(out,dpi=150,bbox_inches="tight")
    print("wrote",out)

if __name__=="__main__":
    main()
