#!/usr/bin/env python3
"""
Aggregate the pyCurves ensemble (legacy + curvesplus JSONs) into tidy tables.

Per structure it extracts:
  legacy convention:
    overall_bend_uu, overall_bend_pp, shortening_percent  (global_axis_bending_summary)
    mean_minor_width, mean_major_width                     (groove, over all levels/sublevels)
  curvesplus convention:
    mean_inclin, mean_tip, mean_xdisp, mean_ydisp          (curvesplus_base_pair_axis)
    mean_minor_width, mean_major_width                     (groove)

Outputs:
  <out>_perstructure.csv   one row per (tf, cond, state) with legacy + curvesplus cols merged
  <out>_summary.csv        per (tf, cond) medians + IQRs of the headline metrics

Usage: python aggregate_pycurves.py <pycurves_root> <out_prefix>
"""
import json, glob, os, sys, csv, statistics as st
root, outpref = sys.argv[1], sys.argv[2]

def groove_means(d):
    g = d["dataframes"].get("groove", {}).get("data", {})
    mn, mj = [], []
    for lvl, entry in g.items():
        subs = entry.get("sub_levels", {}) if isinstance(entry, dict) else {}
        for s, vals in subs.items():
            mw, Mw = vals.get("minor_width"), vals.get("major_width")
            if isinstance(mw, (int, float)): mn.append(mw)
            if isinstance(Mw, (int, float)): mj.append(Mw)
    return (st.mean(mn) if mn else None), (st.mean(mj) if mj else None)

def legacy_metrics(path):
    d = json.load(open(path))
    df = d["dataframes"]
    s = df.get("global_axis_bending_summary")
    bend_uu = bend_pp = short = None
    if s:
        r0 = s[0]
        bend_uu = r0.get("overall_bend_uu"); bend_pp = r0.get("overall_bend_pp")
        short = r0.get("shortening_percent")
    mn, mj = groove_means(d)
    return dict(bend_uu=bend_uu, bend_pp=bend_pp, shortening=short,
                minor_w=mn, major_w=mj)

def curvesplus_metrics(path):
    d = json.load(open(path))
    df = d["dataframes"]
    bp = df.get("curvesplus_base_pair_axis", [])
    def m(key):
        v = [r[key] for r in bp if isinstance(r.get(key), (int, float))]
        return st.mean(v) if v else None
    mn, mj = groove_means(d)
    return dict(inclin=m("inclin"), tip=m("tip"), xdisp=m("xdisp"), ydisp=m("ydisp"),
                minor_w_cp=mn, major_w_cp=mj)

# discover all legacy files -> (tf, cond, state); match curvesplus partner
rows = []
for lf in sorted(glob.glob(f"{root}/*/*_legacy.json")):
    tf = os.path.basename(os.path.dirname(lf))
    base = os.path.basename(lf)[:-len("_legacy.json")]   # e.g. relaxed_state_001 or crystal
    cf = lf[:-len("_legacy.json")] + "_curvesplus.json"
    # parse cond/state from base
    if base == "crystal":
        cond, state = "crystal", "crystal"
    else:
        parts = base.split("_state_")
        cond = parts[0]; state = parts[1] if len(parts) > 1 else ""
    row = dict(tf=tf, cond=cond, state=state)
    try: row.update(legacy_metrics(lf))
    except Exception as e: row["legacy_err"] = str(e)[:40]
    if os.path.exists(cf):
        try: row.update(curvesplus_metrics(cf))
        except Exception as e: row["cp_err"] = str(e)[:40]
    rows.append(row)

cols = ["tf","cond","state","bend_uu","bend_pp","shortening","minor_w","major_w",
        "inclin","tip","xdisp","ydisp","minor_w_cp","major_w_cp"]
with open(f"{outpref}_perstructure.csv","w",newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore"); w.writeheader()
    for r in rows:
        w.writerow({k: (round(r[k],4) if isinstance(r.get(k),(int,float)) else r.get(k,"")) for k in cols})

# per (tf, cond) summary
def stats(vals):
    v = [x for x in vals if isinstance(x,(int,float))]
    if not v: return (None,None,None,len(v))
    v.sort()
    med = st.median(v)
    q1 = v[len(v)//4]; q3 = v[(3*len(v))//4]
    return (round(med,3), round(q1,3), round(q3,3), len(v))
from collections import defaultdict
groups = defaultdict(list)
for r in rows: groups[(r["tf"],r["cond"])].append(r)
with open(f"{outpref}_summary.csv","w",newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["tf","cond","n","bend_uu_med","bend_uu_q1","bend_uu_q3",
                "bend_pp_med","minor_w_med","major_w_med","inclin_med","tip_med"])
    for (tf,cond) in sorted(groups):
        g = groups[(tf,cond)]
        buu = stats([r.get("bend_uu") for r in g])
        bpp = stats([r.get("bend_pp") for r in g])
        mnw = stats([r.get("minor_w") for r in g])
        mjw = stats([r.get("major_w") for r in g])
        inc = stats([r.get("inclin") for r in g])
        tip = stats([r.get("tip") for r in g])
        w.writerow([tf,cond,len(g),buu[0],buu[1],buu[2],bpp[0],mnw[0],mjw[0],inc[0],tip[0]])
print(f"wrote {outpref}_perstructure.csv ({len(rows)} rows) + {outpref}_summary.csv")
# quick console summary: bend by cond per TF
print("\nTF        cond      n   bend_uu_med  minor_w_med")
for (tf,cond) in sorted(groups):
    g=groups[(tf,cond)]
    buu=stats([r.get("bend_uu") for r in g]); mnw=stats([r.get("minor_w") for r in g])
    print(f"{tf:9s} {cond:8s} {len(g):3d}   {str(buu[0]):>10s}   {str(mnw[0]):>10s}")
