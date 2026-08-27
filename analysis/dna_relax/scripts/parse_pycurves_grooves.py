#!/usr/bin/env python3
"""
Parse minor/major groove widths and overall bend from pyCurves .txt output.
pyCurves must be run first (see run_pycurves.sh). Usage:
  python parse_pycurves_grooves.py <pycurves_txt_dir>
Expects files *_crystal.txt / *_frozen.txt / *_relax.txt in the dir.
Reports per-condition minor/major groove width stats + overall axis bend (UU/PP).
"""
import sys,re,glob,os,statistics as st
def parse_grooves(path):
    L=open(path).read().splitlines()
    start=None
    for i,l in enumerate(L):
        if "Minor groove" in l and "Major groove" in l: start=i+2; break
    minor=[]; major=[]; seen=False
    if start is None: return minor,major
    for l in L[start:]:
        if not l.strip():
            if seen: break
            else: continue
        if re.match(r'\s*\|[A-Z]\|', l) or l.strip().startswith("---"): break
        seen=True
        nums=[t for t in l.replace("*"," ").split() if re.match(r'^-?\d+\.\d+$',t)]
        if len(nums)>=1:
            mw=float(nums[0])
            if 3<mw<20: minor.append(mw)
        if len(nums)>=4:
            Mw=float(nums[3])
            if 8<Mw<30: major.append(Mw)
    return minor,major
def parse_bend(path):
    for l in open(path):
        if "Overall axis bend" in l:
            m=re.findall(r'UU=\s*([\d.]+).*PP=\s*([\d.]+)', l)
            if m: return float(m[0][0]), float(m[0][1])
    return None,None
def s(x): return f"n={len(x)} mean={st.mean(x):.2f} med={st.median(x):.2f} range=[{min(x):.1f},{max(x):.1f}]" if x else "none"
d=sys.argv[1]
for f in sorted(glob.glob(f"{d}/*.txt")):
    mn,mj=parse_grooves(f); uu,pp=parse_bend(f)
    print(f"{os.path.basename(f)}: bend UU={uu} PP={pp} | MINOR {s(mn)} | MAJOR {s(mj)}")
