
import sys,re,statistics as st
def parse(path):
    L=open(path).read().splitlines()
    start=None
    for i,l in enumerate(L):
        if "Minor groove" in l and "Major groove" in l: start=i+2; break
    minor=[]; major=[]; seen=False
    for l in L[start:]:
        if not l.strip():
            if seen: break     # blank AFTER data = end of table
            else: continue     # leading blank(s) = skip
        # stop if we hit next section header (starts with letter block like '|L|')
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
for tag,f in [("crystal","tbp_crystal.txt"),("frozen","tbp_state_002_frozen.txt"),("relax","tbp_state_002_relax.txt")]:
    mn,mj=parse(f"{sys.argv[1]}/{f}")
    def s(x): return f"n={len(x)} mean={st.mean(x):.2f} med={st.median(x):.2f} range=[{min(x):.1f},{max(x):.1f}]" if x else "none"
    print(f"{tag}: MINOR {s(mn)} | MAJOR {s(mj)}")
