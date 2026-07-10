"""
DNA-relaxation smoke-test analysis.

Given the Stage-2 docked input and three Stage-3 minimized outputs
(A=frozen default, B=dna_k=same-as-protein, C=dna_k=0), quantify how far the
DNA moved and check for structural blow-up (NaN, backbone fraying).

Metrics per output:
  - DNA heavy-atom RMSD vs the docked input (raw; protein backbone is pinned to
    absolute coords in every run, so no rigid-body alignment is needed and raw
    RMSD is directly comparable across runs). Reported both all-DNA-heavy and
    backbone (P, C1').
  - Protein-CA RMSD vs input (sanity: should be small in all runs).
  - NaN check; DNA backbone (P-P adjacent) bond-length fraying check.

Compat assertion: A vs B DNA-RMSD must be ~identical (< 0.05 A) -> the split
refactor is inert when DNA k == protein k.
"""
import sys, numpy as np

AA = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS",
      "MET","PHE","PRO","SER","THR","TRP","TYR","VAL","HID","HIE","HIP"}
DNA = {"DA","DG","DC","DT","DA5","DG5","DC5","DT5","DA3","DG3","DC3","DT3"}

def parse(path):
    """Return dict (chain,resseq,icode,atomname)->xyz, plus resname map."""
    atoms = {}
    with open(path) as fh:
        for ln in fh:
            if not (ln.startswith("ATOM") or ln.startswith("HETATM")):
                continue
            name = ln[12:16].strip()
            resn = ln[17:20].strip()
            chain = ln[21]
            resseq = ln[22:26].strip()
            icode = ln[26]
            x = float(ln[30:38]); y = float(ln[38:46]); z = float(ln[46:54])
            atoms[(chain, resseq, icode, name)] = (resn, np.array([x, y, z]))
    return atoms

def common_rmsd(a, b, sel):
    """RMSD over atoms present in both a and b passing sel(resn, name)."""
    xs_a, xs_b = [], []
    for key, (resn, xyz) in a.items():
        name = key[3]
        if not sel(resn, name):
            continue
        if key in b:
            xs_a.append(xyz); xs_b.append(b[key][1])
    if not xs_a:
        return np.nan, 0
    A = np.array(xs_a); B = np.array(xs_b)
    return float(np.sqrt(((A - B) ** 2).sum(axis=1).mean())), len(xs_a)

def has_nan(a):
    return any(np.isnan(xyz).any() for _, xyz in a.values())

def dna_bb_sel(resn, name):
    return resn in DNA and name in ("P", "C1'")
def dna_heavy_sel(resn, name):
    return resn in DNA and not name.startswith("H")
def ca_sel(resn, name):
    return resn in AA and name == "CA"

def max_pp_gap(a):
    """Max adjacent P-P distance within each DNA chain (fraying proxy).
    B-DNA P-P ~ 6-7 A; a large value flags a broken/frayed backbone."""
    by_chain = {}
    for (chain, resseq, icode, name), (resn, xyz) in a.items():
        if resn in DNA and name == "P":
            by_chain.setdefault(chain, []).append((int(resseq), xyz))
    worst = 0.0
    for chain, lst in by_chain.items():
        lst.sort()
        for (r1, p1), (r2, p2) in zip(lst, lst[1:]):
            if r2 == r1 + 1:
                worst = max(worst, float(np.linalg.norm(p2 - p1)))
    return worst

if __name__ == "__main__":
    inp = sys.argv[1]
    labels = sys.argv[2::2]
    paths = sys.argv[3::2]
    ref = parse(inp)
    print(f"# reference (docked input): {inp}")
    print(f"# ref DNA-heavy P-P max gap: {max_pp_gap(ref):.2f} A")
    print(f"{'label':<20}{'DNA_heavy_RMSD':>16}{'DNA_bb_RMSD':>14}{'CA_RMSD':>10}{'maxPPgap':>10}{'NaN':>6}")
    results = {}
    for lab, p in zip(labels, paths):
        try:
            o = parse(p)
        except FileNotFoundError:
            print(f"{lab:<20}{'MISSING - run failed':>16}")
            continue
        dh, ndh = common_rmsd(ref, o, dna_heavy_sel)
        db, ndb = common_rmsd(ref, o, dna_bb_sel)
        ca, nca = common_rmsd(ref, o, ca_sel)
        gap = max_pp_gap(o)
        nan = has_nan(o)
        results[lab] = (dh, db, ca, gap, nan)
        print(f"{lab:<20}{dh:>16.3f}{db:>14.3f}{ca:>10.3f}{gap:>10.2f}{str(nan):>6}")
    # compat assertion A vs B
    if "A_frozen" in results and "B_dnak_same" in results:
        da = results["A_frozen"][0]; db = results["B_dnak_same"][0]
        diff = abs(da - db)
        print(f"\n# COMPAT A vs B: DNA_heavy_RMSD diff = {diff:.4f} A "
              f"({'PASS <0.05' if diff < 0.05 else 'REVIEW >=0.05'})")
    if "C_dnak0" in results:
        print(f"# RELAX C: DNA_heavy_RMSD = {results['C_dnak0'][0]:.3f} A "
              f"vs frozen {results.get('A_frozen',[float('nan')])[0]:.3f} A")
