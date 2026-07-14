#!/usr/bin/env python3
"""
interface_rmsd.py — label-robust protein-DNA binding-fidelity metrics.

Single pair OR batch over a monomer_chains tree. Scores EVERY model in a
multi-MODEL docked PDB (one MODEL per BioEmu conformation).

Metrics per state:
  iRMSD_global  : superpose on all DNA-contacting protein residues (backbone), RMSD
  iRMSD_seg_max : superpose each contiguous interface segment alone; worst unit
  iRMSD_seg_mean: mean over segments
  fnat          : native protein-DNA contacts recovered, vs REFERENCE DNA (fixed target)

Protein matched ref<->model by sequence alignment (numbering-agnostic).

Single:  python interface_rmsd.py --ref X.cif --model X_docked.pdb --tag X
Batch :  python interface_rmsd.py --batch /path/to/monomer_chains --out metrics.csv
Requires: biopython>=1.79, numpy
"""
import argparse, sys, numpy as np
from pathlib import Path
from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.Polypeptide import is_aa
from Bio.SVDSuperimposer import SVDSuperimposer
from Bio.Align import PairwiseAligner

DNA_RES = {"DA","DC","DG","DT","DI","DU","A","C","G","T","U"}
BACKBONE = ["N","CA","C","O"]
AA3TO1 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
 'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
 'THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
COLS = "pdb_id,state,n_iface_res,n_segments,seq_ident,iRMSD_global,iRMSD_seg_max,iRMSD_seg_mean,fnat"

def load_models(path):
    p = MMCIFParser(QUIET=True) if str(path).lower().endswith(".cif") else PDBParser(QUIET=True)
    return list(p.get_structure("s", str(path)))

def heavy(res): return [a for a in res if a.element != "H"]

def ordered(model):
    prot, dna = [], []
    for ch in model:
        for r in sorted((r for r in ch), key=lambda r: r.id[1]):
            rn = r.get_resname().strip()
            if rn in DNA_RES: dna.append(r)
            elif is_aa(r, standard=True) and rn in AA3TO1: prot.append(r)
    seq = "".join(AA3TO1[r.get_resname().strip()] for r in prot)
    return prot, seq, dna

def min_dist(a, b):
    A = np.array([x.coord for x in heavy(a)]); B = np.array([x.coord for x in heavy(b)])
    if not len(A) or not len(B): return np.inf
    return np.sqrt(((A[:,None,:]-B[None,:,:])**2).sum(-1)).min()

def align_map(seq_ref, seq_mod):
    al = PairwiseAligner(); al.mode="global"
    al.match_score=2; al.mismatch_score=-1; al.open_gap_score=-5; al.extend_gap_score=-0.5
    aln = al.align(seq_ref, seq_mod)[0]
    m = {}
    for (r0,r1),(m0,m1) in zip(*aln.aligned):
        for o in range(r1-r0): m[r0+o] = m0+o
    ident = sum(seq_ref[r]==seq_mod[m[r]] for r in m)/max(1,len(m))
    return m, ident

def backbone_pairs(rr, mr):
    R, M = [], []
    rd = {a.get_name(): a.coord for a in rr}; md = {a.get_name(): a.coord for a in mr}
    for an in BACKBONE:
        if an in rd and an in md: R.append(rd[an]); M.append(md[an])
    return R, M

def rmsd(R, M):
    if len(R) < 3: return float("nan")
    s = SVDSuperimposer(); s.set(np.array(R), np.array(M)); s.run()
    return float(s.get_rms())

def _select_bound_monomer(rp, rseq, rd, iface_cut, contact_cut):
    """When the reference holds >1 protein chain (a second biological assembly
    or a crystal-packing partner), keep only the best-bound monomer: the protein
    chain making the most DNA contacts, plus the DNA chains that chain contacts.
    Single-protein-chain references return unchanged (byte-identical path)."""
    prot_ids = []
    for r in rp:
        cid = r.get_parent().id
        if cid not in prot_ids: prot_ids.append(cid)
    if len(prot_ids) <= 1:
        return rp, rseq, rd
    best_id, best_n, best_dna = prot_ids[0], -1, set()
    for pid in prot_ids:
        chain_res = [r for r in rp if r.get_parent().id == pid]
        n = 0; dna_hit = set()
        for pr in chain_res:
            for nr in rd:
                if min_dist(pr, nr) <= iface_cut:
                    n += 1; dna_hit.add(nr.get_parent().id)
        if n > best_n:
            best_id, best_n, best_dna = pid, n, dna_hit
    rp2 = [r for r in rp if r.get_parent().id == best_id]
    rd2 = [r for r in rd if r.get_parent().id in best_dna]
    seq2 = "".join(AA3TO1[r.get_resname().strip()] for r in rp2)
    sys.stderr.write(f"[fnat] multi-assembly reference: kept protein chain {best_id} "
                     f"+ DNA {sorted(best_dna)} ({len(rp2)} prot res, {len(rd2)} DNA res) "
                     f"of {len(prot_ids)} protein chains\n")
    return rp2, seq2, rd2


def _chain_ids_by_index(model, indices):
    """Map positional chain indices (0-based file order, identical to mdtraj
    chainid and to Stage-2's PROTEIN_CHAIN/DNA_CHAINS) to Bio.PDB chain ids."""
    chains = list(model.get_chains())
    out = []
    for i in indices:
        if i < 0 or i >= len(chains):
            raise ValueError(f"chain index {i} out of range (model has {len(chains)} chains)")
        out.append(chains[i].id)
    return out


def _select_by_config(rp, rseq, rd, prot_id, dna_ids):
    """Select the reference protein + DNA by EXPLICIT chain ids resolved from the
    pilot config's PROTEIN_CHAIN / DNA_CHAINS. Use this instead of the
    contact-count auto-pick so the gate scores the ensemble against the SAME
    protein copy and duplex that Stage-2 docked onto. Essential for multi-copy
    crystals (e.g. 1hjc: two Runt copies each on its own duplex) where the
    auto-pick can choose the other copy and report a spurious fnat=0."""
    rp2 = [r for r in rp if r.get_parent().id == prot_id]
    rd2 = [r for r in rd if r.get_parent().id in dna_ids]
    if not rp2:
        raise ValueError(f"config protein chain id {prot_id!r} matched no reference residues")
    if not rd2:
        raise ValueError(f"config DNA chain ids {sorted(dna_ids)} matched no reference residues")
    seq2 = "".join(AA3TO1[r.get_resname().strip()] for r in rp2)
    sys.stderr.write(f"[fnat] config-selected reference: protein chain {prot_id} "
                     f"+ DNA {sorted(dna_ids)} ({len(rp2)} prot res, {len(rd2)} DNA res)\n")
    return rp2, seq2, rd2


def ref_side(model, iface_cut, contact_cut, gap, protein_chain=None, dna_chains=None):
    rp, rseq, rd = ordered(model)
    if protein_chain is not None and dna_chains is not None:
        # Config-driven: honor the redock's chain selection (positional indices).
        prot_id = _chain_ids_by_index(model, [protein_chain])[0]
        dna_ids = set(_chain_ids_by_index(model, dna_chains))
        rp, rseq, rd = _select_by_config(rp, rseq, rd, prot_id, dna_ids)
    else:
        # Legacy auto-pick: best-bound monomer by DNA-contact count.
        rp, rseq, rd = _select_bound_monomer(rp, rseq, rd, iface_cut, contact_cut)
    iface = [i for i,pr in enumerate(rp) if any(min_dist(pr,nr)<=iface_cut for nr in rd)]
    segs, cur = [], []
    for i in iface:
        if cur and (rp[i].get_parent().id != rp[cur[-1]].get_parent().id
                    or rp[i].id[1]-rp[cur[-1]].id[1] > gap):
            segs.append(cur); cur=[]
        cur.append(i)
    if cur: segs.append(cur)
    native = [(i,j) for i,pr in enumerate(rp) for j,nr in enumerate(rd)
              if min_dist(pr,nr)<=contact_cut]
    return dict(rp=rp, rseq=rseq, rd=rd, iface=iface, segs=segs, native=native)

def score(R, mod_model, contact_cut, use_model_dna=False):
    mp, mseq, md = ordered(mod_model)
    ref2mod, ident = align_map(R['rseq'], mseq); rp = R['rp']
    def over(idxs):
        A,B=[],[]
        for i in idxs:
            if i in ref2mod:
                a,b = backbone_pairs(rp[i], mp[ref2mod[i]]); A+=a; B+=b
        return rmsd(A,B)
    g = over(R['iface'])
    sv = [v for v in (over(s) for s in R['segs']) if v==v]
    smax = max(sv) if sv else float("nan"); smean = float(np.mean(sv)) if sv else float("nan")
    nat = R['native']
    # fnat target: reference DNA (default) or model's own DNA (frame cross-check)
    dna = md if (use_model_dna and len(md) == len(R['rd'])) else R['rd']
    rec = sum(1 for i,j in nat
              if i in ref2mod and j < len(dna) and min_dist(mp[ref2mod[i]], dna[j]) <= contact_cut)
    fnat = rec/len(nat) if nat else float("nan")
    return dict(n_iface=len(R['iface']), n_seg=len(R['segs']), ident=ident,
                g=g, smax=smax, smean=smean, fnat=fnat)

def fmt(x):
    try: x=float(x)
    except (TypeError,ValueError): return "NA"
    return "NA" if np.isnan(x) else "%.3f"%x

def row(pdb, state, m):
    return ",".join([pdb, str(state), str(m['n_iface']), str(m['n_seg']),
        fmt(m['ident']), fmt(m['g']), fmt(m['smax']), fmt(m['smean']), fmt(m['fnat'])])

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref"); ap.add_argument("--model"); ap.add_argument("--tag", default="")
    ap.add_argument("--batch"); ap.add_argument("--out", default="interface_metrics.csv")
    ap.add_argument("--fail", default="interface_fails.log")
    ap.add_argument("--iface_cutoff", type=float, default=5.0)
    ap.add_argument("--contact_cutoff", type=float, default=4.5)
    ap.add_argument("--gap", type=int, default=4)
    ap.add_argument("--min_ident", type=float, default=0.9)
    ap.add_argument("--use_model_dna", action="store_true",
                    help="measure fnat vs the model's own DNA (coordinate-frame cross-check)")
    # NB: per-directory state scoring + the single fnat quality gate live in
    # score_stage3.py / fnat_gate/fnat_gate.sh (Stage 3 -> Stage 4, B7). This
    # module is the scorer only — it never moves or drops state files.
    a = ap.parse_args()

    if a.batch:
        root = Path(a.batch)
        with open(a.out,"w") as out, open(a.fail,"w") as flog:
            out.write(COLS+"\n")
            for cd in sorted(root.glob("*_chains")):
                pdb = cd.name[:-7]  # strip "_chains"
                ref = cd/f"{pdb}.cif"; conf = cd/f"{pdb}_conformations"
                if not ref.exists() or not conf.exists():
                    flog.write(f"{pdb}\tMISSING_REF_OR_CONF\n"); continue
                try: R = ref_side(load_models(ref)[0], a.iface_cutoff, a.contact_cutoff, a.gap)
                except Exception as e: flog.write(f"{pdb}\tREF_ERROR\t{e}\n"); continue
                docked = sorted(conf.glob("*docked*.pdb"))
                if not docked: flog.write(f"{pdb}\tNO_DOCKED\n"); continue
                n = 0
                for df in docked:
                    try: models = load_models(df)
                    except Exception as e: flog.write(f"{pdb}\t{df.name}\tLOAD_ERROR\t{e}\n"); continue
                    for k, mm in enumerate(models, 1):
                        try: m = score(R, mm, a.contact_cutoff, a.use_model_dna)
                        except Exception as e: flog.write(f"{pdb}\tstate{k}\tSCORE_ERROR\t{e}\n"); continue
                        if m['ident'] < a.min_ident:
                            flog.write(f"{pdb}\tstate{k}\tLOW_IDENT\t{m['ident']:.2f}\n")
                        out.write(row(pdb, k, m)+"\n"); n += 1
                print(f"{pdb}: {n} states", file=sys.stderr)
        print(f"done -> {a.out} (warnings/failures in {a.fail})", file=sys.stderr)
        return

    R = ref_side(load_models(a.ref)[0], a.iface_cutoff, a.contact_cutoff, a.gap)
    tag = (a.tag.replace(",","_") or Path(a.model).stem)
    print(COLS)
    for k, mm in enumerate(load_models(a.model), 1):
        print(row(tag, k, score(R, mm, a.contact_cutoff, a.use_model_dna)))

if __name__ == "__main__":
    main()