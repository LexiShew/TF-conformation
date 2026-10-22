#!/usr/bin/env python
"""
survey_source_chains.py - size every structures/source_chains/<pdb>_chains/<pdb>.cif
so the BioEmu2 batch can be ordered small-first and oversized systems flagged
BEFORE 698 GPU jobs are queued.

Why a survey at all
-------------------
sample_complex.py builds its system with System.from_cif(), which imports the
whole deposited assembly. BioEmu2's model card states a ~10k-token ceiling, and
runtime/VRAM both scale with token count, so the deposited size is the single
best predictor of whether a task will run in minutes, run in hours, or OOM.
Surveying first costs ~1 min of login-node CPU and turns an unknown 698-job
queue into an ordered worklist with the risky tail identified.

Token proxy
-----------
n_tokens_proxy = protein residues + nucleotides. BioEmu2 tokenizes per residue
for polymers (plus one per ligand atom group), so this is a lower bound on the
real token count, not an exact figure - it is used for ORDERING and for a
coarse over-ceiling flag, nothing quantitative.

Env: needs gemmi, which on this cluster lives in the pycurves env only.
    /project2/rohs_102/shewchuk/conda/envs/pycurves/bin/python

Usage
-----
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
      /project2/rohs_102/shewchuk/conda/envs/pycurves/bin/python \
      bioemu2/survey_source_chains.py \
      --out bioemu2/worklist/source_chains_survey.tsv

Output: TSV, one row per structure, sorted ascending by n_tokens_proxy:
    pdb_id  cif  n_protein_res  n_dna_nt  n_rna_nt  n_poly_chains  ligands  n_tokens_proxy  status
`status` is ok / no_dna / no_protein / parse_error - the batch driver skips
anything that is not ok, because a protein-only or DNA-less system is not the
complex this track exists to sample.
"""
import argparse
import sys
from pathlib import Path

import gemmi

_HERE = Path(__file__).resolve().parent
_TFCONF = _HERE.parent
SOURCE_CHAINS = _TFCONF / "structures" / "source_chains"

PEPTIDE = {gemmi.PolymerType.PeptideL, gemmi.PolymerType.PeptideD}
DNA = {gemmi.PolymerType.Dna}
RNA = {gemmi.PolymerType.Rna}
HYBRID = {gemmi.PolymerType.DnaRnaHybrid}


def survey_one(cif):
    """Residue/nucleotide counts for one deposited assembly."""
    st = gemmi.read_structure(str(cif))
    st.setup_entities()
    st.remove_waters()
    prot = dna = rna = 0
    npoly = 0
    ligands = set()
    for ch in st[0]:
        poly = ch.get_polymer()
        n = len(poly)
        if n:
            npoly += 1
            pt = poly.check_polymer_type()
            if pt in PEPTIDE:
                prot += n
            elif pt in DNA or pt in HYBRID:
                dna += n
            elif pt in RNA:
                rna += n
        # non-polymer residues in this chain = heterogens / metals
        for res in ch:
            if poly and res in poly:
                continue
            if res.het_flag == "H":
                ligands.add(res.name)
    return prot, dna, rna, npoly, sorted(ligands)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", default=str(_HERE / "worklist" / "source_chains_survey.tsv"))
    a = p.parse_args()

    rows = []
    dirs = sorted(SOURCE_CHAINS.glob("*_chains"))
    for i, d in enumerate(dirs, 1):
        pdb_id = d.name.replace("_chains", "")
        cif = d / f"{pdb_id}.cif"
        if not cif.exists():
            rows.append((pdb_id, str(cif), 0, 0, 0, 0, "", 0, "missing_cif"))
            continue
        try:
            prot, dna, rna, npoly, ligs = survey_one(cif)
        except Exception as e:                                    # noqa: BLE001
            rows.append((pdb_id, str(cif), 0, 0, 0, 0, "", 0,
                         f"parse_error:{type(e).__name__}"))
            continue
        status = "ok"
        if prot == 0:
            status = "no_protein"
        elif dna == 0:
            status = "no_dna"
        rows.append((pdb_id, str(cif), prot, dna, rna, npoly,
                     ";".join(ligs), prot + dna + rna, status))
        if i % 100 == 0:
            print(f"  ... {i}/{len(dirs)}", file=sys.stderr, flush=True)

    rows.sort(key=lambda r: (r[8] != "ok", r[7]))
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        fh.write("pdb_id\tcif\tn_protein_res\tn_dna_nt\tn_rna_nt\tn_poly_chains"
                 "\tligands\tn_tokens_proxy\tstatus\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    ok = [r for r in rows if r[8] == "ok"]
    print(f"surveyed {len(rows)} structures -> {out}")
    print(f"  ok            : {len(ok)}")
    for st in sorted({r[8] for r in rows if r[8] != "ok"}):
        print(f"  {st:<14}: {sum(1 for r in rows if r[8] == st)}")
    if ok:
        tok = sorted(r[7] for r in ok)
        import statistics
        print(f"  tokens_proxy  : min={tok[0]} median={int(statistics.median(tok))} "
              f"max={tok[-1]}")
        for ceil in (2000, 4000, 6000, 10000):
            print(f"    > {ceil:<6}: {sum(1 for t in tok if t > ceil)}")


if __name__ == "__main__":
    main()
