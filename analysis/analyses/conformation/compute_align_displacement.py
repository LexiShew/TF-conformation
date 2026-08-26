#!/usr/bin/env python3
"""
compute_align_displacement.py — quantify how differently interface-Cα vs global
all-Cα alignment places the carried DNA, per Stage-1 conformer, per pilot.

For each pilot it re-docks every Stage-1 frame twice with stage2_redock.py
(--align-mode interface and --align-mode all), then measures the per-state DNA-P
RMSD between the two docked outputs. Because both modes keep the BioEmu protein
in the same coordinate frame (only the carried reference DNA is transformed by
the mode-specific Kabsch fit), this displacement IS the alignment-induced DNA
misplacement — no crystal superposition, so it is immune to chain-mapping
ambiguity. Writes analysis/analyses/align_compare/dna_displacement_interface_vs_global.csv
(cols: pilot, state, dna_disp), consumed by make_aligncompare.py.

Docking is cheap (~seconds/pilot, CPU); no minimization is needed for this
metric. Env: bioemu (mdtraj + stage2_redock deps) for docking, biopython for
the P-atom read.

Usage:
  conda activate bioemu
  python compute_align_displacement.py --tfs ets1 tbp lef1 [--nstates 100]

Pilot chain config (protein_chain / dna_chains positional indices) is read from
config/pilots/<tf>.sh, matching what the pipeline uses.
"""
import argparse
import glob
import os
import re
import subprocess
import sys

import numpy as np

BASE = os.environ.get("TFCONF_DIR", "/project2/rohs_102/shewchuk/TF-conformation")
DNA_RES = {"DA", "DC", "DG", "DT", "DI", "DU", "A", "C", "G", "T", "U"}


def pilot_cfg(tf):
    """Read PDB_ID / BINDING_CHAIN / PROTEIN_CHAIN / DNA_CHAINS from the pilot config."""
    txt = open(f"{BASE}/config/pilots/{tf}.sh").read()
    def grab(key):
        m = re.search(rf'{key}="?([^"\n]+)"?', txt)
        return m.group(1).strip() if m else None
    return dict(pdb=grab("PDB_ID"), bc=grab("BINDING_CHAIN"),
                pc=grab("PROTEIN_CHAIN"), dc=grab("DNA_CHAINS"))


def dock(tf, cfg, mode, outdir):
    pdb, bc, pc, dc = cfg["pdb"], cfg["bc"], cfg["pc"], cfg["dc"]
    s1 = f"{BASE}/output/stage1_bioemu/{pdb}_chain{bc}_conformations"
    ref = f"{BASE}/structures/source_chains/{pdb}_chains/{pdb}.cif"
    os.makedirs(outdir, exist_ok=True)
    subprocess.run([
        sys.executable, f"{BASE}/stage2_redock/stage2_redock.py",
        "--pdb-id", pdb, "--ref", ref,
        "--traj", f"{s1}/samples_sidechain_rec.xtc",
        "--top", f"{s1}/samples_sidechain_rec.pdb",
        "--out-dir", outdir, "--protein-chain", str(pc),
        "--dna-chains", dc, "--align-mode", mode,
    ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def dna_P(path):
    from Bio.PDB import PDBParser
    model = list(PDBParser(QUIET=True).get_structure("s", path))[0]
    xs = [a.coord for ch in model for r in ch
          if r.resname.strip() in DNA_RES for a in r if a.name == "P"]
    return np.array(xs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tfs", nargs="+", default=["ets1", "tbp", "lef1"])
    ap.add_argument("--nstates", type=int, default=100)
    ap.add_argument("--out", default=f"{BASE}/analysis/analyses/align_compare")
    a = ap.parse_args()
    work = os.path.join(a.out, "_work")
    rows = []
    for tf in a.tfs:
        cfg = pilot_cfg(tf)
        pdb = cfg["pdb"]
        di_dir = os.path.join(work, tf, "interface_docked")
        da_dir = os.path.join(work, tf, "all_docked")
        dock(tf, cfg, "interface", di_dir)
        dock(tf, cfg, "all", da_dir)
        n = 0
        for fi in sorted(glob.glob(f"{di_dir}/{pdb}_state_*.pdb"))[:a.nstates]:
            fa = os.path.join(da_dir, os.path.basename(fi))
            if not os.path.exists(fa):
                continue
            pi, pa = dna_P(fi), dna_P(fa)
            m = min(len(pi), len(pa))
            if m < 3:
                continue
            disp = float(np.sqrt(((pi[:m] - pa[:m]) ** 2).sum(1).mean()))
            rows.append((tf, os.path.basename(fi), round(disp, 3)))
            n += 1
        med = np.median([r[2] for r in rows if r[0] == tf]) if n else float("nan")
        print(f"{tf}: {n} states, median displacement {med:.2f} Å")
    out_csv = os.path.join(a.out, "dna_displacement_interface_vs_global.csv")
    os.makedirs(a.out, exist_ok=True)
    with open(out_csv, "w") as fh:
        fh.write("pilot,state,dna_disp\n")
        for tf, st, d in rows:
            fh.write(f"{tf},{st},{d}\n")
    print("wrote", out_csv, len(rows), "rows")


if __name__ == "__main__":
    main()
