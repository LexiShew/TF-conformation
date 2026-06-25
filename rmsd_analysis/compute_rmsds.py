#!/usr/bin/env python
"""
Compute backbone, sidechain, and total heavy-atom RMSDs for every conformation
at every applicable stage of the DeepPBS BioEmu augmentation pipeline, for
each augmentation variant.

For each (TF, state, stage, variant) we compute two flavors of RMSD:

  1. vs_reference: alignment + RMSD against the reference crystal protein chain
  2. delta_prev:   alignment + RMSD against the same state at the previous stage
                   (e.g. stage 2 delta = state at stage 2 vs same state at stage 1)

Variant handling: stages 0–2 are shared across augmentation variants (BioEmu →
HPacker → redock), so those rows are emitted once with variant="shared". Stage
3 (minimization) diverges per variant — "metal_cage" reads stage3_min/,
"legacy" reads stage3_min_legacy/ — and produces variant-keyed rows.

Atom subsets per pair:

  - backbone:  N, CA, C, O atoms
  - sidechain: heavy atoms NOT in {N, CA, C, O} (residue-name filter excludes DNA)
  - total:     all protein heavy atoms
  - interface: heavy atoms of residues listed in the per-pilot interface JSON
               (built by identify_interface_residues.py). Only populated when
               --interface-dir is provided and a matching {tf}.json is found.

Stage 0 (BioEmu raw) is backbone-only, so sidechain/total RMSDs for stage 0
are reported as NaN. Stage 2 carries reference DNA across by Kabsch transform
but the BioEmu protein backbone is unchanged from Stage 1, so we expect Stage
2 backbone-vs-Stage-1 delta to be essentially zero modulo rotation+translation.

Outputs:
  - per_state_rmsds.csv: long-format (tf, pdb_id, state, stage, variant,
                        comparison, rmsd_backbone_A, rmsd_sidechain_A,
                        rmsd_total_A, n_backbone_atoms, n_sidechain_atoms,
                        n_total_atoms)
  - summary_stats.csv:   per-TF / per-stage / per-variant / per-comparison
                        median/IQR/max
  - Both written to --output-dir (default: cwd)

Usage:
    python compute_rmsds.py --tfs egr1 dux4 tbp --output-dir ./rmsd_analysis
    python compute_rmsds.py --variants metal_cage          # only metal_cage
"""
import argparse
import glob
import os
import sys
import warnings
from pathlib import Path
import numpy as np

# mdtraj imports late so the script can show --help without env activation
def import_mdtraj():
    global md
    import mdtraj as md_mod
    md = md_mod


# ---------------------------------------------------------------------------
# Paths and constants
# ---------------------------------------------------------------------------
# Defaults: cluster layout. Override with CLI flags or env vars when running
# locally (see --bioemu-root / --conf-root / --pilots-dir).
# Defaults follow the current repo layout (self-located); override via env/flags.
#   BIOEMU_ROOT — source chains:   <root>/<pdb>_chains/<pdb>.cif
#   STAGE1_OUT  — Stage 1 library: <root>/<pdb>_chain<X>_conformations/
#                 (topology.pdb, samples.xtc, samples_sidechain_rec.{pdb,xtc})
#   CONF_ROOT   — pipeline output: <root>/<stage>/<tf>/<pdb>_state_NNN.pdb
#   PILOTS_DIR  — pilot configs:   <root>/<tf>.sh
_TFCONF = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIOEMU_ROOT = os.environ.get("DEEPPBS_BIOEMU_ROOT", f"{_TFCONF}/structures/source_chains")
STAGE1_OUT  = os.environ.get("DEEPPBS_STAGE1_OUT",  f"{_TFCONF}/output/stage1_bioemu")
CONF_ROOT   = os.environ.get("DEEPPBS_CONF_ROOT",   f"{_TFCONF}/output")
PILOTS_DIR  = os.environ.get("DEEPPBS_PILOTS_DIR",  f"{_TFCONF}/config/pilots")

BACKBONE_ATOMS = {"N", "CA", "C", "O"}
PROTEIN_RESNAMES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "HID", "HIE", "HIP", "HSD", "HSE", "HSP",
}

# Augmentation variants → stage-3 directory name under <tf>/.
# Stages 0–2 are shared between variants.
VARIANTS = {
    "metal_cage": "stage3_min",
    "legacy":     "stage3_min_legacy",
}


# ---------------------------------------------------------------------------
# Pilot config (parsed from bash files)
# ---------------------------------------------------------------------------
def load_pilot_config(tf_name):
    """Parse the bash pilot config to extract PDB_ID, PROTEIN_CHAIN, N_FRAMES."""
    path = os.path.join(PILOTS_DIR, f"{tf_name}.sh")
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Pilot config not found: {path}")
    cfg = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("export ") and "=" in line:
                key_val = line[len("export "):].split("=", 1)
                key = key_val[0].strip()
                val = key_val[1].strip().strip('"').split('#')[0].strip().strip('"')
                cfg[key] = val
    return cfg


# ---------------------------------------------------------------------------
# Loaders for each stage
# ---------------------------------------------------------------------------
def load_reference_full(pdb_id):
    """Load the full reference structure (we'll restrict to the protein chain
    at RMSD-compute time via the chain_b argument)."""
    cif_path = f"{BIOEMU_ROOT}/{pdb_id}_chains/{pdb_id}.cif"
    if not os.path.isfile(cif_path):
        raise FileNotFoundError(f"Reference CIF not found: {cif_path}")
    return md.load(cif_path)


def stage1_conf_dir(pdb_id):
    """The Stage 1 per-chain conformations dir for this PDB. stage1_bioemu writes
    <pdb>_chain<X>_conformations/ (one per sampled chain); take the first match
    (the binding chain, if only it was sampled)."""
    hits = sorted(glob.glob(f"{STAGE1_OUT}/{pdb_id}_chain*_conformations"))
    return hits[0] if hits else None


def load_stage0(pdb_id):
    """Load BioEmu raw samples (backbone-only XTC + topology)."""
    s1d = stage1_conf_dir(pdb_id)
    if s1d is None:
        return None
    topo_path = f"{s1d}/topology.pdb"
    xtc_path  = f"{s1d}/samples.xtc"
    if not (os.path.isfile(topo_path) and os.path.isfile(xtc_path)):
        return None
    traj = md.load(xtc_path, top=topo_path)
    # Strip Hs (BioEmu shouldn't have them but be safe)
    sel = traj.topology.select("not element H")
    return traj.atom_slice(sel)


def load_stage1(pdb_id):
    """Load HPACKER-relaxed samples (sidechain-restored)."""
    s1d = stage1_conf_dir(pdb_id)
    if s1d is None:
        return None
    topo_path = f"{s1d}/samples_sidechain_rec.pdb"
    xtc_path  = f"{s1d}/samples_sidechain_rec.xtc"
    if not (os.path.isfile(topo_path) and os.path.isfile(xtc_path)):
        return None
    traj = md.load(xtc_path, top=topo_path)
    sel = traj.topology.select("not element H")
    return traj.atom_slice(sel)


def load_stage_per_state(pdb_id, tf_name, stage_dir_name, n_states):
    """Load per-state PDBs (Stage 2 or 3) as a list, indexed by state number."""
    stage_dir = f"{CONF_ROOT}/{stage_dir_name}/{tf_name}"
    out = [None] * (n_states + 1)  # 1-indexed
    if not os.path.isdir(stage_dir):
        return out
    for i in range(1, n_states + 1):
        path = f"{stage_dir}/{pdb_id}_state_{i:03d}.pdb"
        if os.path.isfile(path):
            try:
                t = md.load(path)
                # Select protein heavy atoms only (DNA and Zn dropped here)
                sel_list = []
                for atom in t.topology.atoms:
                    if atom.element is None or atom.element.symbol == "H":
                        continue
                    if atom.residue.name not in PROTEIN_RESNAMES:
                        continue
                    sel_list.append(atom.index)
                if len(sel_list) > 0:
                    out[i] = t.atom_slice(np.array(sel_list, dtype=int))
            except Exception as e:
                print(f"  WARN: failed to load {path}: {e}", file=sys.stderr)
    return out


def infer_tf_from_pdb(pdb_id):
    """Map a PDB id to its pilot/TF name via the pilot configs' PDB_ID field."""
    for cfg_path in sorted(glob.glob(f"{PILOTS_DIR}/*.sh")):
        tf = os.path.basename(cfg_path)[:-len(".sh")]
        try:
            if load_pilot_config(tf).get("PDB_ID", "").lower() == pdb_id.lower():
                return tf
        except Exception:
            continue
    return pdb_id  # fallback


# ---------------------------------------------------------------------------
# Atom-set matching
# ---------------------------------------------------------------------------
def protein_heavy_atom_order(traj, protein_chain=None):
    """Walk the trajectory and return ordered atom indices for protein heavy
    atoms only, in a canonical order: by residue order (skipping non-protein
    residues), and within each residue, by atom name (alphabetic).

    If `protein_chain` is set (int), restrict to that chain index.

    Returns:
        indices (np.ndarray): atom indices into traj
        per_atom_info (list): list of (residue_position, resname, atom_name)
            with residue_position being the 0-based ordinal of this residue
            among protein residues (so position 0 = first protein residue
            regardless of resSeq).
    """
    indices = []
    info = []
    res_pos = 0
    for chain in traj.topology.chains:
        if protein_chain is not None and chain.index != protein_chain:
            continue
        for res in chain.residues:
            if res.name not in PROTEIN_RESNAMES:
                continue
            # Sort atoms within this residue by name for deterministic ordering
            atoms_sorted = sorted(
                [a for a in res.atoms
                 if a.element is not None and a.element.symbol != "H"],
                key=lambda a: a.name,
            )
            for atom in atoms_sorted:
                indices.append(atom.index)
                info.append((res_pos, res.name, atom.name))
            res_pos += 1
    return np.array(indices, dtype=int), info


def common_atom_indices(traj_a, traj_b, chain_a=None, chain_b=None):
    """Find protein heavy atoms present in both trajectories. Matches by
    (residue_position, residue_name, atom_name) where residue_position is
    the 0-based ordinal among protein residues — NOT resSeq, which differs
    between BioEmu (0-N) and crystal (real PDB numbering).

    Returns (idx_a, idx_b, info) such that traj_a.xyz[:, idx_a] aligns
    atom-for-atom with traj_b.xyz[:, idx_b]; `info` is the per-atom
    (residue_position, resname, atom_name) tuples for the surviving pairs,
    used by callers that want to filter to e.g. interface residues.
    """
    idx_a_full, info_a = protein_heavy_atom_order(traj_a, protein_chain=chain_a)
    idx_b_full, info_b = protein_heavy_atom_order(traj_b, protein_chain=chain_b)

    map_b = {}
    for i, key in enumerate(info_b):
        if key not in map_b:
            map_b[key] = i

    idx_a_out, idx_b_out, info_out = [], [], []
    for i, key in enumerate(info_a):
        if key in map_b:
            idx_a_out.append(idx_a_full[i])
            idx_b_out.append(idx_b_full[map_b[key]])
            info_out.append(key)
    return (np.array(idx_a_out, dtype=int),
            np.array(idx_b_out, dtype=int),
            info_out)


def subset_indices(traj, indices, atom_set):
    """Filter `indices` to just those whose atom name is in atom_set (or NOT in it)."""
    if atom_set == "backbone":
        keep = [i for i in indices if traj.topology.atom(i).name in BACKBONE_ATOMS]
    elif atom_set == "sidechain":
        keep = [i for i in indices if traj.topology.atom(i).name not in BACKBONE_ATOMS]
    elif atom_set == "total":
        keep = list(indices)
    else:
        raise ValueError(atom_set)
    return np.array(keep, dtype=int)


# ---------------------------------------------------------------------------
# Kabsch RMSD
# ---------------------------------------------------------------------------
def kabsch_rmsd(coords_a, coords_b):
    """Aligned RMSD between two sets of corresponding 3D coordinates (in nm).
    Returns RMSD in Ångström. Both arrays must be (N, 3) and same length."""
    if coords_a.shape != coords_b.shape:
        return float("nan")
    if len(coords_a) == 0:
        return float("nan")
    # Center
    c_a = coords_a - coords_a.mean(axis=0)
    c_b = coords_b - coords_b.mean(axis=0)
    # Kabsch
    H = c_a.T @ c_b
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    a_aligned = c_a @ R.T
    diff = a_aligned - c_b
    rmsd_nm = np.sqrt((diff ** 2).sum() / len(coords_a))
    return rmsd_nm * 10.0  # nm → Å


_NAN_RESULT = (float("nan"), float("nan"), float("nan"), float("nan"),
               0, 0, 0, 0)


def compute_rmsd_triple(traj_a, frame_a, traj_b, frame_b, chain_a=None,
                        chain_b=None, interface_positions=None):
    """For one frame pair, return (backbone, sidechain, total, interface)
    RMSDs in Å, plus their atom counts. `interface_positions` (set of int
    res-positions) selects which residues count toward the interface RMSD;
    NaN/0 if not provided or empty.

    Uses common-atom intersection by (residue_position, resname, atom_name)
    after restricting each trajectory to its protein chain (if specified).
    """
    if traj_a is None or traj_b is None:
        return _NAN_RESULT
    if frame_a >= traj_a.n_frames or frame_b >= traj_b.n_frames:
        return _NAN_RESULT

    idx_a, idx_b, info = common_atom_indices(traj_a, traj_b,
                                              chain_a=chain_a, chain_b=chain_b)
    if len(idx_a) == 0:
        return _NAN_RESULT

    coords_a_full = traj_a.xyz[frame_a, idx_a, :]
    coords_b_full = traj_b.xyz[frame_b, idx_b, :]

    # Backbone subset (atoms whose name is in BACKBONE_ATOMS in traj_a)
    bb_mask = np.array([
        traj_a.topology.atom(i).name in BACKBONE_ATOMS for i in idx_a
    ])
    if bb_mask.any():
        rmsd_bb = kabsch_rmsd(coords_a_full[bb_mask], coords_b_full[bb_mask])
        n_bb = int(bb_mask.sum())
    else:
        rmsd_bb = float("nan")
        n_bb = 0

    # Sidechain subset (heavy atoms not in BACKBONE_ATOMS)
    sc_mask = ~bb_mask
    if sc_mask.any():
        rmsd_sc = kabsch_rmsd(coords_a_full[sc_mask], coords_b_full[sc_mask])
        n_sc = int(sc_mask.sum())
    else:
        rmsd_sc = float("nan")
        n_sc = 0

    # Total
    rmsd_total = kabsch_rmsd(coords_a_full, coords_b_full)

    # Interface subset (residue-position membership)
    if interface_positions:
        if_mask = np.array([k[0] in interface_positions for k in info])
        if if_mask.any():
            rmsd_if = kabsch_rmsd(coords_a_full[if_mask],
                                  coords_b_full[if_mask])
            n_if = int(if_mask.sum())
        else:
            rmsd_if = float("nan")
            n_if = 0
    else:
        rmsd_if = float("nan")
        n_if = 0

    return rmsd_bb, rmsd_sc, rmsd_total, rmsd_if, n_bb, n_sc, len(idx_a), n_if


# ---------------------------------------------------------------------------
# Main per-TF analysis
# ---------------------------------------------------------------------------
def load_interface_positions(tf_name, interface_dir):
    """Try to load <interface_dir>/<tf>.json. Returns (set_of_positions, json_path)
    or (None, None) if not present."""
    if not interface_dir:
        return None, None
    path = os.path.join(interface_dir, f"{tf_name}.json")
    if not os.path.isfile(path):
        return None, path
    import json
    with open(path) as f:
        data = json.load(f)
    return set(data.get("interface_residue_positions", [])), path


def analyze_tf(tf_name, output_rows, variants, interface_dir=None):
    print(f"\n=== {tf_name} ===")
    cfg = load_pilot_config(tf_name)
    pdb_id = cfg["PDB_ID"]
    protein_chain = int(cfg["PROTEIN_CHAIN"])
    n_frames = int(cfg["N_FRAMES"])
    print(f"  PDB_ID={pdb_id}, PROTEIN_CHAIN={protein_chain}, N_FRAMES={n_frames}")

    # Load reference (full structure; we'll constrain to protein chain at compare time)
    try:
        ref = load_reference_full(pdb_id)
        n_prot_in_ref = sum(
            1 for atom in ref.topology.atoms
            if atom.residue.chain.index == protein_chain
            and atom.residue.name in PROTEIN_RESNAMES
            and (atom.element is None or atom.element.symbol != "H")
        )
        print(f"  Reference: {ref.n_atoms} atoms total, {n_prot_in_ref} protein heavy atoms in chain {protein_chain}")
    except Exception as e:
        print(f"  ERROR loading reference: {e}", file=sys.stderr)
        return

    # Shared stages: 0 (BioEmu), 1 (HPACKER), 2 (docked).
    stage0 = load_stage0(pdb_id)
    print(f"  Stage 0 (BioEmu raw):  {('n_frames=' + str(stage0.n_frames)) if stage0 else 'NOT AVAILABLE'}")
    stage1 = load_stage1(pdb_id)
    print(f"  Stage 1 (HPACKER):     {('n_frames=' + str(stage1.n_frames)) if stage1 else 'NOT AVAILABLE'}")
    stage2_list = load_stage_per_state(pdb_id, tf_name, "stage2_docked", n_frames)
    n2 = sum(1 for t in stage2_list if t is not None)
    print(f"  Stage 2 (docked):      {n2} states present")

    # Per-variant stage 3.
    variant_stage3_lists = {}
    for vname in variants:
        sdir = VARIANTS[vname]
        s3 = load_stage_per_state(pdb_id, tf_name, sdir, n_frames)
        n3 = sum(1 for t in s3 if t is not None)
        print(f"  Stage 3 ({vname:10s}): {n3} states present (from {sdir}/)")
        variant_stage3_lists[vname] = s3

    # Interface residues (optional).
    interface_positions, if_path = load_interface_positions(tf_name, interface_dir)
    if interface_positions is not None:
        print(f"  Interface:             {len(interface_positions)} residues "
              f"(from {if_path})")
    elif if_path is not None:
        print(f"  Interface:             [skipped — {if_path} not found]")

    def _row(stage_n, variant_name, comparison, *, traj_a, frame_a, traj_b,
             frame_b, chain_a=None, chain_b=None):
        bb, sc, tot, ifr, nbb, nsc, ntot, nif = compute_rmsd_triple(
            traj_a, frame_a, traj_b, frame_b,
            chain_a=chain_a, chain_b=chain_b,
            interface_positions=interface_positions,
        )
        return {
            "tf": tf_name, "pdb_id": pdb_id, "state": state_i,
            "stage": stage_n, "variant": variant_name,
            "comparison": comparison,
            "rmsd_backbone_A": bb, "rmsd_sidechain_A": sc,
            "rmsd_total_A": tot, "rmsd_interface_A": ifr,
            "n_backbone_atoms": nbb, "n_sidechain_atoms": nsc,
            "n_total_atoms": ntot, "n_interface_atoms": nif,
        }

    # Iterate over states
    for state_i in range(1, n_frames + 1):
        # Shared stages present for this state
        shared_stages = {}
        if stage0 is not None and state_i - 1 < stage0.n_frames:
            shared_stages[0] = (stage0, state_i - 1)
        if stage1 is not None and state_i - 1 < stage1.n_frames:
            shared_stages[1] = (stage1, state_i - 1)
        if stage2_list[state_i] is not None:
            shared_stages[2] = (stage2_list[state_i], 0)

        # Per-variant stage 3 entries for this state
        variant_stage3 = {}
        for vname in variants:
            if variant_stage3_lists[vname][state_i] is not None:
                variant_stage3[vname] = (variant_stage3_lists[vname][state_i], 0)

        # vs_reference for shared stages 0–2 (variant="shared")
        for stage_n, (traj, frame_i) in sorted(shared_stages.items()):
            output_rows.append(_row(
                stage_n, "shared", "vs_reference",
                traj_a=traj, frame_a=frame_i,
                traj_b=ref, frame_b=0, chain_b=protein_chain,
            ))

        # vs_reference for variant stage 3 (variant=<vname>)
        for vname, (traj, frame_i) in variant_stage3.items():
            output_rows.append(_row(
                3, vname, "vs_reference",
                traj_a=traj, frame_a=frame_i,
                traj_b=ref, frame_b=0, chain_b=protein_chain,
            ))

        # delta_prev within shared stages (0→1, 1→2): variant="shared"
        shared_keys = sorted(shared_stages.keys())
        for i in range(1, len(shared_keys)):
            prev_s = shared_keys[i - 1]
            curr_s = shared_keys[i]
            traj_p, frame_p = shared_stages[prev_s]
            traj_c, frame_c = shared_stages[curr_s]
            output_rows.append(_row(
                curr_s, "shared", f"delta_stage{prev_s}_to_stage{curr_s}",
                traj_a=traj_c, frame_a=frame_c,
                traj_b=traj_p, frame_b=frame_p,
            ))

        # delta stage 2 → stage 3 per variant: variant=<vname>
        if 2 in shared_stages:
            traj_p, frame_p = shared_stages[2]
            for vname, (traj_c, frame_c) in variant_stage3.items():
                output_rows.append(_row(
                    3, vname, "delta_stage2_to_stage3",
                    traj_a=traj_c, frame_a=frame_c,
                    traj_b=traj_p, frame_b=frame_p,
                ))


# ---------------------------------------------------------------------------
# Summary stats
# ---------------------------------------------------------------------------
def compute_summary(rows):
    """Group by (tf, stage, variant, comparison, atom_subset) and compute
    median/IQR/max."""
    from collections import defaultdict
    groups = defaultdict(list)
    for r in rows:
        for subset in ("backbone", "sidechain", "total", "interface"):
            val = r.get(f"rmsd_{subset}_A", float("nan"))
            if not (isinstance(val, float) and np.isnan(val)):
                key = (r["tf"], r["stage"], r["variant"], r["comparison"], subset)
                groups[key].append(val)

    out = []
    for (tf, stage, variant, comp, subset), vals in sorted(groups.items()):
        arr = np.array(vals)
        out.append({
            "tf": tf, "stage": stage, "variant": variant, "comparison": comp,
            "atom_subset": subset,
            "n_states": len(arr),
            "median_A": float(np.median(arr)),
            "q25_A": float(np.percentile(arr, 25)),
            "q75_A": float(np.percentile(arr, 75)),
            "iqr_A": float(np.percentile(arr, 75) - np.percentile(arr, 25)),
            "max_A": float(arr.max()),
            "min_A": float(arr.min()),
        })
    return out


def write_csv(rows, path, columns):
    import csv
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        for r in rows:
            writer.writerow({k: r.get(k, "") for k in columns})


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    global BIOEMU_ROOT, STAGE1_OUT, CONF_ROOT, PILOTS_DIR
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tfs", nargs="+", default=["egr1", "dux4", "tbp"],
                        help="TF names to analyze (default: egr1 dux4 tbp)")
    parser.add_argument("--variants", nargs="+",
                        default=list(VARIANTS),
                        choices=list(VARIANTS),
                        help="Augmentation variants to compute stage-3 RMSDs for "
                             "(default: both metal_cage and legacy).")
    parser.add_argument("--output-dir", default=".", help="Where to write CSVs")
    parser.add_argument("--bioemu-root", default=BIOEMU_ROOT,
                        help=f"Source-chain dir: <pdb>_chains/<pdb>.cif "
                             f"(default: {BIOEMU_ROOT}; env: DEEPPBS_BIOEMU_ROOT).")
    parser.add_argument("--stage1-out", default=STAGE1_OUT,
                        help=f"Stage 1 library: <pdb>_chain<X>_conformations/ "
                             f"(default: {STAGE1_OUT}; env: DEEPPBS_STAGE1_OUT).")
    parser.add_argument("--conf-root", default=CONF_ROOT,
                        help=f"Pipeline output root: <stage>/<tf>/ "
                             f"(default: {CONF_ROOT}; env: DEEPPBS_CONF_ROOT).")
    parser.add_argument("--pilots-dir", default=PILOTS_DIR,
                        help=f"Dir holding <tf>.sh pilot configs "
                             f"(default: {PILOTS_DIR}; env: DEEPPBS_PILOTS_DIR).")
    parser.add_argument("--interface-dir", default=None,
                        help="Optional dir holding per-pilot interface JSONs "
                             "(<tf>.json) produced by "
                             "identify_interface_residues.py. When set, every "
                             "row also carries rmsd_interface_A / "
                             "n_interface_atoms restricted to those residues.")
    args = parser.parse_args()

    # Apply CLI overrides to module-level path globals.
    BIOEMU_ROOT = args.bioemu_root
    STAGE1_OUT  = args.stage1_out
    CONF_ROOT   = args.conf_root
    PILOTS_DIR  = args.pilots_dir

    import_mdtraj()
    warnings.filterwarnings("ignore", category=UserWarning, module="mdtraj")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Variants: {args.variants}")
    if args.interface_dir:
        print(f"Interface dir: {args.interface_dir}")
    all_rows = []
    for tf in args.tfs:
        try:
            analyze_tf(tf, all_rows, args.variants,
                       interface_dir=args.interface_dir)
        except Exception as e:
            print(f"  ERROR analyzing {tf}: {e}", file=sys.stderr)
            import traceback; traceback.print_exc()

    print(f"\nWriting outputs ({len(all_rows)} rows total)")

    per_state_csv = output_dir / "per_state_rmsds.csv"
    write_csv(
        all_rows, per_state_csv,
        columns=["tf", "pdb_id", "state", "stage", "variant", "comparison",
                 "rmsd_backbone_A", "rmsd_sidechain_A", "rmsd_total_A",
                 "rmsd_interface_A",
                 "n_backbone_atoms", "n_sidechain_atoms", "n_total_atoms",
                 "n_interface_atoms"],
    )
    print(f"  Wrote {per_state_csv}")

    summary_rows = compute_summary(all_rows)
    summary_csv = output_dir / "summary_stats.csv"
    write_csv(
        summary_rows, summary_csv,
        columns=["tf", "stage", "variant", "comparison", "atom_subset",
                 "n_states", "median_A", "q25_A", "q75_A", "iqr_A",
                 "min_A", "max_A"],
    )
    print(f"  Wrote {summary_csv}")

    # Pretty-print summary to stdout
    print("\n=== Summary (median RMSD, Å) ===")
    print(f"{'tf':6s} {'stage':5s} {'variant':11s} {'comparison':25s} "
          f"{'subset':10s}  {'n':>4s} {'median':>7s} {'IQR':>7s} {'max':>7s}")
    for r in summary_rows:
        print(f"{r['tf']:6s} {r['stage']:5d} {r['variant']:11s} "
              f"{r['comparison']:25s} {r['atom_subset']:10s}  "
              f"{r['n_states']:4d} {r['median_A']:7.3f} {r['iqr_A']:7.3f} "
              f"{r['max_A']:7.3f}")


if __name__ == "__main__":
    main()
