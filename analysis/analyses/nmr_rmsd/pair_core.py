#!/usr/bin/env python
"""
pair_core.py -- the shared comparison core for an apo/holo NMR pair.

The per-entry script (compute_nmr_state_rmsd.py) aligns each BioEmu ensemble to
its own reference independently, so the apo and holo numbers of a pair land on
different atom sets (57 CA for 1iv6 vs 67 residues in 1ity, ...).  Those are not
comparable, and any apo-vs-holo distribution built from them measures construct
length as much as conformation.

This module builds ONE atom set per pair, shared by all four members

    apo NMR bundle  /  holo NMR bundle  /  BioEmu(apo seq)  /  BioEmu(holo seq)

so every RMSD in the pair is computed on the same residues.  Two trims are
applied on top of the 4-way sequence intersection:

  1. expression tags        -- regex on each member's own sequence (His-tags,
     enterokinase sites).  The intersection alone does not remove these when
     BOTH members carry the tag (5zux/5zuz share a C-term LEHHHHHH).
  2. disordered termini     -- contiguous positions trimmed inward from each end
     while the intra-bundle spread in EITHER NMR bundle exceeds a threshold.
     Only the ends are trimmed, so the core is never punched full of holes.

Alignment helpers (residue_sequence, select_reference_chain, reference_cif,
ATOM_SETS) are imported from compute_nmr_state_rmsd so the two scripts cannot
drift apart.
"""
import csv
import difflib
import re
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
from compute_nmr_state_rmsd import (          # noqa: E402  - reuse, do not duplicate
    ATOM_SETS,
    load_reference,
    reference_cif,
    residue_sequence,
    select_reference_chain,
)

_TFCONF = _HERE.parents[2]
NMR_DIR = _TFCONF / "structures" / "nmr"
CIF_SEARCH_DIRS = [
    NMR_DIR,
    NMR_DIR / "rcsb_downloads",
    _TFCONF / "structures" / "other_proteins" / "rcsb_downloads",
    _TFCONF / "structures" / "source_chains",
]

# Purification tags / cloning scars. Matched against each member's own
# one-letter sequence; every residue inside a hit is dropped from the core.
TAG_PATTERNS = (
    r"H{4,}",          # poly-His
    r"DDDDK",          # enterokinase site
    r"ENLYFQ[GS]?",     # TEV site
)

MEMBERS = ("holo_ref", "apo_ref", "holo_ens", "apo_ens")   # holo_ref is the anchor


# ---------------------------------------------------------------------------
# pair metadata
# ---------------------------------------------------------------------------
def load_pairs(path, pair_ids=None, statuses=("ready",)):
    """Read pairs.csv (# comments allowed) -> list of dicts.

    Rows carry a `status`: "ready" (both members have a BioEmu ensemble on disk)
    or "needs_ensembles" (pair found by PDB search, nothing generated yet). Only
    `statuses` are returned, so the candidate backlog can live in the same file
    without every script skipping 30-odd rows on each run.

    Naming pairs explicitly with pair_ids overrides the status filter -- that is
    how you analyse a pair the moment its ensembles land, before editing the CSV.
    """
    with open(path) as fh:
        rows = list(csv.DictReader(l for l in fh if not l.lstrip().startswith("#")))
    if pair_ids:
        return [r for r in rows if r["pair_id"] in set(pair_ids)]
    if statuses is not None:
        rows = [r for r in rows if r.get("status", "ready") in set(statuses)]
    return rows


def ensemble_paths(pdb_id, chain):
    """(topology, xtc) for <pdb>_chain<X>_conformations/, matching either naming."""
    stem = f"{pdb_id}_chain{chain}"
    d = NMR_DIR / f"{stem}_conformations"
    top = next((p for p in (d / f"{stem}_topology.pdb", d / "topology.pdb") if p.is_file()), None)
    xtc = next((p for p in (d / f"{stem}_samples.xtc", d / "samples.xtc") if p.is_file()), None)
    if top is None or xtc is None:
        raise FileNotFoundError(f"{d}: missing topology.pdb or samples.xtc")
    return top, xtc


# ---------------------------------------------------------------------------
# alignment
# ---------------------------------------------------------------------------
def anchor_map(anchor_seq, other_seq):
    """{anchor position -> other position} from difflib's matching blocks."""
    matcher = difflib.SequenceMatcher(a=anchor_seq, b=other_seq, autojunk=False)
    out = {}
    for i, j, n in matcher.get_matching_blocks():
        for k in range(n):
            out[i + k] = j + k
    return out


def tagged_positions(seq):
    """Positions of seq lying inside any TAG_PATTERNS hit."""
    bad = set()
    for pat in TAG_PATTERNS:
        for m in re.finditer(pat, seq):
            bad.update(range(m.start(), m.end()))
    return bad


def build_core(member_res, atom_names):
    """4-way intersection of residues + the atom indices realising it.

    member_res : {member -> [(residue, one-letter), ...]}, anchored on MEMBERS[0].

    Returns (core, atom_idx) where
      core     : [{member -> residue}, ...]  one entry per surviving position
      atom_idx : {member -> (n_pos, n_atoms_per_pos) int array}
    A position survives only if every member has it AND every member's residue
    carries every atom in atom_names -- so the per-position atom count is
    constant and the arrays are rectangular (needed for the per-position spread).
    """
    anchor = MEMBERS[0]
    seqs = {m: "".join(aa for _, aa in res) for m, res in member_res.items()}
    maps = {m: anchor_map(seqs[anchor], seqs[m]) for m in MEMBERS if m != anchor}
    maps[anchor] = {i: i for i in range(len(seqs[anchor]))}

    masked = {m: tagged_positions(seqs[m]) for m in MEMBERS}

    core, idx = [], {m: [] for m in MEMBERS}
    for pos in range(len(seqs[anchor])):
        if any(pos not in maps[m] for m in MEMBERS):
            continue
        if any(maps[m][pos] in masked[m] for m in MEMBERS):
            continue
        res = {m: member_res[m][maps[m][pos]][0] for m in MEMBERS}
        atoms = {}
        for m in MEMBERS:
            by_name = {a.name: a.index for a in res[m].atoms}
            if not all(n in by_name for n in atom_names):
                atoms = None
                break
            atoms[m] = [by_name[n] for n in atom_names]
        if atoms is None:
            continue
        core.append(res)
        for m in MEMBERS:
            idx[m].append(atoms[m])
    return core, {m: np.array(v, dtype=int) for m, v in idx.items()}


# ---------------------------------------------------------------------------
# terminal trimming
# ---------------------------------------------------------------------------
def position_spread(traj, atom_idx):
    """Per-position RMS deviation about the bundle mean, in Angstrom.

    atom_idx : (n_pos, n_atoms_per_pos). Superposes the bundle on the whole core
    first, so the spread is the internal heterogeneity of the deposited models.
    """
    flat = atom_idx.reshape(-1)
    traj = traj.slice(range(traj.n_frames), copy=True)
    traj.superpose(traj, frame=0, atom_indices=flat, ref_atom_indices=flat)
    xyz = traj.xyz[:, flat].reshape(traj.n_frames, atom_idx.shape[0], atom_idx.shape[1], 3)
    dev = xyz - xyz.mean(axis=0, keepdims=True)
    return np.sqrt((dev ** 2).sum(-1).mean(axis=(0, 2))) * 10.0


def trim_flexible_termini(spreads, threshold):
    """Contiguous inward trim from both ends while max spread > threshold.

    spreads : (n_members_to_check, n_pos). Returns a boolean keep mask.
    """
    worst = np.max(spreads, axis=0)
    n = len(worst)
    lo = 0
    while lo < n and worst[lo] > threshold:
        lo += 1
    hi = n
    while hi > lo and worst[hi - 1] > threshold:
        hi -= 1
    keep = np.zeros(n, dtype=bool)
    keep[lo:hi] = True
    return keep


# ---------------------------------------------------------------------------
# assembly
# ---------------------------------------------------------------------------
def load_pair(pair, atoms="ca", spread_cutoff=3.0, fetch=True, verbose=True):
    """Load all four members of a pair and build their shared core.

    Returns a dict with the trajectories, the core atom indices (post-trim),
    the per-position audit table, and the bookkeeping the callers report.
    """
    import mdtraj as md

    atom_names = ATOM_SETS[atoms]
    cache = NMR_DIR / "cif"
    out = {"pair": pair, "atoms": atoms, "atom_names": atom_names}

    traj, res = {}, {}
    for role, pdb_key, chain_key in (("holo", "holo_pdb", "holo_chain"),
                                     ("apo", "apo_pdb", "apo_chain")):
        pdb_id, chain = pair[pdb_key], pair[chain_key]
        top_path, xtc_path = ensemble_paths(pdb_id, chain)
        ens = md.load(str(xtc_path), top=str(top_path))
        ens_res = residue_sequence(ens.top)
        ens_seq = "".join(aa for _, aa in ens_res)

        cif = reference_cif(pdb_id, CIF_SEARCH_DIRS, cache, fetch=fetch)
        ref = load_reference(cif, verbose=verbose)
        ref_res, ref_chain, identity = select_reference_chain(ref.top, chain, ens_seq)

        traj[f"{role}_ens"], res[f"{role}_ens"] = ens, ens_res
        traj[f"{role}_ref"], res[f"{role}_ref"] = ref, ref_res
        out[f"{role}_ref_chain"] = ref_chain
        out[f"{role}_ref_identity"] = identity
        out[f"{role}_n_states"] = ref.n_frames
        out[f"{role}_n_frames"] = ens.n_frames
        if verbose:
            print(f"  {role:4s} {pdb_id} chain {chain}: {ref.n_frames} NMR states "
                  f"(ref chain {ref_chain}, id {identity:.2f}), {ens.n_frames} BioEmu frames")

    core, idx = build_core(res, atom_names)
    if len(core) < 4:
        raise ValueError(f"only {len(core)} positions survive the 4-way intersection")
    n_intersect = len(core)

    # spread is judged on the two NMR bundles only -- the BioEmu ensembles are
    # what we are testing, so they must not define the core.
    spreads = np.vstack([position_spread(traj[m], idx[m]) for m in ("apo_ref", "holo_ref")])
    keep = trim_flexible_termini(spreads, spread_cutoff)
    if keep.sum() < 4:
        raise ValueError(f"only {int(keep.sum())} positions survive the terminal trim")

    rows = []
    for k, res_map in enumerate(core):
        rows.append({
            "position": k,
            "kept": int(keep[k]),
            "apo_spread_A": round(float(spreads[0, k]), 3),
            "holo_spread_A": round(float(spreads[1, k]), 3),
            **{f"{m}_resid": res_map[m].resSeq for m in MEMBERS},
            **{f"{m}_resname": res_map[m].name for m in MEMBERS},
        })

    out["traj"] = traj
    out["atom_idx"] = {m: idx[m][keep].reshape(-1) for m in MEMBERS}
    out["core_residues"] = rows
    out["n_intersect"] = n_intersect
    out["n_core"] = int(keep.sum())
    out["n_trimmed"] = n_intersect - int(keep.sum())
    out["n_atoms"] = int(out["atom_idx"]["holo_ref"].size)
    if verbose:
        print(f"  core: {n_intersect} positions from the 4-way intersection, "
              f"{out['n_trimmed']} terminal positions trimmed (> {spread_cutoff} A "
              f"bundle spread) -> {out['n_core']} positions / {out['n_atoms']} {atoms} atoms")
    return out


def rmsd_matrix(mobile, reference, mob_idx, ref_idx):
    """(n_ref_frames, n_mobile_frames) optimally-superposed RMSD in Angstrom."""
    import mdtraj as md
    m = np.empty((reference.n_frames, mobile.n_frames))
    for f in range(reference.n_frames):
        m[f] = md.rmsd(mobile, reference, frame=f,
                       atom_indices=mob_idx, ref_atom_indices=ref_idx) * 10.0
    return m
