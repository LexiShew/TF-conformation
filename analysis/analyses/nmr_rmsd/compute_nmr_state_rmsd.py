#!/usr/bin/env python
"""
RMSD between every NMR state of a reference mmCIF and every BioEmu conformation
of the matching chain.

For each <PDB>_chain<X>_conformations/ directory under --nmr-dir we load

  reference : <PDB>.cif            -- all NMR models, chain <X>
  ensemble  : topology.pdb + samples.xtc  (BioEmu, backbone+CB only)

pair up residues by sequence (difflib, so ragged termini / numbering offsets do
not matter), and compute the optimally-superposed RMSD of every
(state, conformation) pair.  Reference mmCIFs missing locally are downloaded
from RCSB unless --no-fetch is given.

Outputs (in --output-dir):
  nmr_state_rmsd_pairs_<atoms>.csv  pdb_id,chain,ref_chain,state,conformation,
                            n_atoms,rmsd_A -- one row per (state, conformation)
  nmr_state_rmsd_best_<atoms>.csv   pdb_id,chain,ref_chain,conformation,
                            best_state,rmsd_A,n_states,n_atoms
                            -- per conformation, its closest NMR state

Usage:
    python compute_nmr_state_rmsd.py                      # all entries, CA atoms
    python compute_nmr_state_rmsd.py --atoms backbone     # N,CA,C,O
    python compute_nmr_state_rmsd.py --pdb-ids 1nk2 2stt
"""
import argparse
import csv
import difflib
import os
import sys
import urllib.request
from pathlib import Path

import numpy as np

_TFCONF = Path(__file__).resolve().parents[3]

ATOM_SETS = {
    "ca": ("CA",),
    "backbone": ("N", "CA", "C", "O"),
}

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEP": "S", "TPO": "T", "PTR": "Y", "HYP": "P",
    "HID": "H", "HIE": "H", "HIP": "H", "HSD": "H", "HSE": "H", "HSP": "H",
}


def load_reference(path, verbose=True):
    """Load a multi-model NMR reference, tolerating ragged depositions.

    A few entries (1osl among them) deposit models with DIFFERENT atom counts,
    which both the mmCIF and the PDB reader reject outright:
        "Atom N for model 11 does not match the order of atoms for model 1"
    That is a property of the file, not of the parser. Fall back to reading the
    models one at a time and keeping only the atoms present in EVERY model,
    ordered as model 1 has them, so the bundle becomes rectangular and every
    downstream RMSD still compares like with like.
    """
    import mdtraj as md
    try:
        return md.load(str(path))
    except (ValueError, IndexError) as exc:
        if verbose:
            print(f"  note: {Path(path).name} is ragged ({str(exc)[:60]}...); "
                  f"rebuilding from the atoms common to all models", file=sys.stderr)

    pdb = _as_pdb(path)
    models = _split_models(pdb)
    if len(models) < 2:
        raise ValueError(f"{path}: unreadable and not splittable into models")

    import tempfile
    trajs = []
    for text in models:
        with tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=False) as fh:
            fh.write(text)
            tmp = fh.name
        try:
            trajs.append(md.load(tmp))
        except Exception:                      # noqa: BLE001 - drop unreadable models
            pass
        finally:
            os.unlink(tmp)
    if not trajs:
        raise ValueError(f"{path}: no model could be read")

    def keys(t):
        return [(a.residue.chain.chain_id, a.residue.resSeq, a.residue.name, a.name)
                for a in t.top.atoms]
    common = set(keys(trajs[0]))
    for t in trajs[1:]:
        common &= set(keys(t))
    if len(common) < 4:
        raise ValueError(f"{path}: models share only {len(common)} atoms")

    order = [k for k in keys(trajs[0]) if k in common]
    sliced = []
    for t in trajs:
        idx = {k: i for i, k in enumerate(keys(t))}
        sliced.append(t.atom_slice([idx[k] for k in order]))
    out = sliced[0].join(sliced[1:]) if len(sliced) > 1 else sliced[0]
    if verbose:
        print(f"  note: kept {out.n_frames}/{len(models)} models, "
              f"{out.n_atoms} atoms common to all", file=sys.stderr)
    return out


def _as_pdb(path):
    """Text of the legacy PDB for this entry, fetching it if only the cif is local."""
    path = Path(path)
    cand = path.with_suffix(".pdb")
    if cand.is_file():
        return cand.read_text()
    pdb_id = path.stem.lower()
    dest = path.parent / f"{pdb_id}.pdb"
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    urllib.request.urlretrieve(url, dest)
    return dest.read_text()


def _split_models(text):
    """Split PDB text on MODEL/ENDMDL into standalone single-model PDB strings."""
    head = [l for l in text.splitlines(True)
            if l.startswith(("CRYST1", "SCALE", "ORIGX", "SEQRES", "HETATM_TEMPLATE"))]
    models, cur, inside = [], [], False
    for line in text.splitlines(True):
        if line.startswith("MODEL "):
            cur, inside = [], True
            continue
        if line.startswith("ENDMDL"):
            if cur:
                models.append("".join(head + cur) + "END\n")
            cur, inside = [], False
            continue
        if inside and line.startswith(("ATOM", "HETATM", "TER", "ANISOU")):
            cur.append(line)
    return models



def find_ensembles(nmr_dir, pdb_ids=None):
    """Yield (pdb_id, chain, topology_path, xtc_path) for each conformation dir."""
    for d in sorted(Path(nmr_dir).glob("*_chain*_conformations")):
        stem = d.name[: -len("_conformations")]          # <pdb>_chain<X>
        pdb_id, chain = stem.split("_chain", 1)
        if pdb_ids and pdb_id not in pdb_ids:
            continue
        # Accept both the original names and the <pdb>_chain<X>_ prefixed ones.
        top = next((p for p in (d / f"{stem}_topology.pdb", d / "topology.pdb") if p.is_file()), None)
        xtc = next((p for p in (d / f"{stem}_samples.xtc", d / "samples.xtc") if p.is_file()), None)
        if top is None or xtc is None:
            print(f"  SKIP {d.name}: missing topology.pdb or samples.xtc", file=sys.stderr)
            continue
        yield pdb_id, chain, top, xtc


def reference_cif(pdb_id, search_dirs, cache_dir, fetch=True):
    """Locate <pdb_id>.cif, downloading it from RCSB into cache_dir if needed."""
    for d in search_dirs:
        for cand in (Path(d) / f"{pdb_id}.cif", Path(d) / f"{pdb_id}_chains" / f"{pdb_id}.cif"):
            if cand.is_file():
                return cand
    dest = Path(cache_dir) / f"{pdb_id}.cif"
    if dest.is_file() and dest.stat().st_size:
        return dest
    if not fetch:
        raise FileNotFoundError(f"no local {pdb_id}.cif and --no-fetch was given")
    dest.parent.mkdir(parents=True, exist_ok=True)
    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    print(f"  fetching {url}")
    tmp = dest.with_suffix(".cif.part")
    urllib.request.urlretrieve(url, tmp)
    tmp.rename(dest)
    return dest


def residue_sequence(chain_or_top):
    """[(residue, one-letter)] for the protein residues of a chain/topology."""
    out = []
    for res in chain_or_top.residues:
        aa = THREE_TO_ONE.get(res.name.upper())
        if aa is None:                     # DNA/RNA, waters, ligands
            continue
        out.append((res, aa))
    return out


def select_reference_chain(top, wanted_chain, mob_seq):
    """Pick the reference protein chain to compare against.

    mdtraj reports mmCIF label_asym_id, which need not equal the auth chain id
    our FASTA/directory names use, and NMR entries are often homodimers whose
    chains match the ensemble equally well.  So: score every protein chain by
    sequence similarity to the ensemble, keep the best score, and among ties
    prefer the chain whose id is the one we asked for.

    Returns (residues, chain_label, similarity).
    """
    scored = []
    for ch in top.chains:
        res = residue_sequence(ch)
        if not res:
            continue
        seq = "".join(aa for _, aa in res)
        ratio = difflib.SequenceMatcher(a=seq, b=mob_seq, autojunk=False).ratio()
        scored.append((ratio, getattr(ch, "chain_id", None) or str(ch.index), res))
    if not scored:
        raise ValueError("no protein chain in the reference structure")
    best = max(r for r, _, _ in scored)
    tied = [s for s in scored if best - s[0] < 1e-9]
    for ratio, label, res in tied:
        if label == wanted_chain:
            return res, label, ratio
    ratio, label, res = tied[0]
    return res, label, ratio


def paired_atom_indices(ref_res, mob_res, atom_names):
    """Sequence-align two residue lists; return matched atom index arrays."""
    ref_seq = "".join(aa for _, aa in ref_res)
    mob_seq = "".join(aa for _, aa in mob_res)
    matcher = difflib.SequenceMatcher(a=ref_seq, b=mob_seq, autojunk=False)
    ref_idx, mob_idx = [], []
    for i, j, n in matcher.get_matching_blocks():
        for k in range(n):
            r, m = ref_res[i + k][0], mob_res[j + k][0]
            r_atoms = {a.name: a.index for a in r.atoms}
            m_atoms = {a.name: a.index for a in m.atoms}
            # all-or-nothing: a residue contributes every atom or none, so the
            # index arrays stay rectangular and can be reshaped per residue.
            if not all(n in r_atoms and n in m_atoms for n in atom_names):
                continue
            for name in atom_names:
                ref_idx.append(r_atoms[name])
                mob_idx.append(m_atoms[name])
    identity = matcher.ratio()
    return np.array(ref_idx, dtype=int), np.array(mob_idx, dtype=int), identity


def trim_termini(ref, ref_idx, mob_idx, n_per_res, cutoff):
    """Drop terminal residues whose spread within the reference bundle exceeds cutoff.

    Same rule the pair analysis uses (pair_core.trim_flexible_termini), applied
    here to a single bundle. Without it the RMSD of an entry with long floppy
    tails measures the tails: 1j5n reads ~6 A untrimmed and ~1 A on its core.
    Only contiguous positions at each end are removed, so interior mobile
    regions stay in.
    """
    r = ref_idx.reshape(-1, n_per_res)
    m = mob_idx.reshape(-1, n_per_res)
    flat = r.reshape(-1)
    t = ref.slice(range(ref.n_frames), copy=True)
    t.superpose(t, frame=0, atom_indices=flat, ref_atom_indices=flat)
    xyz = t.xyz[:, flat].reshape(t.n_frames, r.shape[0], n_per_res, 3)
    dev = xyz - xyz.mean(axis=0, keepdims=True)
    spread = np.sqrt((dev ** 2).sum(-1).mean(axis=(0, 2))) * 10.0

    n = len(spread)
    lo = 0
    while lo < n and spread[lo] > cutoff:
        lo += 1
    hi = n
    while hi > lo and spread[hi - 1] > cutoff:
        hi -= 1
    if hi - lo < 3:                      # refuse to trim away the whole thing
        return ref_idx, mob_idx, n, 0
    return r[lo:hi].reshape(-1), m[lo:hi].reshape(-1), n, n - (hi - lo)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nmr-dir", default=_TFCONF / "structures" / "nmr",
                    help="directory holding <PDB>_chain<X>_conformations/ (default: %(default)s)")
    ap.add_argument("--cif-dir", action="append", default=None,
                    help="extra directory to search for <PDB>.cif (repeatable)")
    ap.add_argument("--cache-dir", default=None,
                    help="where fetched mmCIFs are stored (default: <nmr-dir>/cif)")
    ap.add_argument("--output-dir", default=Path(__file__).resolve().parent)
    ap.add_argument("--atoms", choices=sorted(ATOM_SETS), default="ca",
                    help="atom subset used for superposition + RMSD (default: %(default)s); "
                         "BioEmu output is backbone+CB only, so sidechains are never comparable")
    ap.add_argument("--pdb-ids", nargs="*", default=None, help="restrict to these PDB IDs")
    ap.add_argument("--no-fetch", action="store_true", help="never download from RCSB")
    ap.add_argument("--trim-termini", type=float, default=None, metavar="A",
                    help="drop terminal residues whose spread within the NMR bundle exceeds "
                         "this many Angstrom (3.0 matches the pair analysis). Off by default; "
                         "without it the RMSD of a floppy-tailed entry measures its tails")
    ap.add_argument("--label", default="", help="suffix added to output filenames")
    args = ap.parse_args()

    import mdtraj as md   # imported late so --help works without the env

    nmr_dir = Path(args.nmr_dir)
    cache_dir = Path(args.cache_dir) if args.cache_dir else nmr_dir / "cif"
    search_dirs = [nmr_dir, _TFCONF / "structures" / "source_chains", cache_dir]
    search_dirs += [Path(d) for d in (args.cif_dir or [])]
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    atom_names = ATOM_SETS[args.atoms]

    pair_rows, best_rows = [], []
    for pdb_id, chain, top_path, xtc_path in find_ensembles(nmr_dir, args.pdb_ids):
        print(f"{pdb_id} chain {chain}")
        try:
            cif = reference_cif(pdb_id, search_dirs, cache_dir, fetch=not args.no_fetch)
        except Exception as exc:                      # noqa: BLE001 - report and continue
            print(f"  SKIP: {exc}", file=sys.stderr)
            continue

        ref = load_reference(cif)
        ens = md.load(str(xtc_path), top=str(top_path))

        mob_res = residue_sequence(ens.top)
        mob_seq = "".join(aa for _, aa in mob_res)
        try:
            ref_res, ref_chain, identity = select_reference_chain(ref.top, chain, mob_seq)
        except ValueError as exc:
            print(f"  SKIP: {exc}", file=sys.stderr)
            continue
        if ref_chain != chain:
            print(f"  note: matched reference chain {ref_chain} (asked for {chain})",
                  file=sys.stderr)

        ref_idx, mob_idx, _ = paired_atom_indices(ref_res, mob_res, atom_names)
        n_trimmed = 0
        if args.trim_termini is not None and len(ref_idx) >= 3 * len(atom_names):
            ref_idx, mob_idx, n_res, n_trimmed = trim_termini(
                ref, ref_idx, mob_idx, len(atom_names), args.trim_termini)
            print(f"  trimmed {n_trimmed}/{n_res} terminal residues "
                  f"(> {args.trim_termini} A bundle spread)")
        if len(ref_idx) < 3:
            print(f"  SKIP: only {len(ref_idx)} matched atoms", file=sys.stderr)
            continue
        print(f"  {ref.n_frames} states x {ens.n_frames} conformations, "
              f"{len(ref_idx)} {args.atoms} atoms matched (seq identity {identity:.2f})")

        # md.rmsd superposes on the given atoms; one call per state covers all
        # conformations at once.  Returns nm -> convert to Angstrom.
        rmsd = np.empty((ref.n_frames, ens.n_frames))
        for state in range(ref.n_frames):
            rmsd[state] = md.rmsd(ens, ref, frame=state,
                                  atom_indices=mob_idx, ref_atom_indices=ref_idx) * 10.0

        for state in range(ref.n_frames):
            for conf in range(ens.n_frames):
                pair_rows.append({
                    "pdb_id": pdb_id, "chain": chain, "ref_chain": ref_chain,
                    "state": state + 1, "conformation": conf,
                    "n_atoms": len(ref_idx), "rmsd_A": f"{rmsd[state, conf]:.4f}",
                })
        best_state = rmsd.argmin(axis=0)
        for conf in range(ens.n_frames):
            best_rows.append({
                "pdb_id": pdb_id, "chain": chain, "ref_chain": ref_chain,
                "conformation": conf,
                "best_state": int(best_state[conf]) + 1,
                "rmsd_A": f"{rmsd[best_state[conf], conf]:.4f}",
                "n_states": ref.n_frames, "n_atoms": len(ref_idx),
            })

    if not pair_rows:
        sys.exit("No (state, conformation) pairs computed.")

    tag = f"{args.atoms}{('_' + args.label) if args.label else ''}"
    pairs_csv = out_dir / f"nmr_state_rmsd_pairs_{tag}.csv"
    best_csv = out_dir / f"nmr_state_rmsd_best_{tag}.csv"
    for path, rows in ((pairs_csv, pair_rows), (best_csv, best_rows)):
        with open(path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        print(f"wrote {path}  ({len(rows)} rows)")


if __name__ == "__main__":
    main()
