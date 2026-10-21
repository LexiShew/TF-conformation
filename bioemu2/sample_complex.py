#!/usr/bin/env python
"""
sample_complex.py — sample a protein-DNA COMPLEX ensemble with BioEmu2.

Why this is not a Stage 1 replacement
------------------------------------
Stage 1 (BioEmu 1) samples one PROTEIN CHAIN at a time, backbone-only, from its
sequence. Stage 2 then Kabsch-docks each frame onto the crystal DNA, and Stage 3
minimizes with the DNA backbone restrained. That three-stage dance exists because
BioEmu 1 cannot model nucleic acids.

BioEmu2 can (entity_type "dna"), and it is all-atom. So this script does in one
step what stages 1+2 do in two, and it removes the assumption that item #1 of
docs/concerns.txt flags: that every frame is paired with the same crystal-pose
DNA, making the result "a protein ensemble around a fixed binding pose" rather
than an ensemble of binding modes.

It therefore writes to its OWN tree (output/bioemu2/) and does NOT feed stage2.
This is a parallel track for evaluating whether that reframing holds up, not a
drop-in swap. See bioemu2/README.md.

What it reuses instead of reimplementing
---------------------------------------
* System definition      BioEmu2's System.from_cif() reads the crystal .cif this
                         repo already keeps at
                         structures/source_chains/<PDB>_chains/<PDB>.cif, so
                         protein chains, both DNA strands and structural metals
                         all come across in one call. No YAML to build.
* PDB topology            BioEmu2 hands back an mdtraj.Trajectory, so writing the
                         topology.pdb that this repo's downstream tools expect is
                         traj[0].save_pdb(). No cif->pdb converter needed.
* Structure QC            BioEmu2's own check_physicality().
* Interface scoring       Frames are written as <PDB>_state_NNN.pdb, which is
                         exactly the glob fnat_gate/score_stage3.py already
                         expects, so the existing fnat gate scores these with no
                         new code. That naming is deliberate.

Usage
-----
    python bioemu2/sample_complex.py --tf-name ets1 --num-samples 100
    python bioemu2/sample_complex.py --ref structures/source_chains/1k79_chains/1k79.cif \
        --out output/bioemu2/1k79_complex --num-samples 100

    # drop the structural metals, to A/B whether modelling Zn matters
    python bioemu2/sample_complex.py --tf-name egr1 --exclude-ligands ZN

Outputs (in --out)
------------------
    topology.cif            BioEmu2 native (save_frames)
    samples.xtc             all frames
    topology.pdb            same topology as PDB, for mdtraj/3DNA/PyMOL here
    states/<PDB>_state_NNN.pdb   per-frame, named for fnat_gate compatibility
    physicality.json        BioEmu2's own per-frame validity report
    sample_manifest.json    what was run: checkpoint, seed, counts, timings
"""
import argparse
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

_HERE = Path(__file__).resolve().parent
_TFCONF = _HERE.parent

PILOTS_DIR = _TFCONF / "config" / "pilots"
SOURCE_CHAINS = _TFCONF / "structures" / "source_chains"
DEFAULT_OUT_ROOT = _TFCONF / "output" / "bioemu2"


def pilot_config(tf_name):
    """Read a pilot config by sourcing it, the way validate_pilot.py does.

    Sourcing rather than regex-parsing because configs are generated and may
    expand values; see scripts/pipeline/validate_pilot.py:load_config.
    """
    cfg = PILOTS_DIR / f"{tf_name}.sh"
    if not cfg.exists():
        sys.exit(f"ERROR: no pilot config at {cfg}\n"
                 f"       available: "
                 f"{', '.join(sorted(p.stem for p in PILOTS_DIR.glob('*.sh')))}")
    out = subprocess.run(["bash", "-c", f'set -a; source "{cfg}" >/dev/null 2>&1; env'],
                         capture_output=True, text=True, cwd=_TFCONF)
    env = {}
    for line in out.stdout.splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            env[k] = v
    return env


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--tf-name", help="pilot name; resolves PDB_ID and the "
                                       "reference .cif from config/pilots/")
    src.add_argument("--ref", help="path to a crystal .cif directly")
    p.add_argument("--out", help="output dir (default: output/bioemu2/<pdb>_complex)")
    p.add_argument("--num-samples", type=int, default=100)
    p.add_argument("--seed", type=int, default=1,
                   help="RNG seed. Set it: unseeded runs are what finding F6 was "
                        "about (docs/codebase_review_2026-10-15.md).")
    p.add_argument("--checkpoint", default=None,
                   help="BioEmu2 checkpoint name; default is the package default")
    p.add_argument("--exclude-ligands", nargs="*", default=None,
                   help="CCD codes to drop from the system, e.g. ZN. Use to A/B "
                        "whether modelling a structural metal changes the "
                        "ensemble - this repo currently approximates that with "
                        "Stage 3's coordination cage instead.")
    p.add_argument("--no-ligands", action="store_true",
                   help="drop ALL heterogens (closest to the stages 1-3 behaviour, "
                        "where removeHeterogens strips everything)")
    p.add_argument("--no-states", action="store_true",
                   help="skip the per-frame states/ PDBs (they are what makes the "
                        "existing fnat gate reusable, so normally keep them)")
    a = p.parse_args()

    # ---- resolve inputs -------------------------------------------------
    if a.tf_name:
        env = pilot_config(a.tf_name)
        pdb_id = (env.get("PDB_ID") or "").lower()
        if not pdb_id:
            sys.exit(f"ERROR: {a.tf_name} config has no PDB_ID")
        ref = SOURCE_CHAINS / f"{pdb_id}_chains" / f"{pdb_id}.cif"
        label = a.tf_name
    else:
        ref = Path(a.ref)
        pdb_id = ref.stem.lower()
        label = pdb_id
    if not ref.exists():
        sys.exit(f"ERROR: reference cif not found: {ref}\n"
                 f"       Stage 1's source-chain library should hold it; see "
                 f"stage1_bioemu/README.md for the expected layout.")

    out = Path(a.out) if a.out else DEFAULT_OUT_ROOT / f"{pdb_id}_complex"
    out.mkdir(parents=True, exist_ok=True)

    # ---- import late, so --help works without the bioemu2 env -----------
    try:
        import mdtraj as md
        from bioemu2 import EmuSession, System, load_model, save_frames
    except ImportError as e:
        sys.exit(f"ERROR: {e}\n"
                 f"       BioEmu2 needs its own env (Python 3.12+); the repo's "
                 f"'bioemu' env is 3.11. See bioemu2/README.md.")

    print(f"[bioemu2/{label}] reference : {ref}")
    print(f"[bioemu2/{label}] output    : {out}")
    print(f"[bioemu2/{label}] samples   : {a.num_samples}  seed={a.seed}")

    # ---- build the system straight from the crystal ---------------------
    kw = {}
    if a.no_ligands:
        kw["include_ligands"] = False
    elif a.exclude_ligands:
        kw["exclude_ligands"] = list(a.exclude_ligands)
    system = System.from_cif(str(ref), **kw)

    # Report what came across, so a silently-dropped DNA strand or metal is
    # visible rather than discovered three stages later.
    print(f"[bioemu2/{label}] system entities:")
    by_type = {}
    for ent in system.entities:
        et = getattr(ent, "entity_type", type(ent).__name__.lower())
        by_type.setdefault(str(et), []).append(ent)
        seq = getattr(ent, "sequence", None)
        ccd = getattr(ent, "ccd_code", None)
        asym = getattr(ent, "asym_ids", "?")
        detail = (f"len={len(seq)}" if seq else (f"ccd={ccd}" if ccd else ""))
        print(f"    {str(et):<8} asym={asym!s:<10} {detail}")
    if "dna" not in by_type:
        print(f"[bioemu2/{label}] WARNING: no DNA entity in the system. The whole "
              f"point of this track is sampling the COMPLEX; a protein-only "
              f"system is just a slower Stage 1.", file=sys.stderr)
    n_lig = len(by_type.get("ligand", []))
    print(f"[bioemu2/{label}] {len(system.entities)} entities "
          f"({len(by_type.get('protein', []))} protein, "
          f"{len(by_type.get('dna', []))} dna, {n_lig} ligand)")

    # ---- sample ---------------------------------------------------------
    t0 = time.time()
    model = load_model(a.checkpoint) if a.checkpoint else load_model()
    t_load = time.time() - t0

    try:
        import torch
        torch.manual_seed(a.seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(a.seed)
    except ImportError:
        pass

    t1 = time.time()
    emu = EmuSession(system, model=model, mode="equilibrium")
    # EmuSession.sample() takes the seed itself (bioemu2 0.1 preview); pass it
    # there as well as to torch above, so the RNG that actually drives the
    # reverse diffusion is the seeded one.
    samples = emu.sample(num_samples=a.num_samples, seed=a.seed)
    t_sample = time.time() - t1
    print(f"[bioemu2/{label}] sampled {a.num_samples} in {t_sample/60:.1f} min "
          f"(model load {t_load:.0f}s)")

    # ---- physicality: BioEmu2's own QC, not a reimplementation ----------
    report = None
    try:
        from bioemu2 import check_physicality
        report = check_physicality(samples)
        # PhysicalityReport.passed is a numpy bool array, so `x or []` on it
        # raises "truth value of an array ... is ambiguous" — iterate instead.
        n_pass = int(sum(bool(x) for x in report.passed))
        print(f"[bioemu2/{label}] physicality: {n_pass}/{a.num_samples} frames passed "
              f"(mode={report.mode}, checked={report.checked})")
    except Exception as e:                                  # noqa: BLE001
        print(f"[bioemu2/{label}] physicality check unavailable ({e}); continuing",
              file=sys.stderr)

    # ---- write ----------------------------------------------------------
    save_frames(samples, out, **({"physicality": report} if report else {}))
    print(f"[bioemu2/{label}] wrote topology.cif + samples.xtc")

    # topology.pdb: samples is an mdtraj.Trajectory, so this is a one-liner.
    # Needed because this repo's downstream tools (mdtraj in stage2/fnat_gate,
    # 3DNA in stage4, PyMOL in the figure scripts) all read PDB.
    samples[0].save_pdb(str(out / "topology.pdb"))
    print(f"[bioemu2/{label}] wrote topology.pdb (for mdtraj/3DNA/PyMOL here)")

    # Per-frame PDBs named <PDB>_state_NNN.pdb. This exact pattern is what
    # fnat_gate/score_stage3.py globs, so the existing gate scores these
    # unchanged - that is why the naming matches rather than being prettier.
    if not a.no_states:
        sdir = out / "states"
        sdir.mkdir(exist_ok=True)
        for i in range(samples.n_frames):
            samples[i].save_pdb(str(sdir / f"{pdb_id}_state_{i+1:03d}.pdb"))
        print(f"[bioemu2/{label}] wrote {samples.n_frames} states -> {sdir}")
        print(f"[bioemu2/{label}] score them with the EXISTING gate:")
        print(f"    python fnat_gate/score_stage3.py --ref {ref} \\\n"
              f"        --dir {sdir} --floor 0.5")

    manifest = {
        "tf_name": a.tf_name, "pdb_id": pdb_id, "reference_cif": str(ref),
        "num_samples": a.num_samples, "seed": a.seed,
        "checkpoint": a.checkpoint or "(package default)",
        "exclude_ligands": a.exclude_ligands, "no_ligands": a.no_ligands,
        "n_entities": len(system.entities),
        "entities_by_type": {k: len(v) for k, v in by_type.items()},
        "n_frames_written": int(samples.n_frames),
        "seconds_model_load": round(t_load, 1),
        "seconds_sampling": round(t_sample, 1),
        "note": "Protein-DNA COMPLEX ensemble. Does NOT feed stage2_redock; "
                "this track replaces stages 1-2. See bioemu2/README.md.",
    }
    (out / "sample_manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"[bioemu2/{label}] DONE -> {out}")


if __name__ == "__main__":
    main()
