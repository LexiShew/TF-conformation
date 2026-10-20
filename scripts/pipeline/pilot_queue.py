#!/usr/bin/env python3
"""
pilot_queue.py — keep the cluster busy with pilots without hand-babysitting.

The pipeline is 9 SLURM jobs per (pilot, arm) with dependency edges, submitted
by run_pilot.sh. With 12 pilots x 2 arms that is 24 chains, and deciding which
to launch next has been manual. This turns it into one idempotent command that
you can run from cron:

    # see what is done, running, and next — changes nothing
    python scripts/pipeline/pilot_queue.py --status

    # submit work until MAX_INFLIGHT chains are active, then stop
    python scripts/pipeline/pilot_queue.py --fill

    # what --fill would do, without submitting
    python scripts/pipeline/pilot_queue.py --fill --dry-run

    # keep it topped up forever (or from cron every 30 min, which is better)
    */30 * * * * cd /path/to/TF-conformation && python scripts/pipeline/pilot_queue.py --fill >> slurm_output/queue.log 2>&1

Design notes
------------
* State is inferred from the filesystem, not tracked in a database. A stage is
  "done" if its output exists and is non-empty. That makes the queue restartable
  and crash-proof: delete an output dir and it gets redone.
* Nothing is submitted for a (pilot, arm) that already has jobs in squeue —
  matched by the job-name prefix run_pilot.sh uses.
* Pilots that fail validation are skipped loudly, never submitted.
* --fill submits at most one chain per pilot/arm per invocation, so a cron loop
  ramps up smoothly instead of dumping 24 chains at once.

Priority order for what to run next:
  1. arm=frozen before arm=dnarelax     (frozen is the baseline everything pairs to)
  2. pilots closest to completion first (finish work in progress before starting new)
  3. registry order                      (stable tiebreak)
"""
import argparse
import getpass
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)

REGISTRY = os.path.join(ROOT, "config", "pilots.tsv")
OUT = os.path.join(ROOT, "output")

DEFAULT_MAX_INFLIGHT = 4
ARMS = ("frozen", "dnarelax")


def suffix_for(arm):
    return "" if arm == "frozen" else "_dnarelax"


def read_registry():
    import csv
    with open(REGISTRY) as fh:
        lines = [l for l in fh if not l.startswith("#") and l.strip()]
    return [{k: (v or "").strip() for k, v in r.items() if k}
            for r in csv.DictReader(lines, delimiter="\t")]


def nonempty_dir(p):
    return os.path.isdir(p) and any(os.scandir(p))


def count_files(p, pat=None):
    """Count matching entries, INCLUDING symlinks whose target is absent.

    The fnat gate builds STAGE3_PASS_DIR as a symlink mirror of STAGE3_DIR, and
    in a local checkout those links point at cluster paths that do not resolve.
    DirEntry.is_file() follows symlinks, so it returns False for them and the
    gate looked un-run. Count the link itself.
    """
    if not os.path.isdir(p):
        return 0
    n = 0
    for e in os.scandir(p):
        if pat is not None and not e.name.endswith(pat):
            continue
        if e.is_file() or e.is_symlink():
            n += 1
    return n


def stage_state(tf, pdb, arm):
    """Return an ordered list of (stage, done, detail)."""
    sfx = suffix_for(arm)
    s = []

    # Stage 1 is arm-independent (BioEmu library is shared)
    s1 = os.path.join(OUT, "stage1_bioemu")
    s1_hits = [d for d in (os.listdir(s1) if os.path.isdir(s1) else [])
               if d.startswith(f"{pdb}_chain") and d.endswith("_conformations")]
    s.append(("1 bioemu", bool(s1_hits), f"{len(s1_hits)} chain dir(s)"))

    # Stage 2 is also arm-independent (docking is DNA-restraint-independent)
    s2 = os.path.join(OUT, "stage2_docked", tf)
    n2 = count_files(s2, ".pdb")
    s.append(("2 redock", n2 > 0, f"{n2} docked pdb"))

    s3 = os.path.join(OUT, f"stage3_min{sfx}", tf)
    n3 = count_files(s3, ".pdb")
    s.append(("3 minimize", n3 > 0, f"{n3} minimized pdb"))

    s3p = os.path.join(OUT, f"stage3_min{sfx}", f"{tf}_pass")
    n3p = count_files(s3p, ".pdb")
    s.append(("3g fnat gate", n3p > 0, f"{n3p} passed gate"))

    s4 = os.path.join(OUT, f"stage4_npz{sfx}", tf, "output")
    n4 = count_files(s4, ".npz")
    gate_note = f"{n4} npz"
    if n3p and n4:
        lost = n3p - n4
        if lost > 0:
            gate_note += f"  !! {lost} lost ({100*lost/n3p:.0f}%) in stage 4"
    s.append(("4 preprocess", n4 > 0, gate_note))

    s5 = os.path.join(OUT, f"stage5_aug{sfx}", f"combined_assembly_{tf}")
    n5 = count_files(s5, ".npz")
    s.append(("5 build aug", n5 > 0, f"{n5} assembly npz"))

    s6 = os.path.join(OUT, "stage6_train")
    runs = [d for d in (os.listdir(s6) if os.path.isdir(s6) else [])
            if re.match(rf"^(baseline|augmented)_{tf}_fold\d+{sfx}_s\d+$", d)
            and os.path.exists(os.path.join(s6, d, "Model.best.tar"))]
    s.append(("6 train", len(runs) >= 2, f"{len(runs)} checkpoint(s)"))

    ev = os.path.join(OUT, "stage7_eval", f"id_benchmark_{tf}{sfx}.json")
    s.append(("7 eval", os.path.exists(ev) and os.path.getsize(ev) > 0,
              "json present" if os.path.exists(ev) else "missing"))
    return s


def first_incomplete(states):
    for i, (name, done, _) in enumerate(states):
        if not done:
            return i, name
    return None, None


def squeue_names():
    user = os.environ.get("SLURM_JOB_USER") or getpass.getuser()
    try:
        r = subprocess.run(["squeue", "-u", user, "-h", "-o", "%j|%T"],
                           capture_output=True, text=True, timeout=30)
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return None  # no SLURM here
    out = []
    for line in r.stdout.splitlines():
        if "|" in line:
            name, state = line.rsplit("|", 1)
            out.append((name.strip(), state.strip()))
    return out


def inflight_for(jobs, tf, arm):
    """Does this (pilot, arm) already have queued/running jobs?"""
    if jobs is None:
        return False
    sfx = suffix_for(arm)
    for name, _ in jobs:
        # run_pilot.sh job names carry the stage script name; the pilot comes in
        # via --export TF_NAME, so match on the TF token when present.
        if tf in name and (sfx in name or (sfx == "" and "_dnarelax" not in name)):
            return True
    return False


def validate(tf):
    r = subprocess.run([sys.executable,
                        os.path.join(HERE, "validate_pilot.py"), tf, "--quiet"],
                       capture_output=True, text=True)
    return r.returncode == 0, r.stdout + r.stderr


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--status", action="store_true", help="show state, change nothing")
    ap.add_argument("--fill", action="store_true", help="submit until the queue is full")
    ap.add_argument("--dry-run", action="store_true", help="with --fill, print only")
    ap.add_argument("--max-inflight", type=int,
                    default=int(os.environ.get("MAX_INFLIGHT", DEFAULT_MAX_INFLIGHT)),
                    help=f"max concurrent pilot chains (default {DEFAULT_MAX_INFLIGHT})")
    ap.add_argument("--only", help="restrict to this pilot")
    ap.add_argument("--arm", choices=ARMS, help="restrict to this arm")
    a = ap.parse_args()
    if not (a.status or a.fill):
        a.status = True

    rows = read_registry()
    if a.only:
        rows = [r for r in rows if r["tf_name"] == a.only]
    jobs = squeue_names()

    units = []
    for order, row in enumerate(rows):
        tf, pdb = row["tf_name"], row["pdb_id"].lower()
        for arm in (ARMS if not a.arm else (a.arm,)):
            st = stage_state(tf, pdb, arm)
            idx, nm = first_incomplete(st)
            units.append(dict(tf=tf, pdb=pdb, arm=arm, states=st, order=order,
                              next_idx=idx, next_name=nm,
                              done=idx is None,
                              inflight=inflight_for(jobs, tf, arm)))

    # ---------------- status ----------------
    if a.status or a.fill:
        print(f"{'PILOT':<11}{'ARM':<10}{'PROGRESS':<22}{'NEXT STAGE':<16}STATE")
        print("-" * 78)
        for u in units:
            total = len(u["states"])
            ndone = sum(1 for _, d, _ in u["states"] if d)
            progress = "#" * ndone + "." * (total - ndone) + f"  {ndone}/{total}"
            state = ("COMPLETE" if u["done"]
                     else "in flight" if u["inflight"] else "idle")
            print(f"{u['tf']:<11}{u['arm']:<10}{progress:<22}"
                  f"{(u['next_name'] or '-'):<16}{state}")
        n_c = sum(1 for u in units if u["done"])
        n_f = sum(1 for u in units if u["inflight"])
        print("-" * 78)
        print(f"{len(units)} unit(s): {n_c} complete, {n_f} in flight, "
              f"{len(units)-n_c-n_f} idle   (max-inflight={a.max_inflight})")
        if jobs is None:
            print("NOTE: squeue unavailable here — in-flight detection is off, so "
                  "--fill would re-submit. Run this on the cluster.")

        # Surface stage-4 attrition, the silent failure that cost the IRF pilot
        warn = [(u, d) for u in units for (n, _, d) in u["states"]
                if n == "4 preprocess" and "!!" in d]
        if warn:
            print("\nStage-4 attrition detected (see docs/retired_pilots.md — this is "
                  "what made IRF unusable):")
            for u, d in warn:
                print(f"  {u['tf']:<11}{u['arm']:<10}{d}")

    if not a.fill:
        return 0

    # ---------------- fill ----------------
    if jobs is None and not a.dry_run:
        print("\nERROR: squeue not available; refusing to submit blind.", file=sys.stderr)
        return 1

    slots = a.max_inflight - sum(1 for u in units if u["inflight"])
    if slots <= 0:
        print(f"\nQueue full ({a.max_inflight} in flight). Nothing to do.")
        return 0

    # frozen first, then furthest-along, then registry order
    cand = [u for u in units if not u["done"] and not u["inflight"]]
    cand.sort(key=lambda u: (0 if u["arm"] == "frozen" else 1,
                             -u["next_idx"], u["order"]))

    print(f"\n{slots} slot(s) free; {len(cand)} candidate(s)")
    launched = 0
    for u in cand:
        if launched >= slots:
            break
        ok, msg = validate(u["tf"])
        if not ok:
            print(f"  SKIP {u['tf']}/{u['arm']}: failed validation")
            for line in msg.strip().splitlines():
                if "ERROR" in line:
                    print(f"       {line.strip()}")
            continue
        stage_start = u["next_idx"] + 1
        # the "3g fnat gate" and "3 minimize" slots both map to stage 3/4 args
        stage_arg = {1: "1", 2: "2", 3: "3", 4: "4", 5: "4", 6: "5", 7: "6", 8: "7"}
        start = stage_arg.get(stage_start, "1")
        cmd = ["./scripts/pipeline/run_pilot.sh", u["tf"], start, "7"]
        env = dict(os.environ, ARM=u["arm"])
        print(f"  {'DRY-RUN  ' if a.dry_run else 'SUBMIT   '}"
              f"ARM={u['arm']} {' '.join(cmd)}   (resuming at '{u['next_name']}')")
        if not a.dry_run:
            r = subprocess.run(cmd, cwd=ROOT, env=env,
                               capture_output=True, text=True)
            if r.returncode != 0:
                print(f"       FAILED: {r.stderr.strip().splitlines()[-1] if r.stderr.strip() else 'unknown'}")
                continue
            for line in r.stdout.splitlines():
                if "submitted" in line:
                    print(f"       {line.strip()}")
        launched += 1

    print(f"\n{launched} chain(s) {'would be ' if a.dry_run else ''}launched.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
