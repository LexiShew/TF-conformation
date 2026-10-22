#!/usr/bin/env python
"""
warm_msa_cache.py - prefetch ColabFold MSAs for the BioEmu2 batch, SERIALLY.

Why this phase exists
---------------------
BioEmu2's `equilibrium` mode calls `get_msa(system)` before any diffusion, which
hits the PUBLIC ColabFold server (https://api.colabfold.com). That server
rate-limits per client. Running the sampling array at 14 concurrent tasks earned
20 immediate failures:

    requests.exceptions.HTTPError: 429 Client Error: Too Many Requests
        for url: https://api.colabfold.com/ticket/msa

Each failure landed ~31 s in, before the GPU did any work - so the array's real
concurrency ceiling is set by a public web API, not by GPUs. bioemu2's own retry
is stop_after_attempt(n) with wait_fixed(5), which is far too tight for a 429.

The fix is to stop fetching MSAs from inside GPU jobs at all. This script walks
the worklist on ONE cpu, fetches each MSA with real exponential backoff, and
writes them to the shared cache that BIOEMU2_MSA_CACHE points at. Once warm,
every sampling task is a cache hit, never contacts ColabFold, and the array's
concurrency is bounded only by available GPUs.

It is resumable and cheap to re-run: a cached system returns in well under a
second and is not slept on, so re-running after a partial pass costs minutes.

Dedup is free. The cache key is a sha256 over the system's UNIQUE protein
sequence set, so structures sharing a protein (5iy6-5iyd, 6cnb-6cnd, the 1a1*
zinc-finger series) collapse onto one entry and one request.

Usage
-----
    export BIOEMU2_MSA_CACHE=output/bioemu2/_msa_cache
    python bioemu2/warm_msa_cache.py \
        --worklist bioemu2/worklist/source_chains_survey.tsv \
        --report output/bioemu2/_msa_cache/warm_report.tsv

Exit code is 0 even with per-structure failures; read the report. A structure
that fails here will simply fetch its own MSA when its sampling task runs, so a
warm failure degrades throughput rather than correctness.
"""
import argparse
import csv
import sys
import time
from pathlib import Path

_HERE = Path(__file__).resolve().parent

# Backoff schedule for a 429, in seconds. Deliberately long: the limit is on a
# shared public service and being impatient is what caused the failures.
BACKOFF = [30, 60, 120, 300, 600]
# Pause after a genuine (uncached) fetch, to stay a polite client.
POLITE_SLEEP = 3.0
# A call faster than this was served from cache, so no pause is warranted.
CACHE_HIT_SECONDS = 1.0


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--worklist", default=str(_HERE / "worklist" / "source_chains_survey.tsv"),
                   help="TSV with pdb_id, cif and status columns (survey or worklist)")
    p.add_argument("--report", default=None)
    p.add_argument("--limit", type=int, default=None)
    a = p.parse_args()

    from bioemu2 import System
    from bioemu2.msa import get_msa
    from bioemu2.msa.colabfold import msa_cache_dir

    cache = msa_cache_dir()
    print(f"[warm] cache dir : {cache}", flush=True)
    print(f"[warm] worklist  : {a.worklist}", flush=True)

    rows = [r for r in csv.DictReader(open(a.worklist), delimiter="\t")
            if r.get("status", "ok") == "ok"]
    if a.limit:
        rows = rows[:a.limit]
    print(f"[warm] {len(rows)} structures to walk", flush=True)

    results = []
    n_hit = n_fetch = n_fail = 0
    for i, r in enumerate(rows, 1):
        pdb, cif = r["pdb_id"], r["cif"]
        status, detail, elapsed = "ok", "", 0.0
        for attempt, wait in enumerate([0] + BACKOFF):
            if wait:
                print(f"[warm] {pdb}: backing off {wait}s then retrying "
                      f"(attempt {attempt+1})", flush=True)
                time.sleep(wait)
            try:
                t0 = time.time()
                system = System.from_cif(cif)
                get_msa(system)
                elapsed = time.time() - t0
                status, detail = "ok", ""
                break
            except Exception as e:                                   # noqa: BLE001
                detail = f"{type(e).__name__}: {e}"[:300]
                # Only a rate-limit is worth waiting out; a malformed cif is not.
                if "429" not in detail and "Too Many Requests" not in detail:
                    status = "error"
                    break
                status = "rate_limited"
        else:
            status = "rate_limited_gave_up"

        if status == "ok":
            if elapsed < CACHE_HIT_SECONDS:
                n_hit += 1
            else:
                n_fetch += 1
                time.sleep(POLITE_SLEEP)
        else:
            n_fail += 1
            print(f"[warm] {pdb}: {status} - {detail}", file=sys.stderr, flush=True)

        results.append((pdb, status, round(elapsed, 2), detail))
        if i % 25 == 0 or i == len(rows):
            print(f"[warm] {i}/{len(rows)}  cached={n_hit} fetched={n_fetch} "
                  f"failed={n_fail}  entries={len(list(cache.glob('*.json')))}",
                  flush=True)

    if a.report:
        rp = Path(a.report)
        rp.parent.mkdir(parents=True, exist_ok=True)
        with rp.open("w") as fh:
            fh.write("pdb_id\tstatus\tseconds\tdetail\n")
            for row in results:
                fh.write("\t".join(str(x) for x in row) + "\n")
        print(f"[warm] report -> {rp}", flush=True)

    print(f"[warm] DONE cached={n_hit} fetched={n_fetch} failed={n_fail}", flush=True)


if __name__ == "__main__":
    main()
