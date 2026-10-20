#!/usr/bin/env python3
"""
validate_pilot.py — pre-flight correctness checks for a pilot config.

This exists because of the TBP label error (docs/tbp_label_error_downstream.md):
config/pilots/tbp.sh labelled a TBP structure with NDT80's motif, the mislabel
was baked into every npz's Y_pwm, and it survived into the analysis layer where
a hand-curation made it self-consistent and invisible. Check C1 below would
have caught it on day one.

Run it standalone over everything:
    python scripts/pipeline/validate_pilot.py --all
    python scripts/pipeline/validate_pilot.py --all --strict   # warnings fail too

Or for one pilot (this is what run_pilot.sh does before submitting):
    python scripts/pipeline/validate_pilot.py tbp

Exit codes: 0 = clean (or warnings only), 1 = at least one ERROR.

Checks
------
C1  PWM_LABEL must be listed for PDB_ID in config/deeppbs_pdb_pwms.csv.
C2  TEST_PWM_FILTER must compile and match >=1 id.txt entry.
C3  Matched on-target entries must agree with the pilot's own family.
C4  ON_TARGET_USABLE / ON_TARGET_N_EXPECTED must match reality.
C5  The pilot's source PDB must not be in id.txt (test-set leakage).
C6  BINDING_CHAIN must name a real source-chain file.
C7  PROTEIN_CHAIN / DNA_CHAINS must be disjoint and plausibly typed.
C8  Required vars present; N_FRAMES sane.
C9  Annotation invariant: where an entry's own motif and its Pfam agree on a
    family, family_annotation.csv must agree too.
C10 Training-config seeding: every stage6 training config must carry a usable
    random_seed, and driver.py must actually read it. This is finding F6 --
    all 264 historical runs were unseeded because the key was written but
    never read, and the seeding block was gated on a flag that was always
    False. Both halves are checked so neither can regress alone.
"""
import argparse
import csv
import glob
import os
import re
import subprocess
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
PILOTS = os.path.join(ROOT, "config", "pilots")
PWM_CSV = os.path.join(ROOT, "config", "deeppbs_pdb_pwms.csv")
ID_TXT = os.path.join(ROOT, "data", "folds", "id.txt")
FAM_CSV = os.path.join(ROOT, "analysis", "data", "family_annotation.csv")
SRC_CHAINS = os.path.join(ROOT, "structures", "source_chains")

# Pfam -> canonical family stem, for the C9 invariant. Only codes we can assert
# on; anything absent here is skipped rather than guessed at.
PFAM_FAMILY = {
    "PF00352": "tbp", "PF05224": "ndt80", "PF00178": "ets", "PF00853": "runt",
    "PF00046": "homeodomain", "PF00096": "c2h2", "PF00105": "nuclear",
    "PF07716": "bzip", "PF00170": "bzip", "PF00172": "zn2cys6",
    "PF00125": "histone", "PF00319": "mads", "PF01342": "bhlh",
    "PF00010": "bhlh", "PF00554": "rel", "PF00870": "p53", "PF00447": "hsf",
    "PF00505": "hmg", "PF02257": "irf", "PF09271": "csl", "PF09270": "csl",
    "PF00157": "pou", "PF00292": "pax",
}

# Entries whose pfam legitimately disagrees with family because the motif-level
# convention applies: the scored PWM belongs to a different TF than the
# crystal's first chain (co-crystals). Documented in analysis/README.md.
# Keyed by entry prefix; value is a short reason.
MOTIF_LEVEL_EXEMPT = {
    "4l0z_B": "ETS1 motif on a RUNX1-first co-crystal",
    "3wtu_C": "ETS1 motif on a RUNX1-first co-crystal",
    "3wtw_H": "ETS1 motif on a RUNX1-first co-crystal",
    "3wtw_C": "ETS1 motif on a RUNX1-first co-crystal",
    "3wts_C": "ETS1 motif on a RUNX1-first co-crystal",
    "1hjb_F": "RUNX1 motif on a CEBPB chain",
    "7oh9_K": "TBP motif on a histone chain (nucleosome/TFIIA/TBP assembly)",
}


class Result:
    def __init__(self):
        self.errors = []
        self.warnings = []
        self.passed = []

    def err(self, code, msg):
        self.errors.append(f"[{code}] {msg}")

    def warn(self, code, msg):
        self.warnings.append(f"[{code}] {msg}")

    def ok(self, code, msg):
        self.passed.append(f"[{code}] {msg}")


def load_config(tf_name):
    """Source the pilot config in bash and return its exported vars.

    Sourcing (rather than regex-parsing) is deliberate: configs like
    tbp_dnarelax.sh `source` a parent, and several set values via shell
    expansion. Parsing the text would miss both.
    """
    cfg = os.path.join(PILOTS, f"{tf_name}.sh")
    if not os.path.exists(cfg):
        return None
    script = f'set -a; source "{cfg}" >/dev/null 2>&1; env'
    try:
        out = subprocess.run(["bash", "-c", script], capture_output=True,
                             text=True, timeout=30, cwd=ROOT)
    except subprocess.TimeoutExpired:
        return None
    env = {}
    for line in out.stdout.splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            env[k] = v
    return env


def load_pwm_table():
    table = {}
    if not os.path.exists(PWM_CSV):
        return table
    with open(PWM_CSV) as fh:
        for row in csv.DictReader(fh):
            table[row["pdb_id"].strip().lower()] = row
    return table


def load_ids():
    if not os.path.exists(ID_TXT):
        return []
    with open(ID_TXT) as fh:
        return [l.strip() for l in fh if l.strip()]


def load_families():
    """entry -> family, plus pdb -> set(families), from family_annotation.csv."""
    by_entry, by_pdb = {}, defaultdict(set)
    if not os.path.exists(FAM_CSV):
        return by_entry, by_pdb, []
    rows = []
    with open(FAM_CSV) as fh:
        for row in csv.DictReader(fh):
            rows.append(row)
            e = row["entry"].replace(".npz", "")
            by_entry[e] = row["family"]
            by_pdb[row["pdb"].lower()].add(row["family"])
    return by_entry, by_pdb, rows


def norm(s):
    return re.sub(r"[^a-z0-9]", "", (s or "").lower())


def validate(tf_name, pwm_table, ids, fam_by_entry, verbose=False):
    r = Result()
    env = load_config(tf_name)
    if env is None:
        r.err("C8", f"could not source config/pilots/{tf_name}.sh")
        return r

    pdb = (env.get("PDB_ID") or "").lower()
    label = env.get("PWM_LABEL") or ""
    filt = env.get("TEST_PWM_FILTER") or ""
    fname = env.get("TEST_FILTER_NAME") or ""

    # ---------- C8: required vars ----------
    # NB: FOLD / RAMP_STAGES / STEPS_PER_STAGE / RECOVERY_* are NOT required
    # here. Since the config generator landed they live as lib/common.sh
    # defaults rather than being repeated in all 12 configs, so a generated
    # config legitimately omits them.
    missing = [v for v in ("TF_NAME", "PDB_ID", "PWM_LABEL", "BINDING_CHAIN",
                           "PROTEIN_CHAIN", "DNA_CHAINS")
               if not env.get(v)]
    if missing:
        r.err("C8", f"missing required var(s): {', '.join(missing)}")
    else:
        r.ok("C8", "all required vars present")

    nf = env.get("N_FRAMES")
    if nf and nf.isdigit():
        if not (1 <= int(nf) <= 1000):
            r.warn("C8", f"N_FRAMES={nf} outside sane range 1-1000")
    elif nf:
        r.err("C8", f"N_FRAMES={nf!r} is not an integer")

    # ---------- C1: PWM_LABEL vs master table ----------
    # This is the check that would have caught the TBP error.
    row = pwm_table.get(pdb)
    if row is None:
        r.warn("C1", f"PDB {pdb!r} not in deeppbs_pdb_pwms.csv — cannot verify "
                     f"PWM_LABEL={label!r}")
    else:
        allowed = [p for p in row.get("all_pwms", "").split(";") if p]
        if not allowed:
            r.warn("C1", f"deeppbs_pdb_pwms.csv lists no PWMs for {pdb}")
        elif label in allowed:
            r.ok("C1", f"PWM_LABEL={label} is listed for {pdb}")
        else:
            r.err("C1", f"PWM_LABEL={label!r} is NOT a PWM for {pdb}. "
                        f"Master table lists: {', '.join(allowed)}. "
                        f"This is the TBP-class error — the label does not "
                        f"belong to this structure.")

    # ---------- C5: test-set leakage ----------
    leaked = [e for e in ids if e.lower().startswith(pdb + "_")]
    if leaked:
        r.err("C5", f"source PDB {pdb} appears in id.txt as "
                    f"{', '.join(leaked[:3])}{'...' if len(leaked) > 3 else ''} — "
                    f"using it as an augmentation source leaks held-out data "
                    f"into training")
    else:
        r.ok("C5", f"source PDB {pdb} is not in id.txt (no leakage)")

    # ---------- C2/C3/C4: on-target subset ----------
    usable = (env.get("ON_TARGET_USABLE") or "").strip().lower()
    expected = env.get("ON_TARGET_N_EXPECTED")

    if not filt:
        r.warn("C2", "no TEST_PWM_FILTER set — no on-target subset")
        matched = []
    else:
        try:
            pat = re.compile(filt)
        except re.error as e:
            r.err("C2", f"TEST_PWM_FILTER={filt!r} is not a valid regex: {e}")
            return r
        matched = sorted(e for e in ids if pat.search(e))

        if not matched:
            if usable == "no":
                r.ok("C2", f"filter matches 0 entries, and ON_TARGET_USABLE=no "
                           f"declares that")
            else:
                r.err("C2", f"TEST_PWM_FILTER={filt!r} matches 0 of {len(ids)} "
                            f"id.txt entries. The on-target subset is empty, so "
                            f"stage 7 will silently skip it. Either fix the "
                            f"filter or set ON_TARGET_USABLE=no.")
        else:
            npdb = len({e.split("_")[0] for e in matched})
            r.ok("C2", f"filter matches {len(matched)} entries across {npdb} PDB(s)")

            # Pseudoreplication: entries are per-chain, so N chains of one
            # crystal are not N independent observations.
            if npdb < len(matched):
                sev = r.err if npdb == 1 and len(matched) >= 3 else r.warn
                sev("C2", f"on-target subset has {len(matched)} entries but only "
                          f"{npdb} unique PDB(s) — chains of the same crystal are "
                          f"not independent. Effective structural n is {npdb}, "
                          f"not {len(matched)}; paired tests over "
                          f"{len(matched)} will understate p.")

            if len(matched) < 3:
                r.warn("C2", f"n={len(matched)} is too small for a meaningful "
                             f"paired statistic")

            # ---------- C3: family agreement ----------
            fams = {fam_by_entry.get(e.replace(".npz", "")) for e in matched}
            fams.discard(None)
            own = fam_by_entry.get(None)
            # The pilot's own family: look it up from any entry sharing its label
            own_fams = {fam_by_entry[e] for e in fam_by_entry
                        if norm(label) and norm(label) in norm(e)}
            if len(fams) > 1:
                r.warn("C3", f"on-target subset spans {len(fams)} families "
                             f"({', '.join(sorted(fams))}) — this is a "
                             f"family-transfer test, not an on-target one. "
                             f"TEST_FILTER_NAME={fname!r} should say so.")
            elif fams and own_fams and not (fams & own_fams):
                r.err("C3", f"on-target subset family {fams} does not match the "
                            f"pilot's own motif family {own_fams} — the pilot is "
                            f"being scored on a different family (TBP-class error)")
            elif fams:
                r.ok("C3", f"on-target subset is one family: {fams.pop()}")

    # ---------- C4: declared vs actual ----------
    if expected and expected.isdigit():
        if len(matched) != int(expected):
            r.err("C4", f"ON_TARGET_N_EXPECTED={expected} but the filter matches "
                        f"{len(matched)}. Update the config or the filter — a "
                        f"silent drift here is how subsets change unnoticed.")
        else:
            r.ok("C4", f"on-target n={len(matched)} matches declared "
                       f"ON_TARGET_N_EXPECTED")
    if usable and usable not in ("yes", "no"):
        r.err("C4", f"ON_TARGET_USABLE={usable!r} must be 'yes' or 'no'")

    # ---------- C6: BINDING_CHAIN resolves to a real file ----------
    bc = env.get("BINDING_CHAIN")
    if bc:
        want = os.path.join(SRC_CHAINS, f"{pdb}_chains",
                            f"{pdb}_chain{bc}_protein.pdb")
        if os.path.isdir(os.path.join(SRC_CHAINS, f"{pdb}_chains")):
            if os.path.exists(want):
                r.ok("C6", f"BINDING_CHAIN={bc} -> {os.path.basename(want)} exists")
            else:
                avail = sorted(os.path.basename(p) for p in glob.glob(
                    os.path.join(SRC_CHAINS, f"{pdb}_chains",
                                 f"{pdb}_chain*_protein.pdb")))
                r.err("C6", f"BINDING_CHAIN={bc} but {os.path.basename(want)} "
                            f"does not exist. Available: "
                            f"{', '.join(avail) if avail else '(none)'}")
        else:
            r.warn("C6", f"source-chain dir for {pdb} not present locally "
                         f"(expected on the cluster) — cannot verify "
                         f"BINDING_CHAIN={bc}")

    # ---------- C7: chain indices disjoint ----------
    pc, dc = env.get("PROTEIN_CHAIN"), env.get("DNA_CHAINS")
    if pc is not None and dc:
        # PROTEIN_CHAIN may carry a trailing comment in some configs
        pc_clean = pc.split("#")[0].strip()
        try:
            pci = int(pc_clean)
            dci = [int(x) for x in dc.split(",") if x.strip()]
            if pci in dci:
                r.err("C7", f"PROTEIN_CHAIN={pci} also appears in "
                            f"DNA_CHAINS={dc} — a chain cannot be both")
            elif len(dci) != len(set(dci)):
                r.err("C7", f"DNA_CHAINS={dc} has duplicate entries")
            else:
                r.ok("C7", f"PROTEIN_CHAIN={pci} / DNA_CHAINS={dci} disjoint")
            if pc_clean != pc.strip():
                r.warn("C7", f"PROTEIN_CHAIN carries an inline comment "
                             f"({pc.strip()!r}) — move it to its own line so the "
                             f"value is unambiguous")
        except ValueError:
            r.err("C7", f"PROTEIN_CHAIN={pc!r} / DNA_CHAINS={dc!r} not integers")

    return r


def validate_annotation():
    """C9 — the annotation invariant, repo-wide (not per pilot).

    Where an entry's own motif and its Pfam agree on a family, the family
    column must agree too. Co-crystal cases where the scored motif belongs to a
    different TF than the crystal's first chain are exempt (MOTIF_LEVEL_EXEMPT).
    """
    r = Result()
    _, _, rows = load_families()
    if not rows:
        r.warn("C9", "analysis/data/family_annotation.csv not found — skipped")
        return r

    seen, bad = set(), []
    for row in rows:
        entry = row["entry"].replace(".npz", "")
        if entry in seen:
            continue
        seen.add(entry)
        pf, fam = row.get("pfam", ""), row.get("family", "")
        stem = PFAM_FAMILY.get(pf)
        if not stem or not fam:
            continue
        key = "_".join(entry.split("_")[:2])
        if key in MOTIF_LEVEL_EXEMPT:
            continue
        if stem not in norm(fam) and norm(fam)[:4] not in stem:
            bad.append((entry, pf, stem, fam, row.get("family_src", "")))

    if bad:
        for entry, pf, stem, fam, src in bad:
            r.err("C9", f"{entry}: pfam {pf} means '{stem}' but family='{fam}' "
                        f"(family_src={src}). Either the family label is wrong "
                        f"or this entry needs a MOTIF_LEVEL_EXEMPT reason.")
    else:
        r.ok("C9", f"annotation invariant holds across {len(seen)} entries "
                   f"({len(MOTIF_LEVEL_EXEMPT)} documented motif-level exemptions)")
    return r


def validate_seeding():
    """C10 - training configs are seeded AND driver.py reads the seed.

    Checks both halves of finding F6:
      (a) driver.py must reference C["random_seed"], not just hardcoded values.
          Without this, the key below is decorative.
      (b) every output/stage6_train/config_*.json must carry a non-null
          random_seed, or declare no_random/debug true.
    """
    r = Result()
    driver = os.path.join(ROOT, "stage6_train", "driver.py")
    if not os.path.exists(driver):
        r.warn("C10", "stage6_train/driver.py not found - skipped")
        return r

    raw = open(driver).read()
    # Strip full-line comments before scanning. driver.py documents the old
    # broken block verbatim in a comment, so a naive text search finds the
    # legacy literals in prose and reports a shadowing risk that is not real.
    src = "\n".join(l for l in raw.splitlines()
                    if not l.lstrip().startswith("#"))
    reads_seed = 'C.get("random_seed"' in src or 'C["random_seed"]' in src
    if not reads_seed:
        r.err("C10", "stage6_train/driver.py does not read random_seed from the "
                     "config. Any random_seed written by stage5 is decorative "
                     "and every run is unseeded. This is exactly finding F6 - "
                     "do not let it regress.")
    else:
        r.ok("C10", "driver.py reads random_seed from the config")

    # Legacy trap: the old block seeded from literals regardless of the config.
    # It is legitimate for the 8/0 literals to survive INSIDE the no_random
    # fallback branch, so only warn when they appear before the random_seed
    # branch -- i.e. when they could still shadow it.
    if "np.random.seed(8)" in src:
        i_lit = src.index("np.random.seed(8)")
        i_seed = min([src.index(k) for k in ('C.get("random_seed"', 'C["random_seed"]')
                      if k in src] or [len(src)])
        if i_lit < i_seed:
            r.warn("C10", "driver.py's legacy np.random.seed(8) literal appears "
                          "BEFORE the random_seed branch and may shadow it")
        else:
            r.ok("C10", "legacy 8/0 literals are confined to the no_random "
                        "fallback, after the random_seed branch")

    cfg_dir = os.path.join(ROOT, "output", "stage6_train")
    cfgs = sorted(glob.glob(os.path.join(cfg_dir, "config_*.json")))
    if not cfgs:
        r.warn("C10", "no output/stage6_train/config_*.json found - nothing to "
                      "check (expected on a fresh checkout)")
        return r

    import json as _json
    unseeded = []
    for c in cfgs:
        try:
            d = _json.load(open(c))
        except Exception:
            continue
        seed = d.get("random_seed", None)
        if seed is None and not (d.get("no_random") or d.get("debug")):
            unseeded.append(os.path.basename(c))
    if unseeded:
        r.err("C10", f"{len(unseeded)} of {len(cfgs)} training config(s) have no "
                     f"random_seed and no no_random/debug flag, so their runs are "
                     f"unseeded and irreproducible: "
                     f"{', '.join(unseeded[:4])}"
                     f"{'...' if len(unseeded) > 4 else ''}")
    else:
        r.ok("C10", f"all {len(cfgs)} training config(s) carry a random_seed")

    # Paired-seed sanity: baseline_X_sN and augmented_X_sN must share a seed.
    import re as _re
    by_pair = defaultdict(dict)
    for c in cfgs:
        b = os.path.basename(c)
        m = _re.match(r"config_(baseline|augmented(?:_legacy)?)_(.+)\.json$", b)
        if not m:
            continue
        kind, suffix = m.groups()
        try:
            by_pair[suffix][kind] = _json.load(open(c)).get("random_seed")
        except Exception:
            continue
    mism = [(s, d) for s, d in by_pair.items()
            if "baseline" in d and "augmented" in d
            and d["baseline"] != d["augmented"]]
    if mism:
        for s, d in mism[:5]:
            r.err("C10", f"{s}: baseline seed={d['baseline']} but augmented "
                         f"seed={d['augmented']}. A paired comparison requires "
                         f"the SAME seed in both arms.")
    elif by_pair:
        n_ok = sum(1 for d in by_pair.values()
                   if "baseline" in d and "augmented" in d)
        r.ok("C10", f"{n_ok} baseline/augmented pair(s) share a seed")
    return r


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pilot", nargs="?", help="pilot name (omit with --all)")
    ap.add_argument("--all", action="store_true", help="validate every pilot config")
    ap.add_argument("--strict", action="store_true", help="treat warnings as failures")
    ap.add_argument("--quiet", action="store_true", help="only show problems")
    ap.add_argument("--skip-annotation", action="store_true",
                    help="skip the repo-wide C9 annotation invariant")
    ap.add_argument("--skip-seeding", action="store_true",
                    help="skip the repo-wide C10 training-config seeding check")
    a = ap.parse_args()

    if not a.all and not a.pilot:
        ap.error("give a pilot name or --all")

    pwm_table = load_pwm_table()
    ids = load_ids()
    fam_by_entry, _, _ = load_families()

    if a.all:
        names = sorted(os.path.basename(p)[:-3]
                       for p in glob.glob(os.path.join(PILOTS, "*.sh")))
    else:
        names = [a.pilot]

    n_err = n_warn = 0
    for name in names:
        r = validate(name, pwm_table, ids, fam_by_entry)
        n_err += len(r.errors)
        n_warn += len(r.warnings)
        if r.errors or r.warnings or not a.quiet:
            status = "FAIL" if r.errors else ("WARN" if r.warnings else "ok")
            print(f"\n=== {name} [{status}] ===")
            for e in r.errors:
                print(f"  ERROR   {e}")
            for w in r.warnings:
                print(f"  warning {w}")
            if not a.quiet:
                for p in r.passed:
                    print(f"  pass    {p}")

    repo_checks = []
    if not a.skip_annotation:
        repo_checks.append(("annotation invariant", validate_annotation()))
    if not a.skip_seeding:
        repo_checks.append(("training-config seeding", validate_seeding()))
    if True:
        for title, r in repo_checks:
            n_err += len(r.errors)
            n_warn += len(r.warnings)
            status = "FAIL" if r.errors else ("WARN" if r.warnings else "ok")
            print(f"\n=== {title} (repo-wide) [{status}] ===")
            for e in r.errors:
                print(f"  ERROR   {e}")
            for w in r.warnings:
                print(f"  warning {w}")
            for pp in r.passed:
                print(f"  pass    {pp}")

    print(f"\n{'='*66}")
    print(f"{len(names)} pilot(s) checked — {n_err} error(s), {n_warn} warning(s)")
    if n_err:
        print("FAILED. Fix the errors above before running the pipeline.")
        return 1
    if n_warn and a.strict:
        print("FAILED under --strict (warnings present).")
        return 1
    print("PASSED.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
