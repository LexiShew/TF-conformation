#!/usr/bin/env python3
"""
extract_metrics.py — pipeline output → the metric CSVs the figure suite reads.

Auto-discovers every pilot from disk (fig_common.discover_pilots) and writes,
into analysis/data/:

  perstate_metrics.csv    per-state fnat / iRMSD / n_iface_res, Stage 2 + Stage 3
                          (from each pilot's fnat CSVs). Drives F* and I*.
  perentry_accuracy.csv   per-(pilot,entry) baseline vs augmented metrics over
                          the 130-entry benchmark (from stage7 eval JSONs,
                          seed-averaged). Drives P2/P3.
  mechanism_apo_holo.csv  per-pilot ensemble spread / reachability / aug ΔPearson
                          / curated dna_deform. Drives M1.

family_annotation.csv is NOT regenerated here: the 130-entry benchmark is fixed
and its motif-level entry→family map is pilot-independent. We validate it covers
every entry and otherwise leave it untouched (regenerate separately if the
benchmark set itself ever changes).

Fnat CSV schema (per pilot, output/stage{2,3}*/<pilot>/<pdb>_fnat.csv):
  pdb_id,state,n_iface_res,n_segments,seq_ident,iRMSD_global,
  iRMSD_seg_max,iRMSD_seg_mean,fnat
Eval JSON: results[<ckpt>][<entry.npz>] = {pearsonr,spearmanr,auroc,
  ic_weighted_pcc,ic_corr,brier_multi,ic_diff,mae}; ckpt = {baseline|augmented}_
  <pilot>_fold0[_sN].
"""
import os, glob, json, argparse
import numpy as np
import pandas as pd
from fig_common import (TFCONF, DATA_DIR, STAGE3_DIR, EVAL_DIR, discover_pilots,
                        parse_pilot_config, family_of, dna_deform_of)

STAGE2_DIR = os.path.join(TFCONF, "output", "stage2_docked")
STAGE3_RELAX = os.path.join(TFCONF, "output", "stage3_min_dnarelax")
FNAT_METRIC_COLS = ["n_iface_res", "n_segments", "seq_ident",
                    "iRMSD_global", "iRMSD_seg_max", "iRMSD_seg_mean", "fnat"]


def _read_fnat(stage_dir, pilot, pdb, stage_tag):
    """Read one pilot's fnat CSV under stage_dir, tag with pilot+stage."""
    # fnat CSVs are written next to the states; try <stage_dir>/<pilot>/<pdb>_fnat.csv
    for cand in (os.path.join(stage_dir, pilot, f"{pdb}_fnat.csv"),
                 os.path.join(stage_dir, pilot, f"{pdb.lower()}_fnat.csv")):
        if os.path.exists(cand):
            df = pd.read_csv(cand)
            df.insert(0, "pilot", pilot)
            df.insert(1, "stage", stage_tag)
            return df
    return None


def build_perstate(pilots):
    rows = []
    for tf in pilots:
        cfg = parse_pilot_config(tf)
        pdb = (cfg.get("PDB_ID") or "").lower()
        if not pdb:
            print(f"  [perstate] {tf}: no PDB_ID in config, skipping")
            continue
        # Stage 2 fnat (diagnostic; may be absent — Stage 2 carries every state)
        s2 = _read_fnat(STAGE2_DIR, tf, pdb, "stage2")
        s3 = _read_fnat(STAGE3_DIR, tf, pdb, "stage3")
        for df in (s2, s3):
            if df is not None:
                rows.append(df)
        got = [t for t, df in (("s2", s2), ("s3", s3)) if df is not None]
        print(f"  [perstate] {tf} ({pdb}): {'+'.join(got) if got else 'NO fnat CSV'}")
    if not rows:
        return pd.DataFrame()
    out = pd.concat(rows, ignore_index=True)
    keep = ["pilot", "stage", "pdb_id", "state"] + [c for c in FNAT_METRIC_COLS if c in out.columns]
    return out[keep]


def _seed_avg_entry_metrics(res, arm, pilot):
    """Average per-entry metrics across all seeds of one arm for one pilot.
    ckpt names: {arm}_{pilot}_fold0[_sN]. Returns {entry: {metric: mean}}."""
    import re
    per_entry = {}
    pat = re.compile(rf"^{arm}_{re.escape(pilot)}_fold0(_s\d+)?$")
    ckpts = [c for c in res if pat.match(c)]
    for ck in ckpts:
        for entry, m in res[ck].items():
            if not isinstance(m, dict):
                continue
            per_entry.setdefault(entry, []).append(m)
    out = {}
    for entry, ms in per_entry.items():
        keys = ms[0].keys()
        out[entry] = {k: float(np.nanmean([mm.get(k, np.nan) for mm in ms])) for k in keys}
    return out


def build_perentry(pilots):
    rows = []
    for tf in pilots:
        jf = os.path.join(EVAL_DIR, f"id_benchmark_{tf}.json")
        if not os.path.exists(jf):
            print(f"  [perentry] {tf}: no eval JSON, skipping")
            continue
        res = json.load(open(jf)).get("results", {})
        base = _seed_avg_entry_metrics(res, "baseline", tf)
        aug = _seed_avg_entry_metrics(res, "augmented", tf)
        common = sorted(set(base) & set(aug))
        for e in common:
            row = {"pilot": tf, "entry": e}
            for met in ("pearsonr", "spearmanr", "auroc", "ic_weighted_pcc", "mae"):
                row[f"base_{met}"] = base[e].get(met, np.nan)
                row[f"aug_{met}"] = aug[e].get(met, np.nan)
            rows.append(row)
        print(f"  [perentry] {tf}: {len(common)} entries (base∩aug), "
              f"{len(base)} base / {len(aug)} aug ckpt-entries")
    return pd.DataFrame(rows)


def build_mechanism(pilots, perentry, out_dir):
    """Assemble per-pilot mechanism axes for M1.

    Two axes need per-frame COORDINATE computation (ensemble spread vs the
    crystal, and reachability d_min / reach_ratio) that the fnat/eval CSVs do
    not contain. Those are carried forward from an existing
    mechanism_apo_holo.csv when present (so re-running never DESTROYS the
    reachability axis for pilots already computed), and left NaN for a new
    pilot until the coordinate pass is run for it. The columns this function
    can derive from tabular data — `spread` (BioEmu median pairwise), own-family
    `aug_dP`, and curated `dna_deform` — are always (re)computed.

    aug_dP definition (NOTE — differs from the original M1 CSV): here it is the
    mean (aug−base) Pearson over the pilot's OWN-family benchmark entries. The
    original mechanism_apo_holo.csv used a MIXED definition (a per-pilot
    `aug_kind` flag: 'held-out' subset where available, else 'gen-130'
    whole-benchmark), so the two do NOT match numerically — e.g. for foxa the
    original gen-130 value is +0.010 vs the own-family value here −0.006 (opposite
    sign), and lef1's original gen-130 −0.038 becomes NaN here because the HMG-box
    family has no *other* same-family benchmark entries to average. The own-family
    definition is the cleaner "does augmentation help this TF's own family"
    quantity; if you need to reproduce the original M1 numbers exactly, read
    aug_dP from the pre-existing CSV instead of recomputing. Pilots whose family
    has no same-family entries get NaN (make_M.py handles/annotates them).
    """
    div_path = os.path.join(DATA_DIR, "ensemble_diversity.csv")
    fam_path = os.path.join(DATA_DIR, "family_annotation.csv")
    div = pd.read_csv(div_path) if os.path.exists(div_path) else pd.DataFrame()
    fam = pd.read_csv(fam_path) if os.path.exists(fam_path) else pd.DataFrame()
    # Coordinate-derived reachability axes: prefer a fresh reachability.csv
    # (compute_reachability.py, all pilots) if present; else carry forward from
    # the existing mechanism CSV. reachability.csv WINS — it is the current
    # coordinate pass and covers every pilot it was run for.
    reach_path = os.path.join(DATA_DIR, "reachability.csv")
    prev_path = os.path.join(out_dir, "mechanism_apo_holo.csv")
    prev = pd.read_csv(prev_path) if os.path.exists(prev_path) else pd.DataFrame()
    prev_by = prev.set_index("pilot").to_dict("index") if len(prev) else {}
    if os.path.exists(reach_path):
        rdf = pd.read_csv(reach_path)
        for _, rr in rdf.iterrows():
            p = rr["pilot"]
            prev_by.setdefault(p, {})
            for col in ("d_min", "reach_ratio", "rmsf_mean"):
                if col in rr and not pd.isna(rr[col]):
                    prev_by[p][col] = rr[col]

    rows = []
    for tf in pilots:
        spread = np.nan
        if len(div):
            b = div[(div.pilot == tf) & (div.source == "bioemu")]
            if len(b):
                spread = float(b["median_pairwise"].iloc[0])
        aug_dP = np.nan
        if len(perentry) and len(fam):
            fam_of_entry = fam.set_index("entry")["family"].to_dict()
            sub = perentry[perentry.pilot == tf].copy()
            sub["family"] = sub["entry"].map(fam_of_entry)
            same = sub[sub["family"] == family_of(tf)]
            if len(same):
                aug_dP = float((same["aug_pearsonr"] - same["base_pearsonr"]).mean())
        prow = prev_by.get(tf, {})
        rows.append({
            "pilot": tf, "family": family_of(tf),
            "spread": spread,
            # carried forward from a prior coordinate pass; NaN until computed
            "d_min": prow.get("d_min", np.nan),
            "reach_ratio": prow.get("reach_ratio", np.nan),
            "rmsf_mean": prow.get("rmsf_mean", np.nan),
            "aug_dP": aug_dP,
            "dna_deform": dna_deform_of(tf),
            "reach_status": "computed" if not pd.isna(prow.get("d_min", np.nan)) else "needs_coord_pass",
        })
    return pd.DataFrame(rows)


def validate_family_annotation(perentry):
    fam_path = os.path.join(DATA_DIR, "family_annotation.csv")
    if not os.path.exists(fam_path):
        print("  [family] WARNING: family_annotation.csv absent — P2/P3 need it.")
        return
    fam = pd.read_csv(fam_path)
    known = set(fam.entry.unique())
    need = set(perentry.entry.unique()) if len(perentry) else set()
    missing = need - known
    if missing:
        print(f"  [family] WARNING: {len(missing)} benchmark entries lack a family "
              f"mapping (e.g. {sorted(missing)[:3]}). P2/P3 will drop them.")
    else:
        print(f"  [family] OK — all {len(need)} entries have a family mapping.")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pilots", nargs="+", default=None,
                    help="explicit pilot list; default = auto-discover from disk")
    ap.add_argument("--require", nargs="+", default=["stage3"],
                    help="evidence a pilot must have (stage3/eval/config)")
    ap.add_argument("--out", default=DATA_DIR, help="output dir for CSVs")
    args = ap.parse_args()

    pilots = args.pilots or discover_pilots(require=tuple(args.require))
    print(f"Pilots ({len(pilots)}): {', '.join(pilots)}")
    os.makedirs(args.out, exist_ok=True)

    print("Building perstate_metrics.csv …")
    ps = build_perstate(pilots)
    if len(ps):
        ps.to_csv(os.path.join(args.out, "perstate_metrics.csv"), index=False)
        print(f"  wrote {len(ps)} rows, {ps.pilot.nunique()} pilots")

    print("Building perentry_accuracy.csv …")
    pe = build_perentry([t for t in pilots
                         if os.path.exists(os.path.join(EVAL_DIR, f"id_benchmark_{t}.json"))])
    if len(pe):
        pe.to_csv(os.path.join(args.out, "perentry_accuracy.csv"), index=False)
        print(f"  wrote {len(pe)} rows, {pe.pilot.nunique()} pilots")
        validate_family_annotation(pe)

    print("Building mechanism_apo_holo.csv …")
    me = build_mechanism(pilots, pe, args.out)
    me.to_csv(os.path.join(args.out, "mechanism_apo_holo.csv"), index=False)
    print(f"  wrote {len(me)} pilots")
    print("Done.")


if __name__ == "__main__":
    main()
