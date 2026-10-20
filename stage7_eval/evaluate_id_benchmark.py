"""
Generic id.txt benchmark evaluation with statistical robustness improvements.

Improvements over v1:
  - Hard-checks that all conditions evaluated the same set of entries (n match)
  - Computes per-entry paired Δ-metrics across seeds and reports paired t-test
    p-values and bootstrap 95% CIs (since DeepPBS uses paired random seeds)
  - Explicit warnings when --condition names don't pair (baseline_X/augmented_X)

Each --condition is "name=path"; the name is shown in reports and used for
seed-pair matching.
"""
import argparse
import json
import os
import pickle
import re
import sys

import numpy as np
import torch
from torch_geometric.data import DataLoader

# Self-locate the repo so imports + inputs come from TF-conformation, not the
# DeepPBS trees: deeppbs from lib/, models.model_v2 from stage6_train/.
_TFCONF = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_TFCONF, "stage6_train"))  # models.model_v2
sys.path.insert(0, os.path.join(_TFCONF, "lib"))           # deeppbs package
# The project's one correct implementation of paired stats + BH-FDR. This module
# exists because this file and mechanism_analysis.py used to carry two separate
# implementations with two different bugs (findings F9 and F12).
sys.path.insert(0, os.path.join(_TFCONF, "analysis", "common"))
from stats import paired_delta, bh_fdr, effective_n

from deeppbs.nn.utils import loadDataset
from deeppbs.nn import Evaluator
from models.model_v2 import Model


NC = 4
LABELS_KEY = "Y_pwm"
DEFAULT_ID_FILE = os.path.join(_TFCONF, "data", "folds", "id.txt")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--condition", action="append", required=True,
                   help="name=run_dir; can be specified multiple times")
    p.add_argument("--combined-dir", required=True)
    p.add_argument("--id-file", default=DEFAULT_ID_FILE)
    p.add_argument("--filter", required=True,
                   help="Regex matching entries to focus on in subset report")
    p.add_argument("--filter-name", required=True)
    p.add_argument("--output-json", required=True)
    p.add_argument("--bootstrap-iters", type=int, default=10000,
                   help="Number of bootstrap samples for CI estimation")
    return p.parse_args()


def load_conditions(condition_args):
    conditions = {}
    for c in condition_args:
        if "=" not in c:
            print(f"ERROR: --condition must be name=path, got {c!r}", file=sys.stderr)
            sys.exit(1)
        name, path = c.split("=", 1)
        if not os.path.exists(os.path.join(path, "Model.best.tar")):
            print(f"WARNING: Model.best.tar not found in {path}; condition '{name}' skipped",
                  file=sys.stderr)
        conditions[name] = path
    return conditions


def evaluate_condition(cond_name, cond_dir, entries, combined_dir, device):
    if not os.path.exists(os.path.join(cond_dir, "Model.best.tar")):
        return None

    with open(os.path.join(cond_dir, "scaler.pkl"), "rb") as f:
        scaler = pickle.load(f)
    with open(os.path.join(cond_dir, "config.json")) as f:
        train_C = json.load(f)
    condition = train_C.get("condition", "prot_shape")
    readout = train_C.get("readout", "all")

    dataset, _, _, found_files = loadDataset(
        list(entries), NC, LABELS_KEY, combined_dir,
        cache_dataset=False, balance="unmasked", remove_mask=False,
        scale=True, scaler=scaler, pre_transform=None, feature_mask=None,
    )
    DL = DataLoader(dataset, batch_size=1, shuffle=False, pin_memory=True)

    model = Model(13, 14, condition=condition, readout=readout)
    state = torch.load(os.path.join(cond_dir, "Model.best.tar"), map_location=device)
    model.load_state_dict(state["model_state_dict"])
    model.to(device)

    evaluator = Evaluator(model, NC, device=device, post_process=torch.nn.Softmax(dim=-1))
    val_out = evaluator.eval(
        DL, use_mask=False, batchwise=True,
        return_masks=True, return_predicted=False, return_batches=True,
        xtras=None, threshold=0.5, eval_mode=True,
    )

    results = {}
    for i in range(val_out["num_batches"]):
        name = found_files[i]
        metrics = evaluator.getMetrics(
            val_out['y'][i], val_out['output'][i],
            val_out['masks'][i], val_out['out_masks'][i],
            val_out['batches'][i], threshold=0.5,
        )
        results[name] = metrics
    return results


def check_denominators(results, entries):
    """Verify each condition produced metrics for the same set of entries.
    Returns the common set; warns about per-condition missing entries.
    """
    if not results:
        return set()
    per_cond_keys = {cond: set(r.keys()) for cond, r in results.items()}
    intersection = set.intersection(*per_cond_keys.values())
    union = set.union(*per_cond_keys.values())
    if intersection != union:
        print("\n⚠ WARNING: conditions evaluated different entry sets")
        for cond, keys in per_cond_keys.items():
            missing = union - keys
            if missing:
                print(f"  {cond} is missing {len(missing)} entries: "
                      f"{sorted(missing)[:5]}{'...' if len(missing)>5 else ''}")
        print(f"  Restricting all aggregates to the {len(intersection)}-entry intersection.")
    return intersection


def aggregate(results_dict, entry_subset, common_set):
    """Aggregate metrics over the intersection of entry_subset and common_set."""
    if not results_dict or not entry_subset:
        return None
    keys_use = [e for e in entry_subset if e in common_set and e in results_dict]
    if not keys_use:
        return None
    agg = {}
    sample = results_dict[keys_use[0]]
    for k in sample.keys():
        vals = []
        for e in keys_use:
            v = results_dict[e].get(k)
            if v is None: continue
            try:
                vf = float(v)
                if not np.isnan(vf): vals.append(vf)
            except (TypeError, ValueError):
                continue
        if vals:
            agg[k] = {"mean": float(np.mean(vals)), "median": float(np.median(vals)), "n": len(vals)}
    return agg


def compute_paired_stats(baseline_res, augmented_res, subset, common_set, bootstrap_iters):
    """Per-entry paired Δ-metric stats, returned as a dict per metric.

    F9 FIX (2026-10-15). This used to convert its t-statistic through the NORMAL
    CDF:

        p = 2 * (1 - 0.5 * (1 + erf(abs(t) / sqrt(2))))

    At the n=5 most on-target subsets run, that understates p by 2.6x (t=2.0) to
    14.8x (t=3.0), so every "significant" subset p-value this script ever
    printed was overstated. It now delegates to analysis/common/stats.py, which
    uses the t-distribution with df = n-1.

    Two further fixes in the same pass:
      * The bootstrap RNG was re-seeded to 42 INSIDE the per-metric loop, so all
        six metrics resampled on identical indices and their CIs were perfectly
        correlated. Each metric now gets its own deterministic-but-distinct seed.
      * Entries are per-CHAIN, so N chains of one crystal are not N independent
        observations. n_structures is now reported alongside n_pairs; when they
        differ the paired test is pseudoreplicated.

    The output schema is a superset of the old one -- nothing downstream that
    read n_pairs/mean_delta/ci95_lo/ci95_hi/t_stat/p_value/deltas breaks.
    """
    keys_use = [e for e in subset if e in common_set
                and e in baseline_res and e in augmented_res]
    if len(keys_use) < 2:
        return None

    n_entries, n_struct = effective_n(keys_use)

    out = {}
    sample = baseline_res[keys_use[0]]
    for m_i, k in enumerate(sorted(sample.keys())):
        b_vals, a_vals, used = [], [], []
        for e in keys_use:
            try:
                bv = float(baseline_res[e].get(k, float("nan")))
                av = float(augmented_res[e].get(k, float("nan")))
                if not (np.isnan(bv) or np.isnan(av)):
                    b_vals.append(bv); a_vals.append(av); used.append(e)
            except (TypeError, ValueError):
                continue
        if len(b_vals) < 2:
            continue

        # Distinct per-metric seed so the six bootstrap CIs are independent
        # rather than sharing one resample draw.
        res = paired_delta(b_vals, a_vals, n_boot=bootstrap_iters,
                           seed=42 + 1000 * m_i)
        _, k_struct = effective_n(used)

        out[k] = {
            # original schema, preserved
            "n_pairs": res["n_pairs"],
            "mean_delta": res["mean_delta"],
            "ci95_lo": res["ci_lo"],
            "ci95_hi": res["ci_hi"],
            "t_stat": res["t_stat"],
            "p_value": res["p_value"],
            "deltas": (np.asarray(a_vals) - np.asarray(b_vals)).tolist(),
            # added
            "median_delta": res["median_delta"],
            "sd": res["sd"],
            "cohens_dz": res["cohens_dz"],
            "n_structures": k_struct,
            "pseudoreplicated": bool(k_struct < res["n_pairs"]),
            "p_method": "two-sided one-sample t on paired deltas, df=n-1",
        }

    # BH-FDR across the metrics tested for this pair. Reported separately from
    # p_value and deliberately labelled: these metrics are strongly correlated
    # (pearson/spearman/ic_weighted_pcc all measure the same agreement), so this
    # is a floor on the multiplicity problem, not a complete treatment. The
    # multiplicity that matters more -- across pilots -- has to be applied by
    # whoever aggregates these JSONs.
    if out:
        mkeys = sorted(out)
        qs = bh_fdr([out[k]["p_value"] for k in mkeys])
        for k, q in zip(mkeys, qs):
            out[k]["q_value_across_metrics"] = (None if np.isnan(q) else float(q))
        out["_meta"] = {
            "n_entries": n_entries,
            "n_structures": n_struct,
            "pseudoreplicated": bool(n_struct < n_entries),
            "fdr_scope": "BH across the metrics in this pair only; NOT across "
                         "pilots or seeds",
        }
    return out


def serializable(o):
    if isinstance(o, dict):
        return {k: serializable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [serializable(v) for v in o]
    if isinstance(o, (np.floating, np.integer)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return o


def main():
    args = parse_args()
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")

    conditions = load_conditions(args.condition)

    with open(args.id_file) as f:
        entries = [line.strip() for line in f if line.strip()]
    print(f"Loaded {len(entries)} entries from {args.id_file}")

    pat = re.compile(args.filter)
    subset_entries = sorted([e for e in entries if pat.search(e)])
    print(f"{args.filter_name} subset (filter '{args.filter}'): n={len(subset_entries)}")
    for e in subset_entries:
        print(f"  {e}")

    results = {}
    for cond_name, cond_dir in conditions.items():
        print(f"\n===== {cond_name} ({cond_dir}) =====")
        r = evaluate_condition(cond_name, cond_dir, entries, args.combined_dir, device)
        if r is None:
            print("  SKIPPED")
            continue
        results[cond_name] = r
        print(f"  evaluated {len(r)} entries")

    common_set = check_denominators(results, entries)

    KEYS = ['auroc', 'mae', 'pearsonr', 'ic_weighted_pcc', 'spearmanr', 'brier_multi']

    def print_block(title, entry_subset):
        print(f"\n--- {title} ---")
        for cond, r in results.items():
            agg = aggregate(r, entry_subset, common_set)
            if not agg:
                continue
            line = f"  {cond:30s}  "
            for k in KEYS:
                if k in agg:
                    line += f"{k}={agg[k]['mean']:.4f}  "
            ns = set(agg[k]['n'] for k in KEYS if k in agg)
            line += f"(n={ns.pop() if len(ns)==1 else sorted(ns)})"
            print(line)

    print("\n========== RESULTS ==========")
    print_block("All id.txt entries (intersection across conditions)", common_set)
    print_block(f"{args.filter_name} entries only", subset_entries)

    # Per-entry detail
    print(f"\n--- Per-{args.filter_name}-entry detail ---")
    if results and subset_entries:
        cond_names = list(results.keys())
        DETAIL_KEYS = ['auroc', 'mae', 'pearsonr', 'ic_weighted_pcc']
        hdr = f"  {'entry':45s} "
        for k in DETAIL_KEYS:
            for cn in cond_names:
                hdr += f"{cn[:12]}_{k[:6]:>6s} "
        print(hdr)
        for entry in subset_entries:
            if entry not in common_set:
                continue
            row = f"  {entry:45s} "
            for k in DETAIL_KEYS:
                for cn in cond_names:
                    v = results[cn].get(entry, {}).get(k)
                    try: vf = float(v)
                    except (TypeError, ValueError): vf = float("nan")
                    row += f"{vf:>13.4f}  "
            print(row)

    # =========================================================================
    # Paired statistics: for each baseline_X / augmented_X pair, paired Δ stats
    # =========================================================================
    print("\n--- Paired effect-size statistics ---")
    print(f"  (paired bootstrap 95% CI from {args.bootstrap_iters} samples; "
          f"p-value from paired one-sample t-test on Δ)")
    # Pair conditions by suffix; warn if any condition is unpaired.
    # We recognize three kinds: baseline, augmented, augmented_legacy.
    # augmented_legacy_<X> pairs against baseline_<X> (same baselines as the
    # production augmented_<X>). This supports A/B testing two augmenting
    # protocols against a single baseline.
    seed_pairs = {}
    unpaired = []
    for cond in results:
        # Try the longer prefix first so 'augmented_legacy_<X>' isn't
        # mis-matched as 'augmented' with suffix 'legacy_<X>'.
        m = re.match(r"^(augmented_legacy|baseline|augmented)_(.+)$", cond)
        if m:
            kind, suffix = m.groups()
            seed_pairs.setdefault(suffix, {})[kind] = cond
        else:
            unpaired.append(cond)

    if unpaired:
        print(f"  ⚠ {len(unpaired)} unpaired condition(s) excluded from paired stats: "
              f"{unpaired}")

    # Each suffix must have a baseline. An "augmented" or "augmented_legacy"
    # paired with that baseline is a pairwise comparison.
    incomplete_pairs = []
    for suffix, pair in seed_pairs.items():
        has_baseline = "baseline" in pair
        has_any_aug = "augmented" in pair or "augmented_legacy" in pair
        if not (has_baseline and has_any_aug):
            incomplete_pairs.append(suffix)
    if incomplete_pairs:
        print(f"  ⚠ {len(incomplete_pairs)} incomplete pair(s) excluded "
              f"(missing baseline or any augmented variant): {incomplete_pairs}")

    paired_stats_all = {}
    for suffix, pair in seed_pairs.items():
        if "baseline" not in pair:
            continue
        b_name = pair["baseline"]
        # Compute stats for each augmented variant present
        for aug_kind in ("augmented", "augmented_legacy"):
            if aug_kind not in pair:
                continue
            a_name = pair[aug_kind]
            stats = compute_paired_stats(
                results[b_name], results[a_name],
                subset_entries, common_set, args.bootstrap_iters
            )
            if stats is None:
                print(f"  {suffix}/{aug_kind}: too few paired entries; skipped")
                continue
            # Use a key that distinguishes which augmented variant this is
            pair_key = f"{aug_kind}:{suffix}"
            paired_stats_all[pair_key] = stats
            meta = stats.get("_meta", {})
            ps_note = ""
            if meta.get("pseudoreplicated"):
                ps_note = (f"  [PSEUDOREPLICATED: {meta['n_entries']} entries span "
                           f"only {meta['n_structures']} structure(s) — effective "
                           f"n is {meta['n_structures']}]")
            print(f"\n  Pair: {b_name} vs {a_name} ({args.filter_name} subset){ps_note}")
            # F8 FIX (2026-10-15): this print used to sit OUTSIDE the `for k in
            # KEYS` loop (dedented one level), so only the last key in KEYS was
            # ever printed — every stage-7 log in slurm_output/ shows a single
            # Δbrier_multi line per pair and nothing else. The JSON always held
            # all six metrics; only the console report was truncated.
            for k in KEYS:
                if k not in stats:
                    continue
                s = stats[k]
                sig = "*" if s["p_value"] < 0.05 else " "
                q = s.get("q_value_across_metrics")
                qtxt = f"  q={q:.3f}" if q is not None else ""
                print(f"    Δ{k:18s}  mean={s['mean_delta']:+.4f}  "
                      f"95% CI [{s['ci95_lo']:+.4f}, {s['ci95_hi']:+.4f}]  "
                      f"p={s['p_value']:.3f}{sig}{qtxt}  "
                      f"dz={s.get('cohens_dz', float('nan')):+.2f}  "
                      f"(n={s['n_pairs']})")
            print("    p: two-sided t on paired deltas (df=n-1). "
                  "q: BH across these metrics only, which are correlated — "
                  "NOT across pilots.")

    # Cross-seed aggregate, grouped by augmented variant
    # (keys in paired_stats_all are "augmented:<suffix>" or "augmented_legacy:<suffix>")
    if len(paired_stats_all) >= 2:
        from collections import defaultdict
        by_kind = defaultdict(dict)
        for key, stats in paired_stats_all.items():
            kind, suffix = key.split(":", 1)
            by_kind[kind][suffix] = stats

        for kind in sorted(by_kind):
            stats_dict = by_kind[kind]
            if len(stats_dict) < 2:
                continue
            print(f"\n--- Cross-seed aggregate ({kind} vs baseline; mean of seed Δ-means) ---")
            for k in KEYS:
                seed_deltas = [s[k]['mean_delta'] for s in stats_dict.values() if k in s]
                if len(seed_deltas) >= 2:
                    cross_mean = np.mean(seed_deltas)
                    cross_sd = np.std(seed_deltas, ddof=1)
                    print(f"  Δ{k:18s}  mean across seeds={cross_mean:+.4f}  "
                          f"(SD={cross_sd:.4f}, n_seeds={len(seed_deltas)})")

        # If both augmented and augmented_legacy are present, also report
        # the head-to-head Δ (per-seed paired) — does legacy give larger Δ
        # than production at each seed?
        if "augmented" in by_kind and "augmented_legacy" in by_kind:
            shared_suffixes = sorted(set(by_kind["augmented"]) & set(by_kind["augmented_legacy"]))
            if shared_suffixes:
                print(f"\n--- Head-to-head: augmented vs augmented_legacy (per-seed paired Δ-of-Δs) ---")
                print(f"  For each seed, computes (Δ_augmented_legacy - Δ_augmented).")
                print(f"  Positive => legacy augmenting gave larger improvement over baseline.")
                for k in KEYS:
                    diffs = []
                    for suffix in shared_suffixes:
                        if k in by_kind["augmented"][suffix] and k in by_kind["augmented_legacy"][suffix]:
                            d_prod = by_kind["augmented"][suffix][k]['mean_delta']
                            d_leg  = by_kind["augmented_legacy"][suffix][k]['mean_delta']
                            diffs.append(d_leg - d_prod)
                    if len(diffs) >= 2:
                        cross_mean = np.mean(diffs)
                        cross_sd = np.std(diffs, ddof=1)
                        print(f"  Δ{k:18s}  legacy − production mean={cross_mean:+.4f}  "
                              f"(SD={cross_sd:.4f}, n_seeds={len(diffs)})")

    os.makedirs(os.path.dirname(args.output_json), exist_ok=True)
    output_payload = {
        "results": results,
        "paired_stats": paired_stats_all,
        "common_entries": sorted(common_set),
        "subset_entries": subset_entries,
    }
    with open(args.output_json, "w") as f:
        json.dump(serializable(output_payload), f, indent=2)
    print(f"\nFull results saved: {args.output_json}")


if __name__ == "__main__":
    main()