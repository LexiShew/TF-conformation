#!/usr/bin/env python3
"""Regenerate predicted-PWM .npz files for the AB-test scored entries.

Background
----------
The `id_benchmark_<tf>_legacy_ab.json` files were produced by
`evaluate_id_benchmark.py`, which re-runs each checkpoint in-memory over the
`folds/id.txt` benchmark set and keeps only scalar metrics. It never saved
the predicted PWM arrays. The `predictions/` folder inside each run dir is a
*training* side-effect (driver.py dumps predictions over that run's own
validation fold), so it covers a different, smaller set than the AB benchmark
— the exact AB entries' predicted PWMs are nowhere on disk.

This script closes that gap. For every run dir (same scan as
eval_legacy_ab.sh) it loads `Model.best.tar` + `scaler.pkl` + `config.json`,
runs the entries in `--id-file` through the evaluator, and writes
`<entry>_predict.npz` in the exact same key layout driver.py uses
(Y, Y_mask, P, P_mask, L, P_rc, P_rc_mask). Those files are then directly
renderable by scripts/generate_pwm_logos.py.

Must run where the deeppbs package, checkpoints and combined_assembly data
live (the cluster), inside the `deeppbs` conda env.

Example
-------
    python predict_ab_pwms.py \
        --tf tbp \
        --id-file /path/to/id_tbp_ab.txt \
        --outputs-dir "$OUTPUTS_DIR" \
        --combined-dir "$DATA_DIR/combined_assembly_tbp"
"""

import argparse
import glob
import json
import os
import pickle
import sys

import numpy as np
import torch
from torch_geometric.data import DataLoader

# Self-locate the repo: deeppbs from lib/, models.model_v2 from stage6_train/.
_TFCONF = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_TFCONF, "stage6_train"))  # models.model_v2
sys.path.insert(0, os.path.join(_TFCONF, "lib"))           # deeppbs package

from deeppbs.nn.utils import loadDataset          # noqa: E402
from deeppbs.nn import Evaluator                   # noqa: E402
from models.model_v2 import Model                  # noqa: E402

NC = 4
LABELS_KEY = "Y_pwm"
NF_PROT, NF_DNA = 13, 14


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tf", required=True,
                   help="Pilot name (tbp / dux4 / egr1); used for defaults.")
    p.add_argument("--id-file", required=True,
                   help="Pinned AB entry list (one entry per line, .npz kept).")
    p.add_argument("--combined-dir",
                   help="Feature dir holding the entries' npz "
                        "(default: $DATA_DIR/combined_assembly_<tf>).")
    p.add_argument("--outputs-dir", default=os.environ.get("OUTPUTS_DIR"),
                   help="Dir holding <kind>_<tf>_fold<F>_s* run dirs "
                        "(default: $OUTPUTS_DIR).")
    p.add_argument("--fold", default="0")
    p.add_argument("--kinds", nargs="+",
                   default=["baseline", "augmented", "augmented_legacy"])
    p.add_argument("--pred-subdir", default="predictions_ab",
                   help="Subdir name written inside each run dir. Use "
                        "'predictions' to feed generate_pwm_logos.py as-is "
                        "(WARNING: clobbers the training-validation preds).")
    p.add_argument("--threshold", type=float, default=0.5)
    p.add_argument("--overwrite", action="store_true",
                   help="Re-write entries even if the npz already exists.")
    return p.parse_args()


def find_run_dirs(outputs_dir, kinds, tf, fold):
    """Mirror eval_legacy_ab.sh: outer <kind>_<tf>_fold<F>_s*/ then the inner
    subdir that actually contains Model.best.tar."""
    runs = []
    for kind in kinds:
        pattern = os.path.join(outputs_dir, f"{kind}_{tf}_fold{fold}_s*")
        for outer in sorted(glob.glob(pattern)):
            if not os.path.isdir(outer):
                continue
            for sub in sorted(glob.glob(os.path.join(outer, "*", ""))):
                if os.path.isfile(os.path.join(sub, "Model.best.tar")):
                    runs.append(sub.rstrip("/"))
    return runs


def predict_run(run_dir, entries, combined_dir, pred_subdir, threshold,
                overwrite, device):
    with open(os.path.join(run_dir, "scaler.pkl"), "rb") as f:
        scaler = pickle.load(f)
    with open(os.path.join(run_dir, "config.json")) as f:
        train_C = json.load(f)
    condition = train_C.get("condition", "prot_shape")
    readout = train_C.get("readout", "all")

    dataset, _, _, found_files = loadDataset(
        list(entries), NC, LABELS_KEY, combined_dir,
        cache_dataset=False, balance="unmasked", remove_mask=False,
        scale=True, scaler=scaler, pre_transform=None, feature_mask=None,
    )
    DL = DataLoader(dataset, batch_size=1, shuffle=False, pin_memory=True)

    model = Model(NF_PROT, NF_DNA, condition=condition, readout=readout)
    state = torch.load(os.path.join(run_dir, "Model.best.tar"),
                       map_location=device)
    model.load_state_dict(state["model_state_dict"])
    model.to(device)

    evaluator = Evaluator(model, NC, device=device,
                          post_process=torch.nn.Softmax(dim=-1))
    val_out = evaluator.eval(
        DL, use_mask=False, batchwise=True, return_masks=True,
        return_predicted=False, return_batches=True, xtras=None,
        threshold=threshold, eval_mode=True,
    )

    out_dir = os.path.join(run_dir, pred_subdir)
    os.makedirs(out_dir, exist_ok=True)

    n_written = 0
    n_skipped = 0
    for i in range(val_out["num_batches"]):
        name = found_files[i]                       # e.g. 2evi_A_MA0343.1.jaspar.npz
        out_path = os.path.join(out_dir, f"{name}_predict.npz")
        if os.path.exists(out_path) and not overwrite:
            n_skipped += 1
            continue

        # Reshape exactly as stage6_train/driver.py does before saving.
        y = val_out["y"][i]
        prob = val_out["output"][i]
        mask = val_out["masks"][i]
        idx = val_out["indexes"][i]
        out_mask = val_out["out_masks"][i]
        out_idx = val_out["out_idx"][i]
        logits = val_out["logits"][i]

        y = y[idx]
        mask = mask[idx]
        prob_rc = prob[(1 - out_idx).astype(bool)]
        prob = prob[out_idx]
        out_mask_rc = out_mask[(1 - out_idx).astype(bool)]
        out_mask = out_mask[out_idx]
        logits = logits[out_idx]

        np.savez_compressed(
            out_path,
            Y=y, Y_mask=mask, P=prob, P_mask=out_mask,
            L=logits, P_rc=prob_rc, P_rc_mask=out_mask_rc,
        )
        n_written += 1

    return n_written, n_skipped, len(found_files), out_dir


def main():
    args = parse_args()

    if not args.outputs_dir:
        sys.exit("ERROR: --outputs-dir not given and $OUTPUTS_DIR unset")
    # Prefer the per-pilot $COMBINED_ASSEMBLY_DIR set by lib/common.sh (now under
    # the repo's output/stage5_aug/), falling back to the combined dir for --tf.
    combined_dir = args.combined_dir or os.environ.get("COMBINED_ASSEMBLY_DIR")
    if not combined_dir:
        sys.exit("ERROR: --combined-dir not given and $COMBINED_ASSEMBLY_DIR unset")
    if not os.path.isdir(combined_dir):
        sys.exit(f"ERROR: combined-dir does not exist: {combined_dir}")

    with open(args.id_file) as f:
        entries = [ln.strip() for ln in f if ln.strip()]
    print(f"Loaded {len(entries)} AB entries from {args.id_file}")

    runs = find_run_dirs(args.outputs_dir, args.kinds, args.tf, args.fold)
    if not runs:
        sys.exit(f"ERROR: no run dirs matched under {args.outputs_dir} "
                 f"for kinds={args.kinds} tf={args.tf} fold={args.fold}")
    print(f"Found {len(runs)} run dir(s):")
    for r in runs:
        print(f"  {r}")

    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}  |  combined-dir: {combined_dir}\n")

    grand_written = 0
    for run_dir in runs:
        name = os.path.basename(run_dir)
        try:
            w, s, n, out_dir = predict_run(
                run_dir, entries, combined_dir, args.pred_subdir,
                args.threshold, args.overwrite, device)
            grand_written += w
            print(f"[{name}] wrote {w}, skipped {s} (existing), "
                  f"loaded {n} entries -> {out_dir}")
        except Exception as e:           # keep going across runs
            print(f"[{name}] FAILED: {type(e).__name__}: {e}", file=sys.stderr)

    print(f"\nDone. Wrote {grand_written} prediction npz file(s) "
          f"across {len(runs)} run(s).")


if __name__ == "__main__":
    main()
