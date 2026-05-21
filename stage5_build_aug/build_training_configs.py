"""
Build per-TF training config files for baseline and augmented runs.

Both configs are identical except for `data_dir` and `output_path`, ensuring
the comparison is apples-to-apples. Crucially: a paired random seed is
written into both configs so baseline_seedN and augmented_seedN share
parameter init and shuffle order, isolating the augmentation effect.

Example:
  python build_training_configs.py \\
      --tf-name dux4 \\
      --combined-dir /project2/.../combined_assembly_dux4 \\
      --seed 1 \\
      --output-dir   /project2/.../DeepPBS/run

To run paired multi-seed comparisons:
  for s in 1 2 3 4 5; do
      python build_training_configs.py --tf-name X --combined-dir ... --seed $s \\
          --output-dir /tmp/configs_seed$s
  done
"""
import argparse
import json
import os


BASELINE_DATA_DIR = "/project2/rohs_102/shewchuk/DeepPBS_data/deeppbsmar24/data/assembly2024"
OUTPUTS_ROOT      = "/project2/rohs_102/shewchuk/DeepPBS_outputs"


def base_config(data_dir, output_path, random_seed):
    return {
        "data_dir": data_dir,
        "output_path": output_path,
        "nc": 4,
        "labels_key": "Y_pwm",
        "cache_dataset": False,
        "epochs": 50,
        "batch_size": 1,
        "loss": "soft_ce",
        "condition": "prot_shape",
        "ic_loss_weight": 0,
        "mse_loss_weight": 1,
        "remove_zero_class": False,
        "best_state_metric": "mae",
        "best_state_metric_goal": "min",
        "best_state_metric_threshold": 1.0,
        "best_state_metric_dataset": "validation",
        # CRITICAL: paired seed so baseline and augmented share parameter init
        # and shuffle order. The augmentation effect is then isolated from
        # init/shuffle noise.
        "random_seed": random_seed,
        "no_random": False,
        "model": {"transform_args": []},
        "optimizer": {"name": "adam", "kwargs": {"lr": 0.001, "weight_decay": 0.0001}},
        "scheduler": {"name": "", "kwargs": {}},
    }


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tf-name", required=True)
    p.add_argument("--combined-dir", required=True)
    p.add_argument("--fold", type=int, default=0)
    p.add_argument("--seed", type=int, default=42,
                   help="Random seed paired across baseline/augmented (default: 42)")
    p.add_argument("--output-dir", required=True,
                   help="Where to write config JSON files (typically run/)")
    p.add_argument("--seed-suffix", default="",
                   help="Optional suffix to append to output paths (e.g. '_seed3'). "
                        "If empty, no suffix is added.")
    args = p.parse_args()

    suffix = args.seed_suffix
    baseline_cfg = base_config(
        data_dir=BASELINE_DATA_DIR,
        output_path=f"{OUTPUTS_ROOT}/baseline_{args.tf_name}_fold{args.fold}{suffix}",
        random_seed=args.seed,
    )
    augmented_cfg = base_config(
        data_dir=args.combined_dir,
        output_path=f"{OUTPUTS_ROOT}/augmented_{args.tf_name}_fold{args.fold}{suffix}",
        random_seed=args.seed,
    )

    os.makedirs(args.output_dir, exist_ok=True)
    for name, cfg in (("baseline", baseline_cfg), ("augmented", augmented_cfg)):
        path = os.path.join(args.output_dir, f"config_{name}_{args.tf_name}_fold{args.fold}{suffix}.json")
        with open(path, "w") as f:
            json.dump(cfg, f, indent=2)
        print(f"Wrote {path}  (seed={args.seed})")


if __name__ == "__main__":
    main()
