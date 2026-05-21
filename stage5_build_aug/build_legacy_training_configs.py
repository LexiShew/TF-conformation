"""
Build training configs for the legacy A/B test.

Only writes augmented_legacy configs (baseline configs already exist from the
prior multi-seed run and don't depend on the Stage 3 protocol).

Usage:
  python build_legacy_training_configs.py \\
      --tf-name egr1 \\
      --combined-dir /project2/.../combined_assembly_legacy_egr1 \\
      --seeds 1 2 3 4 5 \\
      --output-dir /project2/.../DeepPBS/run
"""
import argparse
import json
import os

OUTPUTS_ROOT = "/project2/rohs_102/shewchuk/DeepPBS_outputs"


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
        # Same seeds as the existing multi-seed run so we can pair against
        # the existing baseline checkpoints exactly.
        "random_seed": random_seed,
        "no_random": False,
        "model": {"transform_args": []},
        "optimizer": {"name": "adam", "kwargs": {"lr": 0.001, "weight_decay": 0.0001}},
        "scheduler": {"name": "", "kwargs": {}},
    }


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tf-name", required=True)
    p.add_argument("--combined-dir", required=True,
                   help="combined_assembly_legacy_<tf>/ directory")
    p.add_argument("--fold", type=int, default=0)
    p.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3, 4, 5],
                   help="Seeds to pair against existing baseline checkpoints")
    p.add_argument("--output-dir", required=True)
    args = p.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    for seed in args.seeds:
        suffix = f"_s{seed}"
        cfg = base_config(
            data_dir=args.combined_dir,
            output_path=f"{OUTPUTS_ROOT}/augmented_legacy_{args.tf_name}_fold{args.fold}{suffix}",
            random_seed=seed,
        )
        path = os.path.join(args.output_dir,
                            f"config_augmented_legacy_{args.tf_name}_fold{args.fold}{suffix}.json")
        with open(path, "w") as f:
            json.dump(cfg, f, indent=2)
        print(f"Wrote {path}  (seed={seed})")


if __name__ == "__main__":
    main()