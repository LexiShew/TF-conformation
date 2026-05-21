"""
Build a per-TF "combined assembly" directory containing symlinks to:
  - All original training npzs (assembly2024/)
  - New augmentation npzs filtered by PWM label

Each pilot gets its own combined dir so they don't interfere with each other.

Example:
  python build_combined_assembly.py \\
      --orig-dir /project2/.../assembly2024 \\
      --stage4-dir /project2/.../conformations/dux4/stage4_npz/output \\
      --pwm-filter DUX4_HUMAN.H11MO.0.A \\
      --out-dir /project2/.../combined_assembly_dux4
"""
import argparse
import os
import sys


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--orig-dir", required=True, help="Original assembly2024 directory")
    p.add_argument("--stage4-dir", required=True, help="Directory with new npz files")
    p.add_argument("--pwm-filter", required=True, help="Substring filter for new entries")
    p.add_argument("--out-dir", required=True, help="Output combined assembly directory")
    args = p.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # Symlink originals
    n_orig = 0
    for f in os.listdir(args.orig_dir):
        if f.endswith(".npz"):
            src = os.path.join(args.orig_dir, f)
            dst = os.path.join(args.out_dir, f)
            if not os.path.exists(dst):
                os.symlink(src, dst)
            n_orig += 1

    # Symlink filtered new entries
    n_new = 0
    n_skipped_existing = 0
    for f in os.listdir(args.stage4_dir):
        if not f.endswith(".npz"):
            continue
        if args.pwm_filter not in f:
            continue
        src = os.path.join(args.stage4_dir, f)
        dst = os.path.join(args.out_dir, f)
        if os.path.exists(dst):
            n_skipped_existing += 1
            continue
        os.symlink(src, dst)
        n_new += 1

    n_total = len(os.listdir(args.out_dir))
    print(f"Built {args.out_dir}")
    print(f"  Original entries (linked):      {n_orig}")
    print(f"  New entries (linked, filtered): {n_new}")
    if n_skipped_existing:
        print(f"  Skipped (already linked):       {n_skipped_existing}")
    print(f"  Total entries in dir:           {n_total}")

    if n_new == 0:
        print(f"WARNING: no new entries matched filter '{args.pwm_filter}'", file=sys.stderr)


if __name__ == "__main__":
    main()
