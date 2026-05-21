"""
Build an augmented training-fold file by appending newly-processed npz entries
to an original train fold. Enforces a trailing newline on the original file
to avoid the line-concatenation bug we hit during the EGR1 pilot.

Filters new entries by a substring (typically the PWM label) to ensure the
augmentation matches the test set's PWM convention.

Example:
  python build_augmented_fold.py \\
      --orig-train /project2/.../folds/train0.txt \\
      --stage4-dir /project2/.../conformations/dux4/stage4_npz/output \\
      --pwm-filter DUX4_HUMAN.H11MO.0.A \\
      --out-train  /project2/.../folds_aug/train0_aug_dux4.txt
"""
import argparse
import os
import sys


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--orig-train", required=True, help="Original train fold file to copy")
    p.add_argument("--stage4-dir", required=True, help="Directory with new npz files")
    p.add_argument("--pwm-filter", required=True,
                   help="Substring filter for selecting which new entries to include "
                        "(typically the PWM label; e.g. 'DUX4_HUMAN.H11MO.0.A')")
    p.add_argument("--out-train", required=True, help="Output augmented train fold")
    p.add_argument("--no-overwrite", action="store_true",
                   help="Refuse to overwrite an existing --out-train")
    args = p.parse_args()

    if args.no_overwrite and os.path.exists(args.out_train):
        print(f"ERROR: {args.out_train} already exists and --no-overwrite was set", file=sys.stderr)
        sys.exit(1)

    # Read original
    with open(args.orig_train) as f:
        orig_content = f.read()

    if not orig_content.endswith("\n"):
        orig_content += "\n"
    orig_lines = [l for l in orig_content.splitlines() if l.strip()]

    # Find new entries, filtered and sorted
    new_entries = sorted(
        f for f in os.listdir(args.stage4_dir)
        if f.endswith(".npz") and args.pwm_filter in f
    )

    if not new_entries:
        print(f"ERROR: no .npz files matched filter '{args.pwm_filter}' in {args.stage4_dir}", file=sys.stderr)
        sys.exit(1)

    # Sanity check: warn if any of the new entries appear in the original
    # (would indicate accidental leakage / double-counting)
    orig_set = set(orig_lines)
    overlap = [e for e in new_entries if e in orig_set]
    if overlap:
        print(f"WARNING: {len(overlap)} new entries already appear in {args.orig_train}", file=sys.stderr)
        print(f"  e.g. {overlap[:3]}", file=sys.stderr)
        print(f"  These will produce duplicates in the output. Consider deduplicating.", file=sys.stderr)

    # Write
    os.makedirs(os.path.dirname(args.out_train), exist_ok=True)
    with open(args.out_train, "w") as f:
        f.write(orig_content)
        for e in new_entries:
            f.write(e + "\n")

    n_orig = len(orig_lines)
    n_new = len(new_entries)
    n_total = n_orig + n_new
    print(f"Wrote {args.out_train}")
    print(f"  Original entries: {n_orig}")
    print(f"  New entries:      {n_new} (matching '{args.pwm_filter}')")
    print(f"  Total:            {n_total}")


if __name__ == "__main__":
    main()
