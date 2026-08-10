#!/usr/bin/env python3
"""
make_clash_trajectory.py — Ch. 9.1 figure: the vdW-ramp clash trajectory.

Parses a Stage-3 minimization log (the engine prints, per representative state,
heavy-atom protein-DNA clash count and min contact distance at: Initial, each
vdW-ramp stage sigma=0.1..1.0, and Final) and plots both series across
Phase 0 -> ramp -> Final. Shows how the sigma-ramp resolves steric clashes
gently instead of in one hard minimization.

Engine lines parsed (grep-able):
  [tag] Initial heavy-atom clashes: N, min_dist: D A
  [tag]   Stage i/n: sigma=s, PE=..., clashes=N, min_d=D A
  [tag]   Final: PE=..., clashes=N, min_d=D A
Phase 0 (H-only min) freezes heavy atoms, so the heavy-atom clash count is
unchanged from Initial and is carried forward as the 'H-min' point.

Usage:  python make_clash_trajectory.py --log <stage3_log.out> [--out F.png]
Env: any with matplotlib (deeppbs). Stdlib parse + matplotlib.
"""
import argparse
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TEAL = "#2E7FA8"
ALARM = "#D1495B"
RAMP_BAND = "#F2F4F6"


def parse_log(path):
    init_c = init_d = None
    stages = []
    final = None
    for ln in open(path):
        m = re.search(r"Initial heavy-atom clashes:\s*(\d+),\s*min_dist:\s*([\d.]+)", ln)
        if m:
            init_c, init_d = int(m.group(1)), float(m.group(2)); continue
        m = re.search(r"Stage\s+\d+/\d+:\s*σ=([\d.]+),.*clashes=(\d+),\s*min_d=([\d.]+)", ln)
        if not m:
            m = re.search(r"Stage\s+\d+/\d+:\s*sigma=([\d.]+),.*clashes=(\d+),\s*min_d=([\d.]+)", ln)
        if m:
            stages.append(("σ=%s" % m.group(1), int(m.group(2)), float(m.group(3)))); continue
        m = re.search(r"Final:.*clashes=(\d+),\s*min_d=([\d.]+)", ln)
        if m:
            final = (int(m.group(1)), float(m.group(2)))
    if init_c is None or not stages or final is None:
        raise SystemExit("could not parse clash trajectory from %s" % path)
    labels = ["init", "H-min"] + [s[0] for s in stages] + ["final"]
    clashes = [init_c, init_c] + [s[1] for s in stages] + [final[0]]
    mind = [init_d, init_d] + [s[2] for s in stages] + [final[1]]
    return labels, clashes, mind


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--log", required=True)
    ap.add_argument("--out", default="F_clash_trajectory.png")
    ap.add_argument("--state-label", default="FOXA 1vtn, state 001")
    a = ap.parse_args()
    labels, clashes, mind = parse_log(a.log)
    xi = np.arange(len(labels))

    fig, ax1 = plt.subplots(figsize=(6.8, 3.4))
    ax1.axvspan(1.5, len(labels) - 1.5, color=RAMP_BAND, zorder=0)
    # band label centered over the early-ramp plateau, above the clash line
    ax1.text(3.0, max(clashes) * 1.11,
             "vdW σ-ramp (Phase 1)", ha="center", va="top", fontsize=6.5, color="#888")
    ax1.plot(xi, clashes, "-o", color=ALARM, lw=1.8, ms=5, zorder=3)
    ax1.set_ylabel("protein–DNA heavy-atom clashes", color=ALARM)
    ax1.tick_params(axis="y", colors=ALARM)
    ax1.set_ylim(-max(clashes) * 0.06, max(clashes) * 1.18)
    ax2 = ax1.twinx()
    ax2.plot(xi, mind, "-s", color=TEAL, lw=1.8, ms=5, zorder=3)
    ax2.set_ylabel("min. protein–DNA distance (Å)", color=TEAL)
    ax2.tick_params(axis="y", colors=TEAL)
    ax2.set_ylim(0, max(mind) * 1.15)
    ax2.spines["top"].set_visible(False)
    ax1.set_xticks(xi)
    ax1.set_xticklabels(labels, fontsize=6.5)
    ax1.set_xlabel("minimization phase")
    ax1.set_title("The vdW ramp resolves clashes gently: %d → %d as atoms grow to full radius"
                  % (clashes[0], clashes[-1]), fontsize=8, loc="left")
    fig.text(0.01, -0.04,
             "One representative state (%s). Phase 0 (H-only min) freezes heavy atoms, so the "
             "heavy-atom clash count is unchanged from init; the σ-ramp (Phase 1) grows vdW radii "
             "0.1→1.0 over five stages, resolving clashes and opening the min contact distance; "
             "Phase 2 is the final full-radius minimization." % a.state_label,
             fontsize=6.0, ha="left", va="top", wrap=True)
    fig.savefig(a.out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote %s (clashes %d→%d, min_d %.2f→%.2f Å)"
          % (a.out, clashes[0], clashes[-1], mind[0], mind[-1]))


if __name__ == "__main__":
    main()
