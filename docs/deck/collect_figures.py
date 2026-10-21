#!/usr/bin/env python
"""Collect the figures the v5 deck needs into one directory under short names.

deck_spec_v5.json refers to images as N1.jpg, pca_rok.jpg and so on. This maps
those names onto their real locations in the repo and writes downscaled JPEGs
(1600 px wide, quality 88) -- about 4 MB total, which keeps the built deck near
10 MB instead of carrying full 300-dpi PNGs.

Run from the repo root:

    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
    /project2/rohs_102/shewchuk/conda/envs/deeppbs/bin/python \
        docs/deck/collect_figures.py --out docs/deck/deckimg

(the thread pin is needed on an endeavour login node, same as make_talk_figs.py)

Regenerate N1-N10 first if the tables changed:
    cd analysis/analyses/nmr_rmsd && python make_talk_figs.py
"""
import argparse
import os

from PIL import Image

NMR = "analysis/analyses/nmr_rmsd"
FIG = "analysis/figures"
MECH = "analysis/analyses/mechanism"

# short name used in deck_spec_v5.json  ->  path relative to the repo root
SOURCES = {
    "N1.jpg": f"{NMR}/plots/talk/N1_three_way_fnat.png",
    "N2.jpg": f"{NMR}/plots/talk/N2_pair_roster.png",
    "N3.jpg": f"{NMR}/plots/talk/N3_shared_core.png",
    "N4.jpg": f"{NMR}/plots/talk/N4_pct_across.png",
    "N5.jpg": f"{NMR}/plots/talk/N5_gate.png",
    "N6.jpg": f"{NMR}/plots/talk/N6_residual_ratio.png",
    "N7.jpg": f"{NMR}/plots/talk/N7_pilot_fnat.png",
    "N8.jpg": f"{NMR}/plots/talk/N8_augmentation.png",
    "N9.jpg": f"{NMR}/plots/talk/N9_three_way_design.png",
    "N10.jpg": f"{NMR}/plots/talk/N10_construct_control.png",
    "B1.jpg": f"{MECH}/figures/B1_bend_decomposition.png",
    "B2.jpg": f"{MECH}/figures/B2_bend_localization.png",
    "B3.jpg": f"{MECH}/figures/B3_bend_reconstruction.png",
    "N11.jpg": f"{NMR}/plots/talk/N11_fnat_vs_separation.png",
    "N12.jpg": f"{NMR}/plots/talk/N12_protein_x_dna.png",
    "pca_nhp6a.jpg": f"{NMR}/plots/pair_pca_nhp6a_ca.png",
    "pca_c20.jpg": f"{NMR}/plots/pair_pca_c20_ca.png",
    "pca_c10.jpg": f"{NMR}/plots/pair_pca_c10_ca.png",
    "pca_rok.jpg": f"{NMR}/plots/pair_pca_rok_ca.png",
    "pca_c05.jpg": f"{NMR}/plots/pair_pca_c05_ca.png",
    "trim_sweep.jpg": f"{NMR}/plots/trim_sweep_ca.png",
    "rmsd_ragged.jpg": f"{NMR}/plots/rmsd_apo_holo_ca.png",
    "rmsd_core.jpg": f"{NMR}/plots/rmsd_apo_holo_ca_core.png",
    "side_by_side.jpg": f"{NMR}/plots/pair_side_by_side_all_ca_core.png",
    "F1_fnat.jpg": f"{FIG}/F1_fnat_distributions.png",
    "D1_diversity.jpg": f"{FIG}/D1_diversity.png",
    "M1_mechanism.jpg": f"{FIG}/M1_apo_holo_mechanism.png",
    "af3_mgwfl.jpg": "analysis/analyses/dna_relax/figures/af3_vs_ensemble_mgwfl.png",
    "struct_pwm.jpg": "analysis/analyses/struct_pwm/figures/struct_pwm_context.png",
}

# Deliberately NOT in the deck, recorded so the next person does not re-add them:
#   pair_pca_all_ca.png -- all 16 pairs stacked vertically; at 16:9 it fits to
#       0.70 in wide and is unreadable. Slide 15 states the variance shares as
#       text instead.
#   M1_apo_holo_mechanism.png is backup only: it is reachability vs augmentation
#       outcome, not the apo-fnat mechanism read-out that slide 22 makes.
#   rmsd_apo_holo_ca{,_core}.png are backup only: they measure ensemble-to-own-
#       bundle RMSD, not apo-holo separation.


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", default=".", help="repo root (default: cwd)")
    ap.add_argument("--out", default="docs/deck/deckimg")
    ap.add_argument("--width", type=int, default=1600)
    ap.add_argument("--quality", type=int, default=88)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    missing, total = [], 0
    for name, rel in sorted(SOURCES.items()):
        src = os.path.join(a.repo, rel)
        if not os.path.exists(src):
            missing.append(rel)
            continue
        im = Image.open(src)
        w, h = im.size
        scale = min(1.0, a.width / w)
        dst = os.path.join(a.out, name)
        im.convert("RGB").resize((int(w * scale), int(h * scale)),
                                 Image.LANCZOS).save(dst, quality=a.quality)
        total += os.path.getsize(dst)
    print("wrote %d/%d images to %s (%.2f MB)"
          % (len(SOURCES) - len(missing), len(SOURCES), a.out, total / 1e6))
    if missing:
        raise SystemExit("MISSING:\n  " + "\n  ".join(missing))


if __name__ == "__main__":
    main()
