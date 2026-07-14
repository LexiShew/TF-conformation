#!/usr/bin/env python3
"""Assemble the six per-pilot two-panel PyMOL renders into a montage.

Layout: 6 rows (one per pilot, in fnat pass-rate order) x 2 columns
(AF3 | BioEmu). Reads analysis/figures/pymol/<pilot>_{af3,bioemu}.png,
writes analysis/figures/pymol/ensembles_montage.png. Panels are downscaled
so the whole montage is a manageable size.

Run in any env with Pillow (e.g. deeppbs). Cap BLAS threads on the login node.
"""
import os
from PIL import Image, ImageDraw, ImageFont

BASE = "/project2/rohs_102/shewchuk/TF-conformation"
PDIR = os.path.join(BASE, "analysis", "figures", "pymol")

# fnat pass-rate order (README)
ORDER = ["ets1", "tbp", "egr1", "engrailed", "foxa", "lef1"]
LABELS = {"ets1": "ETS1", "tbp": "TBP", "egr1": "EGR1",
          "engrailed": "engrailed", "foxa": "FOXA", "lef1": "LEF1"}

PANEL = 460       # scaled panel side (px)
PAD = 16          # gap between panels
TOP = 132         # header band (title + column headers + legend)
ROWLAB = 128      # left strip for pilot names

C_CRYST = (64, 64, 64)
C_AF3 = (185, 167, 214)
C_BIO = (76, 154, 168)
C_DNA = (203, 185, 148)


def font(sz, bold=False):
    cands = (["/usr/share/fonts/dejavu/DejaVuSans-Bold.ttf",
              "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"] if bold
             else ["/usr/share/fonts/dejavu/DejaVuSans.ttf",
                   "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"])
    for p in cands:
        if os.path.isfile(p):
            return ImageFont.truetype(p, sz)
    return ImageFont.load_default()


def load(pilot, tag):
    fp = os.path.join(PDIR, f"{pilot}_{tag}.png")
    if not os.path.isfile(fp):
        return None
    im = Image.open(fp).convert("RGB")
    return im.resize((PANEL, PANEL), Image.LANCZOS)


def main():
    nrow = len(ORDER)
    W = ROWLAB + 2 * PANEL + 3 * PAD
    H = TOP + nrow * PANEL + (nrow + 1) * PAD
    canvas = Image.new("RGB", (W, H), "white")
    draw = ImageDraw.Draw(canvas)

    f_title = font(26, bold=True)
    f_col = font(24, bold=True)
    f_row = font(24, bold=True)
    f_leg = font(18)

    # title centered on the full width, sized to fit
    title = "AlphaFold3 vs BioEmu conformational ensembles (superposed on crystal)"
    tw = draw.textlength(title, font=f_title)
    draw.text(((W - tw) / 2, 14), title, fill="black", font=f_title)

    # legend row, centered below the title
    leg = [("crystal", C_CRYST), ("crystal DNA", C_DNA),
           ("AF3 (10)", C_AF3), ("BioEmu (10)", C_BIO)]
    entries = [(n, col, draw.textlength(n, font=f_leg)) for n, col in leg]
    total = sum(26 + tw2 + 26 for _, _, tw2 in entries)
    x = (W - total) / 2
    for name, col, tw2 in entries:
        draw.rectangle([x, 54, x + 18, 72], fill=col, outline="black")
        draw.text((x + 24, 52), name, fill="black", font=f_leg)
        x += 26 + tw2 + 26

    # column headers
    col_x = [ROWLAB + PAD, ROWLAB + 2 * PAD + PANEL]
    for cx, ctitle in zip(col_x, ["AlphaFold3", "BioEmu"]):
        ctw = draw.textlength(ctitle, font=f_col)
        draw.text((cx + (PANEL - ctw) / 2, TOP - 34), ctitle, fill="black", font=f_col)

    # rows
    for r, p in enumerate(ORDER):
        py = TOP + PAD + r * (PANEL + PAD)
        # row label (vertical center)
        lab = LABELS[p]
        bbox = draw.textbbox((0, 0), lab, font=f_row)
        th = bbox[3] - bbox[1]
        draw.text((PAD, py + (PANEL - th) / 2 - bbox[1]), lab, fill="black", font=f_row)
        for cx, tag in zip(col_x, ["af3", "bioemu"]):
            im = load(p, tag)
            if im is not None:
                canvas.paste(im, (cx, py))
            else:
                draw.rectangle([cx, py, cx + PANEL, py + PANEL], outline="grey")
                draw.text((cx + 20, py + 20), f"missing {p}_{tag}", fill="grey", font=f_leg)

    out = os.path.join(PDIR, "ensembles_montage.png")
    canvas.save(out, dpi=(200, 200))
    print("wrote", out, canvas.size)


if __name__ == "__main__":
    main()
