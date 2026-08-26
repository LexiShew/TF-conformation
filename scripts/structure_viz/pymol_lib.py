"""pymol_lib.py — shared PyMOL rendering primitives for structure_viz.

Consolidates the former scripts/pymol/{color_protein_states, gradient_protein_states,
gradient_protein_split, spectrum_states}.py — four near-duplicate copies of the same
gradient-colouring routine (two were byte-identical; the others differed only in palette
and object-naming) — into one palette set + one function.

Load inside PyMOL, then use the extended commands (back-compat names preserved):
    run scripts/structure_viz/pymol_lib.py
    gradient_states obj=topology, n=97                 # rainbow, splits states internally
    gradient_protein_states obj=topology, n=97         # old name -> rainbow preset
    gradient_protein_split  prefix=1aay_state, n=97    # old name -> pre-split objects ({p}_001)
    spectrum_states obj=topology, n=97                 # old name -> blue->purple->pink preset
"""
from pymol import cmd

# Named state-gradient palettes (hex anchors, interpolated across N states).
RAINBOW9 = ("0xff7979", "0xffaa7f", "0xfff082", "0x88ff80", "0x7af4ff",
            "0x799fff", "0xc87fff", "0xff82f7", "0xff82aa")
BLUE_PURPLE_PINK = ("0x0000FF", "0x800080", "0xFF66B2")


def hex_to_rgb(h):
    """'0xRRGGBB' / '#RRGGBB' / 'RRGGBB' -> (r, g, b) floats in [0, 1]."""
    if h.startswith(("0x", "0X")):
        h = h[2:]
    elif h.startswith("#"):
        h = h[1:]
    return tuple(int(h[i:i + 2], 16) / 255 for i in (0, 2, 4))


def _interp(anchors, t):
    """Piecewise-linear interpolation across a list of RGB anchors at t in [0, 1]."""
    seg = t * (len(anchors) - 1)
    i = min(int(seg), len(anchors) - 2)
    f = seg - i
    a, b = anchors[i], anchors[i + 1]
    return [a[k] + (b[k] - a[k]) * f for k in range(3)]


def gradient_states(obj="topology", n=97, colors=RAINBOW9,
                    split=True, name_fmt="{obj}_{s:04d}", sele="polymer.protein"):
    """Colour per-state objects along a gradient by conformation index.

    split=True  -> run cmd.split_states(obj) first; objects are named {obj}_0001…
                   (the old gradient_protein_states / spectrum_states behaviour).
    split=False -> assume the per-state objects already exist; pass a matching
                   name_fmt (e.g. "{obj}_{s:03d}") — the old gradient_protein_split case.
    sele restricts the colouring (default: the protein only).
    """
    n = int(n)
    if split:
        cmd.split_states(obj)
        cmd.disable(obj)
    anchors = [hex_to_rgb(c) for c in colors]
    for s in range(1, n + 1):
        rgb = _interp(anchors, (s - 1) / (n - 1) if n > 1 else 0.0)
        cname = f"grad_{s}"
        cmd.set_color(cname, rgb)
        cmd.color(cname, f"{name_fmt.format(obj=obj, s=s)} and {sele}")


# Back-compat PyMOL command names, so existing sessions/notebooks keep working.
cmd.extend("gradient_states", gradient_states)
cmd.extend("gradient_protein_states",
           lambda obj="topology", n=97: gradient_states(obj, n, RAINBOW9))
cmd.extend("gradient_protein_split",
           lambda prefix="1aay_state", n=97, sele="polymer.protein":
           gradient_states(prefix, n, RAINBOW9, split=False,
                           name_fmt="{obj}_{s:03d}", sele=sele))
cmd.extend("spectrum_states",
           lambda obj="topology", n=97: gradient_states(obj, n, BLUE_PURPLE_PINK))
