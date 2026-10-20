"""
pymol_common.py — shared PyMOL render settings for the conformation figures.

Consolidates a `base_quality()` that was copy-pasted into 5 scripts under
analysis/analyses/conformation/ with three different cartoon_sampling values
and an inconsistent stick_quality — i.e. figures that were meant to look
identical did not. Import this instead:

    import sys, os
    sys.path.insert(0, os.path.join(TFCONF, "analysis", "common"))
    from pymol_common import base_quality
    base_quality()                      # defaults
    base_quality(cartoon_sampling=8)    # make_S.py's historical value
"""


def base_quality(cartoon_sampling=10, stick_quality=15, bg="white"):
    """Apply the standard figure-quality render settings.

    cartoon_sampling: 10 everywhere except make_S.py, which used 8. Pass it
    explicitly to reproduce an older figure exactly.
    """
    from pymol import cmd
    cmd.set("cartoon_sampling", cartoon_sampling)
    cmd.set("cartoon_fancy_helices", 0)
    cmd.set("cartoon_flat_sheets", 1)
    cmd.set("transparency_mode", 1)
    cmd.set("ray_opaque_background", 0)
    cmd.set("antialias", 2)
    cmd.set("ray_shadows", 0)
    cmd.set("stick_quality", stick_quality)
    cmd.bg_color(bg)
