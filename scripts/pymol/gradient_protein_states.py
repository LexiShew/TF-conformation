from pymol import cmd

def gradient_protein_states(obj="topology", n=97,
                             colors=("0xff7979", "0xffaa7f", "0xfff082", "0x88ff80", "0x7af4ff", "0x799fff", "0xc87fff", "0xff82f7", "0xff82aa")):
    n = int(n)

    cmd.split_states(obj)
    cmd.disable(obj)

    def hex_to_rgb(h):
        if h.startswith("0x") or h.startswith("0X"):
            h = h[2:]
        elif h.startswith("#"):
            h = h[1:]
        return tuple(int(h[i:i+2], 16)/255 for i in (0, 2, 4))

    anchors = [hex_to_rgb(c) for c in colors]

    def interp(t):
        seg = t * (len(anchors) - 1)
        i = min(int(seg), len(anchors) - 2)
        f = seg - i
        a, b = anchors[i], anchors[i+1]
        return [a[k] + (b[k]-a[k])*f for k in range(3)]

    for s in range(1, n + 1):
        rgb = interp((s-1)/(n-1))
        cname = f"grad_{s}"
        cmd.set_color(cname, list(rgb))
        cmd.color(cname, f"{obj}_{s:04d} and polymer.protein")

cmd.extend("gradient_protein_states", gradient_protein_states)