#!/usr/bin/env python3
"""seq_to_kdna.py — map a DNA duplex sequence to per-base-pair-step k_dna values
using the hexABC sequence-dependent stiffness prior.

Usage (as a library, from stage3_minimize.py):
    from seq_to_kdna import load_table, seq_to_kdna
    table = load_table("stiffness_table.json")
    k_per_step = seq_to_kdna(seq, table, k_base=1.5)   # len = len(seq)-1

The per-step k values are then broadcast onto the DNA backbone atoms of the two
nucleotides flanking each step, replacing the scalar k_dna in the Stage-3
CustomExternalForce. Steps that are floppier than average (TA, CA, TG) get a
weaker tether (more freedom to deform); stiffer steps (AT, GC) get a stronger one.
"""
import json
import numpy as np

_COMP = {"A":"T","T":"A","G":"C","C":"G"}

def load_table(path):
    t = json.load(open(path))
    return {s: d["mean"] for s, d in t["stiffness"]["sum"].items()}

def _step_key(a, b, table):
    di = a + b
    wc = _COMP[b] + _COMP[a]
    key = f"{di}/{wc}"
    return key if key in table else None

def seq_to_kdna(seq, table, k_base=1.5):
    """Return per-step k_dna array (length len(seq)-1), mean-normalized to k_base."""
    seq = seq.upper()
    mean_all = float(np.mean(list(table.values())))
    stiff = []
    for i in range(len(seq) - 1):
        key = _step_key(seq[i], seq[i+1], table)
        stiff.append(table.get(key, mean_all) if key else mean_all)
    stiff = np.asarray(stiff, float)
    return (stiff / stiff.mean()) * k_base
