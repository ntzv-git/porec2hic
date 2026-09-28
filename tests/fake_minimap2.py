#!/usr/bin/env python3
"""Minimal minimap2 stand-in (mappy) used by the test: prints query-side PAF columns."""
import os, random, sys, mappy
# FAKE_CHAIN_JITTER=n : shrink alignment ends by 0..n bp to mimic chain-only (no -c) PAF coordinates
JIT = int(os.environ.get("FAKE_CHAIN_JITTER", 0))
random.seed(0)
args = sys.argv[1:]
opt = {}
pos = []
i = 0
while i < len(args):
    a = args[i]
    if a in ("-x", "-t", "-I", "--split-prefix", "-N"):
        opt[a] = args[i + 1]; i += 2
    elif a.startswith("-"):
        i += 1
    else:
        pos.append(a); i += 1
target, query = pos
al = mappy.Aligner(target, preset=opt.get("-x", "map-ont"), best_n=int(opt.get("-N", 5)) + 1)
for name, seq, _ in mappy.fastx_read(query):
    for h in al.map(seq):
        qs, qe = h.q_st + random.randint(0, JIT), h.q_en - random.randint(0, JIT)
        print(f"{name}\t{len(seq)}\t{qs}\t{qe}\t{'+' if h.strand > 0 else '-'}\t{h.ctg}\t{h.ctg_len}\t{h.r_st}\t{h.r_en}\t{h.mlen}\t{h.blen}\t{h.mapq}\ttp:A:{'P' if h.is_primary else 'S'}")
