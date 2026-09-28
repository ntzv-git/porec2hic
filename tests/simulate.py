#!/usr/bin/env python3
"""Simulate genome, HiFi reads and NlaIII Pore-C concatemers with known junctions."""
import random, sys, gzip, json
random.seed(1)
out = sys.argv[1]
G = 400_000
g = "".join(random.choice("ACGT") for _ in range(G))
COMP = str.maketrans("ACGT", "TGCA")
rc = lambda s: s.translate(COMP)[::-1]
cuts = [i + 4 for i in range(G - 3) if g.startswith("CATG", i)]   # NlaIII CATG^

def mut(s, rate):
    s = list(s)
    for i in range(len(s)):
        if random.random() < rate:
            s[i] = random.choice("ACGT".replace(s[i], ""))
    return "".join(s)

def fq(name, s):
    return f"@{name}\n{s}\n+\n{'5' * len(s)}\n"

with gzip.open(f"{out}/hifi.fq.gz", "wt") as o:           # 20x HiFi, 15 kb
    for k in range(G * 20 // 15000):
        p = random.randrange(0, G - 15000); s = g[p:p + 15000]
        o.write(fq(f"hifi{k}", mut(s if random.random() < .5 else rc(s), 0.001)))

truth = {}
with gzip.open(f"{out}/porec.fq.gz", "wt") as o:
    for k in range(3000):
        parts, junc = [], []
        for _ in range(random.randint(1, 6)):
            i = random.randrange(len(cuts) - 6)
            j = i + random.randint(1, 4)          # 1-4 NlaIII fragments (partial digestion -> genomic CATG inside)
            a, b = cuts[i], cuts[j]
            piece = g[a:b] if random.random() < .5 else rc(g[a - 4:b - 4])
            if parts:
                junc.append(sum(map(len, parts)))
            parts.append(piece)
        s = "".join(parts)
        truth[f"porec{k}"] = junc
        o.write(fq(f"porec{k}", mut(s, 0.01)))
json.dump(truth, open(f"{out}/truth.json", "w"))
