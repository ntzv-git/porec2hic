#!/usr/bin/env python3
"""Précision/rappel des coupes par rapport à la vérité simulée (à lancer dans le dossier de test)."""
import collections, gzip, json, sys
out = sys.argv[1] if len(sys.argv) > 1 else "out"
truth = json.load(open("truth.json"))
seen = collections.defaultdict(dict)
for l in gzip.open(f"{out}/sites.final.tsv.gz", "rt"):
    r, ln, ms, me, x, sp, c = l.split()
    seen[r][int(x)] = c
cls_true, cls_gen = collections.Counter(), collections.Counter()
for r, js in truth.items():
    js = set(js)
    for x, c in seen[r].items():
        (cls_true if x in js else cls_gen)[c] += 1
    cls_true["absent (motif muté)"] += sum(j not in seen[r] for j in js)
print("vraies jonctions :", dict(cls_true))
print("sites génomiques :", dict(cls_gen))
tp, fp = cls_true["C"], cls_gen["C"]
print(f"précision={tp / (tp + fp):.4f} rappel={tp / sum(cls_true.values()):.4f}")
