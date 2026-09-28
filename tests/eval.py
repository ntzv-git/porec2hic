#!/usr/bin/env python3
"""Precision/recall of the junction calls against the simulated truth (run inside the test dir)."""
import gzip, json, collections, sys
truth = json.load(open("truth.json"))
tp=fp=fn=0; cls_true=collections.Counter(); cls_gen=collections.Counter(); spans=collections.Counter()
seen=collections.defaultdict(dict)
for l in gzip.open("out/porec_hic_sites.tsv.gz","rt"):
    r,x,sp,cl,cr,c=l.split(); seen[r][int(x)]=c
for r,js in truth.items():
    js=set(js)
    for x,c in seen[r].items():
        if x in js: cls_true[c]+=1
        else: cls_gen[c]+=1
    missing=[j for j in js if j not in seen[r]]
    cls_true["absent(edge/mutated)"]+=len(missing)
print("true junctions :", dict(cls_true)); print("genomic sites  :", dict(cls_gen))
J=cls_true["J"]; FP=cls_gen["J"]
print(f"precision={J/(J+FP):.4f} recall={J/sum(cls_true.values()):.4f}")
