#!/usr/bin/env python3
"""
porec_hifi_split.py -- HiFi-guided in-silico digestion of Pore-C reads and
all-to-all conversion of the resulting monomers into pseudo-Hi-C read pairs.

Single streaming pass (replaces seqkit locate + samtools depth + bedtools
subtract + seqkit subseq + awk of the previous pipeline):

  1. for each Pore-C read, list the restriction sites (IUPAC motif, both
     strands when the motif is not palindromic) and place the cut exactly at
     the enzyme cut position (e.g. NlaIII CATG^ -> offset 4);
  2. classify each site from the Pore-C -> HiFi alignments (query intervals):
       span  = #alignments covering [x-FLANK, x+FLANK]  (HiFi read continuous
               across the site -> genomic CATG, NOT a ligation junction)
       covL  = #alignments overlapping [x-W, x), covR = #alignments overlapping [x, x+W)
       PROTECTED  : span >= MIN_COV
       JUNCTION   : span <  MIN_COV and max(covL, covR) >= MIN_SIDE_COV
                    (HiFi evidence next to the site, but no HiFi read crosses it)
       UNRESOLVED : otherwise (no HiFi information) -> cut only with --cut-unresolved
     JUNCTION sites closer than 2*FLANK are merged: only the one closest to the
     alignment breakpoint is cut (other motifs next to a junction also lack
     spanning reads);
  3. cut the read at JUNCTION sites (optionally duplicating the motif on both
     monomers, as it is shared by the two ligated fragments);
  4. drop monomers shorter than --min-monomer-len and emit all C(n,2)
     monomer pairs as R1/R2.

Inputs
  --fastq  Pore-C reads (FASTQ, plain or .gz), in the same order as used for minimap2
  --aln    compact alignment table (plain or .gz), one line per aligned Pore-C read:
           read_id <TAB> qlen <TAB> qstart1,qstart2,... <TAB> qend1,qend2,...
           (built from minimap2 PAF by porec2hic_hifi.sh; same read order as --fastq)

Parallelism: one reader process, N worker processes; each worker writes its own
gzip shards, concatenated at the end (concatenated gzip members = valid gzip),
so R1/R2 stay in sync and no FASTQ text travels back through IPC.
"""
import argparse
import bisect
import gzip
import os
import queue
import re
import shutil
import subprocess
import sys
import time
import multiprocessing as mp
from collections import Counter

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[CG]", "W": "[AT]", "K": "[GT]", "M": "[AC]",
    "B": "[CGT]", "D": "[AGT]", "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}
COMP = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")

STAT_KEYS = [
    "reads", "reads_with_aln", "bases",
    "sites_total", "sites_edge", "sites_protected", "sites_junction", "sites_junction_adjacent",
    "sites_unresolved",
    "cuts", "monomers_raw", "monomers_short", "monomers_kept",
    "reads_multi", "reads_single", "reads_empty", "reads_capped", "pairs",
]

SITE_STAT = {"P": "sites_protected", "J": "sites_junction", "j": "sites_junction_adjacent",
             "U": "sites_unresolved"}


def revcomp(s):
    return s.translate(COMP)[::-1]


def build_patterns(motif, cut_offset):
    motif = motif.upper()
    if any(c not in IUPAC for c in motif):
        sys.exit(f"[porec_hifi_split] invalid motif: {motif}")
    if not 0 <= cut_offset <= len(motif):
        sys.exit("[porec_hifi_split] --cut-offset must be within [0, len(motif)]")
    pats = []
    for m, off in ((motif, cut_offset), (revcomp(motif), len(motif) - cut_offset)):
        rx = "(?=" + "".join(IUPAC[c] for c in m) + ")"   # lookahead -> overlapping hits
        pats.append((re.compile(rx.encode()), off))
        if revcomp(motif) == motif:                        # palindrome: one strand is enough
            break
    return pats


# ----------------------------------------------------------------------------- worker
CFG = None


def _median(v):
    return v[len(v) // 2] if v else None


def classify_sites(seq, starts, ends, st):
    """Return the list of cuts (x, motif_start, motif_end) for one read."""
    c = CFG
    f, w, L = c["flank"], c["side_window"], len(seq)
    useq = seq.upper()
    sites = {}
    for rx, off in c["patterns"]:
        for mt in rx.finditer(useq):
            s = mt.start()
            sites.setdefault(s + off, s)
    if not sites:
        return []
    m = c["motif_len"]
    res = []                                            # [x, span, covL, covR, class]
    for x in sorted(sites):
        st["sites_total"] += 1
        if x < f or x > L - f:
            st["sites_edge"] += 1
            continue
        # starts/ends are sorted and every interval is >= 2*flank long, hence
        # #(s <= x-f and e >= x+f) = #(s <= x-f) - #(e < x+f)
        span = bisect.bisect_right(starts, x - f) - bisect.bisect_left(ends, x + f)
        # alignments overlapping [x-w, x) and [x, x+w)
        cov_l = bisect.bisect_left(starts, x) - bisect.bisect_right(ends, x - w)
        cov_r = bisect.bisect_left(starts, x + w) - bisect.bisect_right(ends, x)
        if span >= c["min_cov"]:
            cls = "P"
        elif max(cov_l, cov_r) >= c["min_side_cov"]:
            cls = "J"
        else:
            cls = "U"
        res.append([x, span, cov_l, cov_r, cls])

    # Motifs closer than 2*flank to a junction also lack spanning reads: keep one
    # site per cluster, the closest to the alignment breakpoint (midpoint between
    # the median end of the left alignments and the median start of the right ones).
    i, n = 0, len(res)
    while i < n:
        if res[i][4] != "J":
            i += 1
            continue
        j = i
        while j + 1 < n and res[j + 1][4] == "J" and res[j + 1][0] - res[j][0] < 2 * f:
            j += 1
        if j > i:
            x0, x1 = res[i][0], res[j][0]
            e_in = ends[bisect.bisect_left(ends, x0 - w):bisect.bisect_right(ends, x1 + f)]
            s_in = starts[bisect.bisect_left(starts, x0 - f):bisect.bisect_right(starts, x1 + w)]
            me, ms = _median(e_in), _median(s_in)
            b = (me + ms) / 2 if me is not None and ms is not None else (me if ms is None else ms)
            best = i if b is None else min(range(i, j + 1), key=lambda k: abs(res[k][0] - b))
            for k in range(i, j + 1):
                if k != best:
                    res[k][4] = "j"                     # junction-adjacent motif, not cut
        i = j + 1

    cuts, rows = [], c["sites_rows"]
    for x, span, cov_l, cov_r, cls in res:
        if cls == "J" or (cls == "U" and c["cut_unresolved"]):
            cuts.append((x, sites[x], sites[x] + m))
        st[SITE_STAT[cls]] += 1
        if rows is not None:
            rows.append(b"%s\t%d\t%d\t%d\t%d\t%s\n" % (c["_name"], x, span, cov_l, cov_r, cls.encode()))
    return cuts


def process_read(name, seq, qual, aln, st, out):
    c = CFG
    L = len(seq)
    st["reads"] += 1
    st["bases"] += L
    f2 = 2 * c["flank"]
    starts, ends = [], []
    if aln is not None:
        st["reads_with_aln"] += 1
        ss, ee = aln
        for s, e in zip(ss.split(b","), ee.split(b",")):
            s, e = int(s), int(e)
            if e - s >= f2:                 # required by the span formula (see classify_sites)
                starts.append(s)
                ends.append(e)
        starts.sort()
        ends.sort()
    c["_name"] = name
    cuts = classify_sites(seq, starts, ends, st)
    st["cuts"] += len(cuts)

    segs, beg = [], 0
    dup = c["dup_motif"]
    for x, ms, me in cuts:
        end = me if dup else x
        if end > beg:
            segs.append((beg, end))
        beg = ms if dup else x
    if L > beg:
        segs.append((beg, L))
    st["monomers_raw"] += len(segs)

    kept = [(s, e) for s, e in segs if e - s >= c["min_mono_len"]]
    st["monomers_short"] += len(segs) - len(kept)
    st["monomers_kept"] += len(kept)
    n = len(kept)
    c["hist"][n] += 1

    if out["mono"] is not None:
        w = out["mono"]
        for s, e in kept:
            w.write(b"@%s_%d_%d\n%s\n+\n%s\n" % (name, s, e, seq[s:e], qual[s:e]))

    if n == 0:
        st["reads_empty"] += 1
        return
    if n == 1:
        st["reads_single"] += 1
        return
    if c["max_monomers"] and n > c["max_monomers"]:
        st["reads_capped"] += 1
        return
    st["reads_multi"] += 1
    recs = [b"\n%s\n+\n%s\n" % (seq[s:e], qual[s:e]) for s, e in kept]
    r1, r2 = out["r1"], out["r2"]
    buf1, buf2 = [], []
    for i in range(n):
        ri = recs[i]
        for j in range(i + 1, n):
            tag = b"@%s:%d-%d" % (name, i + 1, j + 1)
            buf1.append(tag + b"/1" + ri)
            buf2.append(tag + b"/2" + recs[j])
    r1.write(b"".join(buf1))
    r2.write(b"".join(buf2))
    st["pairs"] += len(buf1)


def worker(wid, cfg, in_q, res_q):
    global CFG
    CFG = cfg
    cfg["hist"] = Counter()
    lvl, tmp = cfg["gz_level"], cfg["tmpdir"]

    def gz(tag):
        return gzip.open(os.path.join(tmp, f"{tag}.{wid:04d}.gz"), "wb", compresslevel=lvl)

    out = {"r1": gz("R1"), "r2": gz("R2"),
           "mono": gz("monomers") if cfg["monomers"] else None}
    site_w = gz("sites") if cfg["sites"] else None
    st = Counter({k: 0 for k in STAT_KEYS})
    try:
        while True:
            batch = in_q.get()
            if batch is None:
                break
            cfg["sites_rows"] = [] if site_w is not None else None
            for name, seq, qual, aln in batch:
                process_read(name, seq, qual, aln, st, out)
            if site_w is not None:
                site_w.write(b"".join(cfg["sites_rows"]))
        for w in (out["r1"], out["r2"], out["mono"], site_w):
            if w is not None:
                w.close()
        res_q.put((wid, dict(st), dict(cfg["hist"]), None))
    except Exception as exc:                     # report instead of dying silently
        import traceback
        res_q.put((wid, None, None, traceback.format_exc()))
        raise exc


# ----------------------------------------------------------------------------- reader
def open_in(path):
    if path.endswith(".gz"):
        if shutil.which("pigz"):
            p = subprocess.Popen(["pigz", "-dc", "-p", "4", path], stdout=subprocess.PIPE, bufsize=1 << 20)
            return p.stdout
        return gzip.open(path, "rb")
    return open(path, "rb", buffering=1 << 20)


def read_batches(fq_path, aln_path, batch_size):
    fq, al = open_in(fq_path), open_in(aln_path)

    def next_aln():
        line = al.readline()
        if not line:
            return None
        p = line.rstrip(b"\r\n").split(b"\t")
        if len(p) < 4:
            sys.exit(f"[porec_hifi_split] malformed alignment line: {line[:200]!r}")
        return p[0], p[2], p[3]

    pending = next_aln()
    batch = []
    while True:
        h = fq.readline()
        if not h:
            break
        if h[:1] != b"@":
            sys.exit(f"[porec_hifi_split] FASTQ parse error near: {h[:200]!r}")
        seq = fq.readline().rstrip(b"\r\n")
        fq.readline()
        qual = fq.readline().rstrip(b"\r\n")
        name = h[1:].split(None, 1)[0]
        aln = None
        if pending is not None and pending[0] == name:
            aln = (pending[1], pending[2])
            pending = next_aln()
        batch.append((name, seq, qual, aln))
        if len(batch) >= batch_size:
            yield batch
            batch = []
    if batch:
        yield batch
    if pending is not None:
        sys.exit("[porec_hifi_split] ERROR: alignment table not exhausted (read "
                 f"{pending[0].decode()!r} never matched) -> FASTQ and alignments are not in the same order")


def concat(parts, dest):
    with open(dest, "wb") as o:
        for p in parts:
            with open(p, "rb") as i:
                shutil.copyfileobj(i, o, 16 << 20)
            os.remove(p)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fastq", required=True)
    ap.add_argument("--aln", required=True)
    ap.add_argument("--prefix", required=True, help="output prefix (-> PREFIX_R1.fastq.gz, PREFIX_R2.fastq.gz, PREFIX.stats.tsv)")
    ap.add_argument("--motif", default="CATG")
    ap.add_argument("--cut-offset", type=int, default=4, help="cut position inside the motif (NlaIII CATG^ = 4, DpnII ^GATC = 0)")
    ap.add_argument("--min-cov", type=int, default=3, help="HiFi alignments spanning a site needed to protect it")
    ap.add_argument("--min-side-cov", type=int, default=None, help="HiFi alignments on one side needed to call a junction (default: --min-cov)")
    ap.add_argument("--flank", type=int, default=25, help="bp an alignment must extend on both sides of a site to span it")
    ap.add_argument("--side-window", type=int, default=100, help="window (bp) on each side of a site where HiFi alignments are counted for junction calling")
    ap.add_argument("--cut-unresolved", action="store_true", help="also cut sites without HiFi information (old behaviour)")
    ap.add_argument("--no-dup-motif", action="store_true", help="cut exactly at the cut offset instead of keeping the motif on both monomers")
    ap.add_argument("--min-monomer-len", type=int, default=50)
    ap.add_argument("--max-monomers", type=int, default=0, help="skip reads with more monomers than this for pairing (0 = no limit)")
    ap.add_argument("--monomers", action="store_true", help="also write PREFIX_monomers.fastq.gz")
    ap.add_argument("--sites", action="store_true", help="also write PREFIX_sites.tsv.gz (read, pos, span, covL, covR, class)")
    ap.add_argument("-t", "--threads", type=int, default=os.cpu_count())
    ap.add_argument("--batch-size", type=int, default=2000)
    ap.add_argument("--gzip-level", type=int, default=3)
    a = ap.parse_args()

    nworkers = max(1, a.threads - 1)
    tmpdir = a.prefix + ".shards"
    os.makedirs(tmpdir, exist_ok=True)
    cfg = dict(
        patterns=build_patterns(a.motif, a.cut_offset), motif_len=len(a.motif),
        min_cov=a.min_cov, min_side_cov=a.min_cov if a.min_side_cov is None else a.min_side_cov,
        flank=a.flank, side_window=a.side_window, cut_unresolved=a.cut_unresolved, dup_motif=not a.no_dup_motif,
        min_mono_len=a.min_monomer_len, max_monomers=a.max_monomers,
        monomers=a.monomers, sites=a.sites, gz_level=a.gzip_level, tmpdir=tmpdir,
    )

    ctx = mp.get_context("fork")
    in_q, res_q = ctx.Queue(maxsize=4 * nworkers), ctx.Queue()
    procs = [ctx.Process(target=worker, args=(i, cfg, in_q, res_q), daemon=True) for i in range(nworkers)]
    for p in procs:
        p.start()

    def put(item):
        while True:
            try:
                in_q.put(item, timeout=5)
                return
            except queue.Full:
                if not all(p.is_alive() for p in procs):
                    sys.exit("[porec_hifi_split] a worker died, aborting")

    t0 = time.time()
    for batch in read_batches(a.fastq, a.aln, a.batch_size):
        put(batch)
    for _ in procs:
        put(None)

    tot, hist = Counter(), Counter()
    for _ in procs:
        wid, st, h, err = res_q.get()
        if err:
            sys.exit(f"[porec_hifi_split] worker {wid} failed:\n{err}")
        tot.update(st)
        hist.update(h)
    for p in procs:
        p.join()

    ids = range(nworkers)
    outs = [("R1", f"{a.prefix}_R1.fastq.gz"), ("R2", f"{a.prefix}_R2.fastq.gz")]
    if a.monomers:
        outs.append(("monomers", f"{a.prefix}_monomers.fastq.gz"))
    if a.sites:
        outs.append(("sites", f"{a.prefix}_sites.tsv.gz"))
    for tag, dest in outs:
        concat([os.path.join(tmpdir, f"{tag}.{i:04d}.gz") for i in ids], dest)
    os.rmdir(tmpdir)

    with open(f"{a.prefix}.stats.tsv", "w") as o:
        for k in STAT_KEYS:
            o.write(f"{k}\t{tot[k]}\n")
        for n in sorted(hist):
            o.write(f"monomers_per_read_{n}\t{hist[n]}\n")
    sys.stderr.write(f"[porec_hifi_split] {tot['reads']} reads -> {tot['pairs']} pairs in {time.time() - t0:.0f}s\n")


if __name__ == "__main__":
    main()
