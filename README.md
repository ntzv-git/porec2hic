# porec2hic

Splits Pore-C reads into monomers at restriction sites, then turns the monomers of each
read into pseudo-Hi-C read pairs (R1/R2), all-to-all by default.

| tool | which sites are cut | speed | needs |
|---|---|---|---|
| `porec2hic_digest.sh` | **every** motif (same as `pore-c-py digest`) | one streaming pass | `mawk`, `pigz` |
| `porec2hic_hifi.sh` | motifs **not** read through by HiFi reads (genomic motifs are kept) | aligns Pore-C on HiFi: hours to days | `minimap2`, `seqkit`, `bedtools` ≥ 2.26, `gawk`, `pigz` |

Both tools write `PREFIX_R1.fastq.gz`, `PREFIX_R2.fastq.gz` and `PREFIX.stats`, name pairs
`@read:i-j/1` and `/2` (i, j = monomer ranks in the read), and take their parameters as
environment variables. Their default output directories differ (`porec2hic_digest_out`,
`porec2hic_out`), so both can run from the same directory.

## `porec2hic_digest.sh`: cut every site

```bash
POREC_FQ=porec.fq.gz THREADS=96 ./porec2hic_digest.sh       # NlaIII, all-to-all, monomers >= 30 bp
./porec2hic_digest.sh --help                                # all parameters
```

`decompress | paste - - - - | split -n r/WORKERS --filter=mawk`:

- **Workers.** Reads are dealt round robin to `WORKERS` mawk processes (default
  `THREADS / 4`). Each worker cuts every motif and drops monomers shorter than
  `MIN_MONO_LEN`. It writes its pairs to its own gzip chunks, through its own `pigz`.
- **Concatenation.** The chunks are concatenated in the same order for R1 and R2, so pairs
  stay in sync. `KEEP_CHUNKS=1` skips this step and prints comma-separated chunk lists,
  which `hifiasm --h1/--h2` accept.
- **Speed.** A worker digests 60–90 MB/s of FASTQ. Decompressing the input (`rapidgzip`
  or `igzip` are used when installed) and compressing the output are the limiting steps.
- **Cut.** NlaIII `CATG^` (`CUT_OFFSET=4`) leaves the motif at the end of the left
  monomer. `DUP_MOTIF=1` keeps it on both monomers.
- **`PAIR_MODE=seq`.** Writes only adjacent monomers (n-1 pairs instead of n(n-1)/2).

HiC-Pro `digest_genome.py` digests a reference genome into a BED of fragments. It does
not split reads, so it cannot replace this step.

## `porec2hic_hifi.sh`: cut only the sites HiFi reads do not read through

```bash
POREC_FQ=porec.fq.gz HIFI_FQ=hifi.fq.gz THREADS=96 ./porec2hic_hifi.sh
# -> porec2hic_out/porec_hic_R1.fastq.gz, porec_hic_R2.fastq.gz,
#    porec_hic_monomers.fastq.gz, sites.final.tsv.gz
```

A CATG motif in a Pore-C read is either a ligation junction or a genomic motif.
Cutting every motif splits genomic fragments into several monomers. HiFi reads
tell the two apart: they read through genomic motifs, never through junctions.

| step | tools | output |
|---|---|---|
| 1. read lengths, windows for long reads | `seqkit fx2tab` + `awk` | `porec.genome`, `porec_windows.fa.gz` |
| 2. Pore-C → HiFi alignment | `minimap2 -c -N` + `gawk` | `hifi_blocks.bed.gz`: continuous HiFi blocks, in Pore-C coordinates |
| 3. restriction sites | `seqkit locate` + `awk` | `sites.tsv` |
| 4. HiFi protection | `bedtools map` | `sites.final.tsv.gz` (`span`, class P/C) |
| 5. cutting, short monomer filter, pairs | `seqkit fx2tab` + `awk` | R1/R2, monomers |

- **Order.** All files stay in Pore-C FASTQ order: no global sort.
- **Resume.** Each step writes `*.tmp` and renames it on success. A re-run skips steps
  whose output exists and reruns every step after the first one that runs. To redo a
  step, delete its output.
- **Log.** Timestamps and per-step durations, kept across runs in `timings.tsv`. Counts
  are cached in `counts.tsv`. Ctrl-C or `kill` stops the step and its child processes.

### Cutting rule

A site at cut position `x` is **protected** (not cut) when at least `MIN_COV` (3)
distinct HiFi reads have one continuous block covering `[x-FLANK, x+FLANK]`. Every other
site is **cut**, including sites without HiFi coverage.

- **Continuous block.** An alignment is split at every indel ≥ `MAX_GAP` (50 bp). Without
  this, a ligation of two fragments a few hundred bp apart would align as one HiFi read
  with a deletion and protect the junction. Blocks of the same HiFi read are merged
  again when colinear with a hole < `MAX_GAP`, which covers sequencing errors on the
  site.
- **`FLANK` (15 bp).** At an NlaIII junction, CATG belongs to both fragments, so HiFi
  alignments run about 4 bp past a junction. `FLANK` must exceed this overshoot
  (`FLANK=5` misses 27% of junctions).
- **Pseudo-monomers.** A motif within ~`FLANK` bp of a junction is cut too, leaving a
  fragment of a few bp. Monomers shorter than `MIN_MONO_LEN` (default `FLANK` + motif
  + 1 = 20 bp) are dropped before pairing.
- **Long reads as windows.** minimap2 caps secondary hits (`-N`) per query, so the
  monomers of a long read would share one cap. Reads longer than `WINDOW` (2 kb) are
  therefore aligned as windows overlapping by 250 bp (+14% aligned bases). This raises
  precision on 25 kb reads from 0.47 to 0.93.

On simulated data (20x HiFi, 1% errors), 99.4% of junctions with an intact motif are
cut. Counting junctions whose motif carries a sequencing error, recall is 95.7%.

### Parameters

| variable | default | role |
|---|---|---|
| `MOTIF` / `CUT_OFFSET` | `CATG` / 4 | enzyme (palindromic, IUPAC allowed) |
| `MIN_COV` | 3 | continuous HiFi reads needed to protect a site |
| `FLANK` / `MAX_GAP` | 15 / 50 | continuity margin, indel that breaks continuity (bp) |
| `MIN_MONO_LEN` | 20 | shorter monomers are dropped before pairing |
| `DUP_MOTIF` | 1 | keep the motif on both monomers |
| `WRITE_MONOMERS` | 1 | write `porec_hic_monomers.fastq.gz` (used by the checks below) |
| `MM2_PRESET` | `map-ont` | `lr:hq` for ONT R10 Q20+ reads |
| `MM2_N` | 200 | secondary hits per window; keep ≥ 5 × HiFi depth |
| `WINDOW` / `WINDOW_OVERLAP` | 2000 / 250 | window alignment (`WINDOW=0`: whole reads) |
| `MM2_BATCH` / `MM2_SPLIT` | 50G / auto | HiFi index batch size, `--split-prefix` only if several batches |
| `MM2_MAX_OCC` | minimap2 default | `minimap2 -f`, see *Cost* |
| `MM2_EXTRA` | – | extra minimap2 options |

### Cost

Step 2 produces about (Pore-C monomers) × (HiFi depth) alignments. HiFi reads are a
redundant reference, so minimap2's repeat threshold gets very high (log line
`mid_occ = ...`): 11,530 on a pentaploid at 18x per copy, where step 2 took 112 h.
`MM2_MAX_OCC=1000` caps it. Keep it well above HiFi depth × copies: on simulated data a
value of 25 lost precision, while 1000 matched the default.

### Checks

**HiFi support** (`sites.final.tsv.gz`, column 6): expect a peak at 0 (junctions) and a
peak near the HiFi depth (genomic sites). Many sites at 1–2 are mostly false cuts: add
HiFi data, raise `MM2_N` or lower `MIN_COV`.

```bash
zcat porec2hic_out/sites.final.tsv.gz | cut -f6 | sort -n | uniq -c
```

**False-cut rate.** A false cut gives two consecutive monomers that map end to end on the
assembly. The fraction below is an upper bound, since it includes re-ligations of
adjacent fragments (1.6% on the simulation).

```bash
minimap2 -x map-ont -t 32 --secondary=no assembly.fa porec2hic_out/porec_hic_monomers.fastq.gz | \
awk -v tol=100 '
  $1 == q { next }
  { q = $1; n = split($1, a, "_"); id = substr($1, 1, length($1) - length(a[n-1]) - length(a[n]) - 2) }
  $12 < 10 { pid = ""; next }
  {
    if (id == pid) {
      tested++
      if ($6 == pt && $5 == ps) {
        gap = ($5 == "+") ? ($8 - $3) - (pte + (pql - pqe)) : (pts - (pql - pqe)) - ($9 + $3)
        if (gap > -tol && gap < tol) contig++
      }
    }
    pid = id; pt = $6; ps = $5; pts = $8; pte = $9; pql = $2; pqe = $4
  }
  END { printf "consecutive monomers: %d, contiguous on the assembly: %d (%.2f%%)\n", tested, contig, 100 * contig / tested }'
```

**Missed cuts.** A junction left inside a monomer shows up as a position that no HiFi
block crosses, with HiFi reads on both sides. Run inside the output directory. The script
finds about 80% of missed junctions, with almost no false positives, so multiply the
count by about 1.25. The expected rate is mostly motifs hit by sequencing errors: about
1 − (1 − e)⁴, so 4–8% at Q17–Q19.

```bash
F=15; M=30; W=100; MC=3        # FLANK, margin from monomer ends, side window, MIN_COV
zcat hifi_blocks.bed.gz | awk -v f=$F 'BEGIN{OFS="\t"} $3 - $2 > 2 * f { print $1, $2 + f, $3 - f }' \
  | bedtools merge -i - > qc_spanned.bed
bedtools complement -i qc_spanned.bed -g porec.genome \
  | awk 'NR == FNR { L[$1] = $2; next } $2 > 0 && $3 < L[$1]' porec.genome - > qc_holes.bed
zcat porec_hic_monomers.fastq.gz | awk -v m=$M 'NR % 4 == 1 { n = split(substr($1, 2), a, "_")
    r = substr($1, 2, length($1) - length(a[n-1]) - length(a[n]) - 3)
    if (a[n] - a[n-1] > 2 * m) print r "\t" a[n-1] + m "\t" a[n] - m }' > qc_interiors.bed
bedtools intersect -sorted -g porec.genome -f 1.0 -u -a qc_holes.bed -b qc_interiors.bed > qc_inside.bed
paste qc_inside.bed \
  <(awk -v w=$W 'BEGIN{OFS="\t"} { s = $2 - w; print $1, (s < 0 ? 0 : s), $2 }' qc_inside.bed \
      | bedtools map -g porec.genome -a - -b hifi_blocks.bed.gz -c 4 -o count_distinct | cut -f4) \
  <(awk -v w=$W 'BEGIN{OFS="\t"} { print $1, $3, $3 + w }' qc_inside.bed \
      | bedtools map -g porec.genome -a - -b hifi_blocks.bed.gz -c 4 -o count_distinct | cut -f4) \
  | awk -v mc=$MC '$4 >= mc && $5 >= mc' > qc_missed.bed
CUTS=$(zcat sites.final.tsv.gz | awk '$7 == "C"' | wc -l); MISSED=$(wc -l < qc_missed.bed)
awk -v c=$CUTS -v m=$MISSED 'BEGIN { printf "missed junctions: %d ; cuts: %d ; missed / (cuts + missed) = %.2f%%\n", m, c, 100 * m / (c + m) }'
```

## Using the pairs (both tools)

**Mapping.**

```bash
bwa mem -5SP -T0 -t 32 assembly.fa porec_hic_R1.fastq.gz porec_hic_R2.fastq.gz | \
  pairtools parse --min-mapq 30 --walks-policy 5unique --max-inter-align-gap 30 \
                  --chroms-path assembly.genome | pairtools sort -o porec_hic.pairs.gz
```

- **No `pairtools dedup`.** Monomers start and end at restriction sites, so independent
  molecules often give identical coordinates. Pore-C has no PCR step anyway.
- **Short-range pairs.** Remove cis pairs closer than ~1 kb:
  `pairtools select '(chrom1 != chrom2) or (abs(pos1 - pos2) >= 1000)'`.

**Highly fragmented reads.** A read with n monomers gives n(n-1)/2 pairs, so a few reads
can dominate. To keep only direct ligations (`j = i + 1`), or to cap `j`, filter on the
name while keeping R1 and R2 in sync:

```bash
paste <(zcat porec_hic_R1.fastq.gz | paste - - - -) <(zcat porec_hic_R2.fastq.gz | paste - - - -) | \
awk -F'\t' '{ split($1, a, ":"); split(a[length(a)], b, /[-\/]/); if (b[2] == b[1] + 1) print }' | \
tee >(cut -f1-4 | tr '\t' '\n' | gzip > direct_R1.fastq.gz) | cut -f5-8 | tr '\t' '\n' | gzip > direct_R2.fastq.gz
```
