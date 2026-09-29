# porec2hic

Splits Pore-C reads into monomers at restriction sites, **keeping the sites that HiFi
reads show to be genomic**, then turns the monomers of each read into **all-to-all**
pseudo-Hi-C read pairs (R1/R2).

```bash
POREC_FQ=porec.fq.gz HIFI_FQ=hifi.fq.gz THREADS=96 ./porec2hic_hifi.sh
# -> porec2hic_out/porec_hic_R1.fastq.gz, porec_hic_R2.fastq.gz,
#    porec_hic_monomers.fastq.gz, sites.final.tsv.gz

# in the background, e.g. with a larger HiFi index batch:
nohup env POREC_FQ=porec.fq.gz HIFI_FQ=hifi.fq.gz MM2_BATCH=200G THREADS=96 \
    bash ./porec2hic_hifi.sh > porec2hic.log 2>&1 &
```

Requirements: `minimap2`, `seqkit`, `bedtools` (≥ 2.26), `gawk`, and optionally `pigz`.
No Python.

## Steps

| step | tools | output |
|---|---|---|
| 1. read lengths | `seqkit fx2tab -n -l` | `porec.genome` (FASTQ order) |
| 2. overlapping windows for long reads | `seqkit fx2tab` + `awk` | `porec_windows.fa.gz` |
| 3. Pore-C → HiFi alignment | `minimap2 -c -N` + `gawk` | `hifi_blocks.bed.gz`: continuous blocks of each HiFi read, in Pore-C coordinates |
| 4. restriction sites | `seqkit locate` + `awk` | exact cut position `x = motif start + CUT_OFFSET` |
| 5. HiFi protection | `bedtools map -f 1.0 -o count_distinct` | `sites.final.tsv.gz` (`span`, class P/C) |
| 6. cutting, pseudo-monomer filter, all-to-all pairs | `seqkit fx2tab` + `awk` | R1/R2, monomers |

- **No global sort:** all files stay in Pore-C FASTQ order, and
  `bedtools map -g porec.genome` works in that order.
- **Resume:** each step writes its output under a temporary name (`*.tmp`) and renames it
  only when the step succeeds. On a re-run, a step whose output file exists is skipped.
  As soon as one step runs again, all later steps run again too, so outputs stay
  consistent. To redo a step, delete its output file.

  | step | output checked |
  |---|---|
  | 1 | `porec.genome` |
  | 2 | `porec_windows.fa.gz` |
  | 3 | `hifi_blocks.bed.gz` |
  | 4 | `sites.tsv` (or `sites.final.tsv.gz`) |
  | 5 | `sites.final.tsv.gz` |
  | 6 | `porec_hic_R1.fastq.gz`, `porec_hic_R2.fastq.gz`, `porec_hic.stats` (and `porec_hic_monomers.fastq.gz`) |

- **Log:** every message is timestamped. The log gives the duration of each step and
  the total time at the end. The total covers the current run only: skipped steps
  count for 0.
- **Stop:** on `kill` or Ctrl-C, the running step and its child processes (minimap2…)
  are stopped at once, and the log says which step to resume.

## When is a HiFi read "continuous" across a site?

Pore-C reads (queries) are aligned on HiFi reads (targets). `-N` keeps the other HiFi
reads of the same locus. Each alignment is then turned into **continuous blocks**, in
Pore-C read coordinates:

1. the alignment is split at every indel ≥ `MAX_GAP`, read from the CIGAR (`-c`);
2. blocks of the **same HiFi read**, on the same strand, are merged when they are
   colinear: the hole between them is less than `MAX_GAP` on the Pore-C read **and** on
   the HiFi read, and their diagonals differ by less than `MAX_GAP`. This also joins the
   overlapping blocks of two consecutive windows (see below).

What this means for a HiFi read aligned on the bases flanking a site:

| situation | meaning | protects the site? |
|---|---|---|
| one block covers `[x-FLANK, x+FLANK]` | the HiFi read reads through the site: genomic motif | yes |
| two blocks of the same HiFi read, colinear, hole < `MAX_GAP` (e.g. a Pore-C sequencing error on the site) | merged, so continuous | yes |
| two blocks of the same HiFi read with a jump ≥ `MAX_GAP`, or on different strands (short-range cis ligation, inversion) | not continuous | no |
| blocks from different HiFi reads, one ending before `x` and the other starting after it | typical junction | no |
| a block that passes `x` by less than `FLANK` bp | alignment overshoot, not evidence of continuity | no |

## Short and long Pore-C reads: alignment by windows

minimap2 caps the number of secondary hits (`-N`) **per query**. With whole reads as
queries, all the monomers of a read share one cap:

- **Short reads:** a few monomers, so the cap is enough.
- **Long reads:** tens or hundreds of monomers. The longest monomers use up the cap, and
  the others keep only their primary alignment, so fewer than `MIN_COV` HiFi reads cross
  their sites, which are then cut by mistake.

Pore-C reads longer than `WINDOW` (2 kb) are therefore aligned as **overlapping
windows** of 2 kb, overlapping by `WINDOW_OVERLAP` (250 bp). The last window is aligned
on the read end. Reads of 2 kb or less stay whole.

- **Per-window cap.** Each window gets its own `-N` cap, so HiFi coverage per monomer no
  longer depends on read length.
- **Back to read coordinates.** Blocks are converted back to read coordinates, and the
  blocks of a HiFi read in two overlapping windows are merged.
- **Every site is evaluated.** Every site lies at least 125 bp from the edge of one
  window, so it is always evaluated away from window boundaries.
- **Cost.** Only the overlaps are aligned twice, about 14% extra bases for reads longer
  than 2 kb.

Simulation: 3,000 short reads (1–6 monomers) and 300 long reads (20–40 monomers,
about 25 kb), HiFi 20x.

| alignment | `MM2_N` | short reads: precision | long reads: precision | long reads: false cuts |
|---|---|---|---|---|
| whole reads | 100 | 0.903 | **0.472** | 9,333 |
| whole reads | 200 | 0.909 | 0.553 | 6,748 |
| windows 1 kb / 250 bp | 100 | 0.921 | 0.931 | 622 |
| **windows 2 kb / 250 bp (default)** | 200 | 0.916 | **0.925** | 679 |
| windows 4 kb / 500 bp | 100 | 0.907 | 0.854 | 1,439 |
| windows 2 kb / 500 bp | 30 | 0.729 | 0.744 | 2,898 |

- **Recall.** Recall on intact motifs is 99.4–99.6% in every configuration.
- **Window size.** 1–2 kb windows give the same result for short and long reads. 4 kb
  windows start to lose precision again.
- **`MM2_N`.** It must stay well above the HiFi depth: 30 is too low at 20x.
- **Recommendation for `MM2_N`.** Aim for at least 5 × the HiFi depth. With 36.9 Gb of
  HiFi, the default of 200 covers genomes of 1 Gb (37x) and larger.

### `MAX_GAP`

`MAX_GAP` (default 50 bp) is the size of a hole, meaning an insertion or deletion
between the Pore-C read and the HiFi read, from which the HiFi read is no longer
considered continuous.

- **Why it is needed.** Take a ligation between two fragments that are close in the
  genome, in the same orientation, for example 300 bp apart. minimap2 aligns it on a
  single HiFi read with a 300 bp deletion. Without `MAX_GAP`, that HiFi read would look
  continuous and would protect the junction.
- **Why 50 bp.** Nanopore and HiFi indel errors are much shorter than this, and 50 bp is
  the usual threshold for a structural variant.
- **Heterozygous structural variants.** A site next to a heterozygous SV of 50 bp or more
  is broken only for reads of the other haplotype. Reads from the same haplotype still
  protect it.

## Cutting rule

| class | condition | action |
|---|---|---|
| **P** protected | `span ≥ MIN_COV`: at least `MIN_COV` (3) **distinct** HiFi reads have a continuous block covering `[x-FLANK, x+FLANK]` entirely | not cut |
| **C** cut | `span < MIN_COV`, including when no HiFi read is present | **cut** |

At Pore-C read ends, the window is truncated to the available sequence.

This rule favours cutting: a site poorly covered by HiFi reads is cut. On simulated data
(400 kb genome, 20x HiFi, 3,000 concatemers with 1% errors):

- **Real junctions:** 99.4% of junctions with an intact motif are cut. Recall is 95.7%
  if junctions whose motif was altered by a sequencing error are counted too.
- **Genomic motifs cut by mistake:** 787, for 7,295 real junctions cut.
  - 552 lie less than 20 bp from a real junction (see *Pseudo-monomers*).
  - 235 lie elsewhere, on sites covered by only 0 to 2 HiFi reads.

## Choice of `FLANK`

At an NlaIII junction, the CATG motif belongs to **both** ligated fragments. A HiFi
alignment therefore runs past the junction by at least 4 bp, plus a few chance matches.
`FLANK` must be larger than this overshoot.

Overshoot of HiFi blocks past a real junction (simulation, minimap2 2.28 `-c`):

| overshoot (bp) | 1 | 2 | 3 | **4** | 5 | 6 | 7 | 8 | 9 | 10 | 11–15 | > 15 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| % of blocks | 14 | 5 | 2 | **56** | 14 | 4 | 1.2 | 1.2 | 0.8 | 0.5 | 1 | < 0.5 |

- **Too small:** real junctions get protected because of the overshoot. With `FLANK=5`,
  27% of junctions are missed.
- **Too large:** more motifs close to a junction become unprotected, which produces
  longer pseudo-monomers.

The default is 15 bp.

## Pseudo-monomers

A motif lying less than about `FLANK` bp from a real junction is not crossed by any HiFi
read with the required margin. It is therefore cut too, and a **pseudo-monomer** of a few
bp appears between the two cuts. Monomer lengths on the simulation:

| length (bp) | 0–4 | 5–9 | 10–14 | 15–19 | 20–24 | 25–29 | 30–49 |
|---|---|---|---|---|---|---|---|
| pseudo-monomers | 4 | 132 | 304 | 131 | 18 | 16 | 69 |
| real monomers | 0 | 24 | 63 | 76 | 43 | 35 | 185 |

- **Filter threshold.** Pseudo-monomers are at most about `FLANK` + motif length long
  (19 bp for CATG with `FLANK=15`). Monomers shorter than `MIN_MONO_LEN` (default
  `FLANK` + motif length + 1 = 20 bp) are therefore dropped **before** pairing.
- **Cost.** Real monomers of that size cannot be mapped uniquely anyway.
- **To disable:** `MIN_MONO_LEN=1`.

## Pairs

- **Cut position.** `CUT_OFFSET` is the cut position inside the motif: NlaIII `CATG^` = 4,
  DpnII `^GATC` = 0.
- **Motif on both sides.** With `DUP_MOTIF=1`, the motif rebuilt by ligation is kept on
  both monomers (`…CATG | CATG…`), so each monomer matches its genomic sequence exactly.
- **All-to-all pairs.** A read with n kept monomers gives n(n-1)/2 pairs, named
  `@read:i-j/1` and `@read:i-j/2`, where i and j are the monomer ranks in the read.
- **Orientation.** R2 is written in the same orientation as the Pore-C read, not reverse
  complemented. This does not matter for contact maps.

## Parameters

| variable | default | role |
|---|---|---|
| `MOTIF` / `CUT_OFFSET` | `CATG` / 4 | enzyme (palindromic motif, IUPAC codes allowed) |
| `MIN_COV` | 3 | continuous HiFi reads needed to protect a site; otherwise the site is cut |
| `FLANK` | 15 | continuity margin (bp) |
| `MAX_GAP` | 50 | indel/hole (bp) that breaks continuity |
| `MIN_MONO_LEN` | `FLANK` + motif + 1 (20) | shorter monomers are dropped before pairing |
| `DUP_MOTIF` | 1 | keep the motif on both monomers |
| `WRITE_MONOMERS` | 1 | write `porec_hic_monomers.fastq.gz` (used by the QC below) |
| `MM2_PRESET` | `map-ont` | `lr:hq` for ONT R10 Q20+ reads |
| `MM2_N` | 200 | minimap2 secondary hits, capped **per window**; aim for ≥ 5 × HiFi depth |
| `WINDOW` / `WINDOW_OVERLAP` | 2000 / 250 | reads longer than `WINDOW` are aligned as overlapping windows (`WINDOW=0`: whole reads) |
| `MM2_BATCH` | 50G | HiFi index batch size (`-I` with `--split-prefix`) |
| `MM2_EXTRA` | – | extra minimap2 options |

## Cost

- **Volume.** Step 3 produces about (number of Pore-C monomers) × (HiFi depth)
  alignments, with CIGAR.
- **HiFi depth.** Time and the size of `hifi_blocks.bed.gz` scale with HiFi depth (up to
  `MM2_N` hits per window).
- **Windows.** They add about 14% of aligned bases for reads longer than 2 kb.

## Recommended post-analysis checks and filters

### 1. HiFi support at the sites (`sites.final.tsv.gz`)

```bash
zcat porec2hic_out/sites.final.tsv.gz | cut -f6 | sort -n | uniq -c   # distribution of span
```

Expect a bimodal distribution:

- **span = 0:** real junctions.
- **Peak near the HiFi depth (≤ `MM2_N`):** genomic sites.
- **Many sites at 1 to 2 (< `MIN_COV`):** these are cut for lack of HiFi support, and most
  of them are false cuts. Add HiFi data, increase `MM2_N`, or lower `MIN_COV`.

### 2. False-cut rate

A false cut splits a genomic fragment into two monomers that map **contiguously** on the
genome. Map the monomers (`porec_hic_monomers.fastq.gz`, in read order) on the assembly
and count consecutive monomers of the same read that map end to end:

```bash
minimap2 -x map-ont -t 32 --secondary=no assembly.fa porec2hic_out/porec_hic_monomers.fastq.gz | \
awk -v tol=100 '
  $1 == q { next }                                   # best line of each monomer only
  { q = $1; n = split($1, a, "_"); id = substr($1, 1, length($1) - length(a[n-1]) - length(a[n]) - 2) }
  $12 < 10 { pid = ""; next }                        # MAPQ < 10: not usable
  {
    if (id == pid) {
      tested++
      if ($6 == pt && $5 == ps) {                    # same contig, same strand, end to end
        gap = ($5 == "+") ? ($8 - $3) - (pte + (pql - pqe)) : (pts - (pql - pqe)) - ($9 + $3)
        if (gap > -tol && gap < tol) contig++
      }
    }
    pid = id; pt = $6; ps = $5; pts = $8; pte = $9; pql = $2; pqe = $4
  }
  END { printf "consecutive monomers: %d, contiguous on the assembly: %d (%.2f%%)\n", tested, contig, 100 * contig / tested }'
```

- **Reading the result.** The contiguous fraction is an upper bound of the false-cut rate,
  because it also includes real re-ligations of adjacent fragments. On the simulation it
  is 1.6%.
- **If it is high:** increase HiFi coverage or `MM2_N`, or lower `MIN_COV`.

### 3. Missed cuts (junctions left inside monomers)

A missed junction leaves a monomer made of two genomic fragments. On the Pore-C read, it
shows up as a **hole**: a position inside a monomer that no HiFi block crosses, while
HiFi reads align on both sides. The check below finds these holes from the pipeline
outputs, without any new alignment. Run it inside the output directory:

```bash
F=15; M=30; W=100; MC=3        # FLANK, margin from monomer ends, side window, MIN_COV
# regions crossed (with FLANK margin) by at least one continuous HiFi block
zcat hifi_blocks.bed.gz | awk -v f=$F 'BEGIN{OFS="\t"} $3 - $2 > 2 * f { print $1, $2 + f, $3 - f }' \
  | bedtools merge -i - > qc_spanned.bed
# holes = parts of the reads crossed by no HiFi block (read ends excluded)
bedtools complement -i qc_spanned.bed -g porec.genome \
  | awk 'NR == FNR { L[$1] = $2; next } $2 > 0 && $3 < L[$1]' porec.genome - > qc_holes.bed
# holes strictly inside a kept monomer (M bp away from its ends)
zcat porec_hic_monomers.fastq.gz | awk -v m=$M 'NR % 4 == 1 { n = split(substr($1, 2), a, "_")
    r = substr($1, 2, length($1) - length(a[n-1]) - length(a[n]) - 3)
    if (a[n] - a[n-1] > 2 * m) print r "\t" a[n-1] + m "\t" a[n] - m }' > qc_interiors.bed
bedtools intersect -sorted -g porec.genome -f 1.0 -u -a qc_holes.bed -b qc_interiors.bed > qc_inside.bed
# keep holes with >= MC HiFi reads on BOTH sides (otherwise it is a HiFi coverage gap)
paste qc_inside.bed \
  <(awk -v w=$W 'BEGIN{OFS="\t"} { s = $2 - w; print $1, (s < 0 ? 0 : s), $2 }' qc_inside.bed \
      | bedtools map -g porec.genome -a - -b hifi_blocks.bed.gz -c 4 -o count_distinct | cut -f4) \
  <(awk -v w=$W 'BEGIN{OFS="\t"} { print $1, $3, $3 + w }' qc_inside.bed \
      | bedtools map -g porec.genome -a - -b hifi_blocks.bed.gz -c 4 -o count_distinct | cut -f4) \
  | awk -v mc=$MC '$4 >= mc && $5 >= mc' > qc_missed.bed
CUTS=$(zcat sites.final.tsv.gz | awk '$7 == "C"' | wc -l); MISSED=$(wc -l < qc_missed.bed)
awk -v c=$CUTS -v m=$MISSED 'BEGIN { printf "missed junctions: %d ; cuts: %d ; missed / (cuts + missed) = %.2f%%\n", m, c, 100 * m / (c + m) }'
```

Validation on the simulation (16,450 real junctions, 745 missed):

- **Specificity:** 594 holes found, of which 593 are real missed junctions.
- **Sensitivity:** 80% of missed junctions are found. The others are the protected
  junctions (0.5%, no hole by definition) and junctions within 30 bp of a monomer end.
- **Reading the result:** multiply the count by about 1.25 to estimate the real number of
  missed junctions.
- **Denominator:** it includes false cuts, so the percentage is slightly underestimated.

`qc_missed.bed` lists the positions of the missed junctions, and its width column is
informative. The check does not depend on the motif, so it also finds ligations that did
not happen at a restriction site.

Expected causes of missed cuts:

| cause | expected rate | why |
|---|---|---|
| motif altered by a sequencing error | ≈ 1 − (1 − e)⁴: **4–8%** for Q17–Q19 reads (e = 1.2–2%) | the cut is only made on an intact motif |
| junction protected by ≥ `MIN_COV` HiFi reads | ≈ 0.5% (simulation) | ligation of genomically adjacent fragments, alignment overshoot; can increase with homeologous copies (`-p 0.5`) |
| ligation outside a restriction site | unknown, measured by this check | no motif to cut |

### 4. Mapping the pairs

```bash
bwa mem -5SP -T0 -t 32 assembly.fa porec_hic_R1.fastq.gz porec_hic_R2.fastq.gz | \
  pairtools parse --min-mapq 30 --walks-policy 5unique --max-inter-align-gap 30 \
                  --chroms-path assembly.genome | pairtools sort -o porec_hic.pairs.gz
```

- **Do not run `pairtools dedup`.** Every monomer starts and ends at a restriction site,
  so independent molecules often give pairs with identical coordinates. Removing them as
  "duplicates" would discard real contacts. Pore-C has no PCR step anyway.
- **MAPQ.** Filter on MAPQ (≥ 30 for contact maps, ≥ 1–10 for scaffolding).
- **Self-ligation and neighbouring fragments.** These are not informative. Remove cis
  pairs closer than about 1 kb:

  ```bash
  pairtools select '(chrom1 != chrom2) or (abs(pos1 - pos2) >= 1000)'
  ```

### 5. Weight of highly fragmented reads

A read with n monomers contributes n(n-1)/2 pairs, so a few highly fragmented reads can
dominate. The read name carries the monomer ranks (`read:i-j`). Pairs can be filtered
before mapping, keeping R1 and R2 in sync. For example, to keep only direct ligations
(`j = i + 1`):

```bash
paste <(zcat porec_hic_R1.fastq.gz | paste - - - -) <(zcat porec_hic_R2.fastq.gz | paste - - - -) | \
awk -F'\t' '{ split($1, a, ":"); split(a[length(a)], b, /[-\/]/); if (b[2] == b[1] + 1) print }' | \
tee >(cut -f1-4 | tr '\t' '\n' | gzip > direct_R1.fastq.gz) | cut -f5-8 | tr '\t' '\n' | gzip > direct_R2.fastq.gz
```

- **Cap on fragmentation.** The same approach can cap the number of monomers per read
  (maximum `j`).
- **Distribution.** The `porec_hic.stats` file and the `monomers/read` line of the log give
  the fragmentation distribution.
