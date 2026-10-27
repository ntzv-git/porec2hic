# porec2hic

Splits Pore-C reads into monomers at restriction sites, then turns the monomers of each
read into pseudo-Hi-C read pairs (R1/R2), all-to-all by default, for Hi-C tools such as
hifiasm or YAHS.

| tool | which sites are cut | speed | needs |
|---|---|---|---|
| `porec2hic_digest.sh` | **every** motif | one streaming pass | `mawk`, `pigz` |
| `porec2hic_hifi.sh` | motifs **not** read through by HiFi reads (genomic motifs are kept) | aligns Pore-C on HiFi: hours to days | `minimap2`, `seqkit`, `bedtools` ≥ 2.26, `gawk`, `pigz` |

Both tools write `PREFIX_R1.fastq.gz`, `PREFIX_R2.fastq.gz` and `PREFIX.stats`, name pairs
`@read:i-j/1` and `/2` (i, j = monomer ranks in the read), and take their parameters as
environment variables (`--help` lists them). Their default output directories differ
(`porec2hic_digest_out`, `porec2hic_out`), so both can run from the same directory.

## Comparison with other digestion tools

| | `porec2hic_digest.sh` | `porec2hic_hifi.sh` | `pore-c-py digest` (wc-pore-c) | HiC-Pro `digest_genome.py` |
|---|---|---|---|---|
| what is digested | Pore-C reads | Pore-C reads | Pore-C reads | reference genome |
| input | FASTQ(.gz) | FASTQ(.gz) + HiFi reads | FASTQ(.gz) or unaligned BAM | genome FASTA |
| sites cut | every motif | motifs not read through by HiFi reads | every motif | every motif |
| enzymes | one motif (IUPAC allowed) | one motif (IUPAC allowed) | one enzyme, by Biopython name | one or more, by name or sequence |
| output | R1/R2 pairs (FASTQ.gz) | R1/R2 pairs (FASTQ.gz) | monomers (BAM) | restriction fragments (BED) |
| short monomer filter | yes (`MIN_MONO_LEN`) | yes (`MIN_MONO_LEN`) | no | – |
| highly fragmented reads | `PAIR_MODE=seq` limits pairs | – | `--max_monomers` drops them | – |
| modified bases (MM/ML tags) | lost | lost | kept on each monomer | – |
| pairing | all-to-all or adjacent | all-to-all | no: a second tool is needed | – |

- **pore-c-py digest.** It cuts at the same position (NlaIII `CATG^`), so with
  `DUP_MOTIF=0 MIN_MONO_LEN=1` `porec2hic_digest.sh` gives the same monomers. It processes reads in
  Python and writes a BAM, which a second pass must decode again to build the pairs.
  `porec2hic_digest.sh` does both steps in one parallel pass, with no intermediate file.
- **HiC-Pro digest_genome.py.** It digests the genome, not the reads. HiC-Pro uses its BED
  to assign aligned Hi-C reads to restriction fragments. It cannot split Pore-C reads.

## `porec2hic_digest.sh`: cut every site

```bash
POREC_FQ=porec.fq.gz THREADS=96 ./porec2hic_digest.sh       # NlaIII, all-to-all, monomers >= 30 bp
nohup env POREC_FQ=porec.fq.gz THREADS=96 bash ./porec2hic_digest.sh > digest.log 2>&1 &
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
  Output is compressed at `GZ_LEVEL=6` (pigz default). `GZ_LEVEL=1` is faster but
  gives about 30% larger files.
- **Cut.** NlaIII `CATG^` (`CUT_OFFSET=4`) leaves the motif at the end of the left
  monomer. `DUP_MOTIF=1` keeps it on both monomers.
- **`PAIR_MODE=seq`.** Writes only adjacent monomers (n-1 pairs instead of n(n-1)/2).
  In all-to-all mode, a few highly fragmented reads can dominate the output.

## `porec2hic_hifi.sh`: cut only the sites HiFi reads do not read through

```bash
POREC_FQ=porec.fq.gz HIFI_FQ=hifi.fq.gz THREADS=96 ./porec2hic_hifi.sh
nohup env POREC_FQ=porec.fq.gz HIFI_FQ=hifi.fq.gz MM2_MAX_OCC=1000 THREADS=96 \
    bash ./porec2hic_hifi.sh > porec2hic.log 2>&1 &
```

A motif in a Pore-C read is either a ligation junction or a genomic motif. Cutting every
motif splits genomic fragments into several monomers. HiFi reads tell the two apart:
they read through genomic motifs, never through junctions.

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
  alignments run a few bp past a junction. `FLANK` must exceed this overshoot.
- **Pseudo-monomers.** A motif within ~`FLANK` bp of a junction is cut too, leaving a
  fragment of a few bp. Monomers shorter than `MIN_MONO_LEN` (default `FLANK` + motif
  + 1 = 20 bp) are dropped before pairing.
- **Long reads as windows.** minimap2 caps secondary hits (`-N`) per query, so the
  monomers of a long read would share one cap and many of their sites would lack HiFi
  support. Reads longer than `WINDOW` (2 kb) are therefore aligned as windows
  overlapping by 250 bp, each with its own cap.
- **Missed cuts.** A junction whose motif carries a sequencing error is not cut, since
  only intact motifs are evaluated.

### Parameters

| variable | default | role |
|---|---|---|
| `MOTIF` / `CUT_OFFSET` | `CATG` / 4 | enzyme (palindromic, IUPAC allowed) |
| `MIN_COV` | 3 | continuous HiFi reads needed to protect a site |
| `FLANK` / `MAX_GAP` | 15 / 50 | continuity margin, indel that breaks continuity (bp) |
| `MIN_MONO_LEN` | 20 | shorter monomers are dropped before pairing |
| `DUP_MOTIF` | 1 | keep the motif on both monomers |
| `WRITE_MONOMERS` | 1 | write `porec_hic_monomers.fastq.gz` |
| `MM2_PRESET` | `map-ont` | `lr:hq` for ONT R10 Q20+ reads |
| `MM2_N` | 200 | secondary hits per window; keep ≥ 5 × HiFi depth |
| `WINDOW` / `WINDOW_OVERLAP` | 2000 / 250 | window alignment (`WINDOW=0`: whole reads) |
| `MM2_BATCH` / `MM2_SPLIT` | 50G / auto | HiFi index batch size, `--split-prefix` only if several batches |
| `MM2_MAX_OCC` | minimap2 default | `minimap2 -f`, see *Cost* |
| `MM2_EXTRA` | – | extra minimap2 options |

### Cost

Step 2 produces about (Pore-C monomers) × (HiFi depth) alignments, and is the longest
step. HiFi reads are a redundant reference, so minimap2's repeat threshold gets very
high (log line `mid_occ = ...`), especially on polyploids, and alignment slows down a
lot. `MM2_MAX_OCC=1000` caps it. Keep it well above HiFi depth × number of copies:
too low a value discards seeds from single-copy sequence.

### Checking the result

`sites.final.tsv.gz` gives, for each site, the number of HiFi reads that read through it
(`span`, column 6). Expect a peak at 0 (junctions) and a peak near the HiFi depth
(genomic sites). Many sites at 1 or 2 are mostly genomic sites cut by mistake: add HiFi
data, raise `MM2_N`, or lower `MIN_COV`.
