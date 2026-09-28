# porec2hic

Splits Pore-C reads into monomers at restriction sites, **keeping the sites that HiFi
reads show to be genomic**, then turns the monomers of each read into **all-to-all**
pseudo-Hi-C read pairs (R1/R2).

```bash
POREC_FQ=porec.fq.gz HIFI_FQ=hifi.fq.gz THREADS=96 ./porec2hic_hifi.sh
# -> porec2hic_out/porec_hic_R1.fastq.gz, porec_hic_R2.fastq.gz,
#    porec_hic_monomers.fastq.gz, sites.final.tsv.gz
```

Requirements: `minimap2`, `seqkit`, `bedtools` (≥ 2.26), `gawk`, and optionally `pigz`.
No Python.

## Steps

| step | tools | output |
|---|---|---|
| 1. read lengths | `seqkit fx2tab -n -l` | `porec.genome` (FASTQ order) |
| 2. Pore-C → HiFi alignment | `minimap2 -c -N` + `gawk` | `hifi_blocks.bed.gz`: continuous blocks of each HiFi read, in Pore-C coordinates |
| 3. restriction sites | `seqkit locate` + `awk` | exact cut position `x = motif start + CUT_OFFSET` |
| 4. HiFi protection | `bedtools map -f 1.0 -o count_distinct` | `sites.final.tsv.gz` (`span`, class P/C) |
| 5. cutting, pseudo-monomer filter, all-to-all pairs | `seqkit fx2tab` + `awk` | R1/R2, monomers |

- **No global sort:** all files stay in Pore-C FASTQ order, and
  `bedtools map -g porec.genome` works in that order.
- **Resume:** each step writes a `stepN.done` file, so a re-run resumes where it stopped.

## When is a HiFi read "continuous" across a site?

Pore-C reads (queries) are aligned on HiFi reads (targets). `-N` keeps the other HiFi
reads of the same locus. Each alignment is then turned into **continuous blocks**, in
Pore-C read coordinates:

1. the alignment is split at every indel ≥ `MAX_GAP`, read from the CIGAR (`-c`);
2. blocks of the **same HiFi read**, on the same strand, are merged when they are
   colinear and separated by less than `MAX_GAP` on the Pore-C read **and** on the HiFi
   read.

What this means for a HiFi read aligned on the bases flanking a site:

| situation | meaning | protects the site? |
|---|---|---|
| one block covers `[x-FLANK, x+FLANK]` | the HiFi read reads through the site: genomic motif | yes |
| two blocks of the same HiFi read, colinear, hole < `MAX_GAP` (e.g. a Pore-C sequencing error on the site) | merged, so continuous | yes |
| two blocks of the same HiFi read with a jump ≥ `MAX_GAP`, or on different strands (short-range cis ligation, inversion) | not continuous | no |
| blocks from different HiFi reads, one ending before `x` and the other starting after it | typical junction | no |
| a block that passes `x` by less than `FLANK` bp | alignment overshoot, not evidence of continuity | no |

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
| `MM2_N` | 100 | minimap2 secondary hits. The cap applies **per Pore-C read**; aim for about `MIN_COV` × number of monomers × 5 |
| `MM2_BATCH` | 50G | HiFi index batch size (`-I` with `--split-prefix`) |
| `MM2_EXTRA` | – | extra minimap2 options |

## Cost

- **Volume.** Step 2 produces about (number of Pore-C monomers) × (HiFi depth)
  alignments, with CIGAR.
- **Subsampling HiFi.** With `MIN_COV=3`, 10–15x of HiFi is enough. Subsampling the HiFi
  reads (`seqkit sample`) reduces minimap2 time and the size of `hifi_blocks.bed.gz`
  accordingly.

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

### 3. Mapping the pairs

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

### 4. Weight of highly fragmented reads

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
- **Distribution.** The `step5.stats` file and the `monomers/read` line of the log give
  the fragmentation distribution.
