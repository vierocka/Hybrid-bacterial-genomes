# Hybrid bacterial genome adjustment

Scripts to apply point mutations (SNPs) from a filtered VCF to a donor reference genome, producing a corrected reference for a bacterial hybrid strain. Used to improve RNA-seq read mappability in newly emerged bacterial hybrids.

This approach was used in:

> Murina V, Kasari M, Takada H, Hinnu M, Saha CK, Grimshaw JW, Seki T, Tozawa Y, Felden B, Hauryliuk V, Atkinson GC, Tenson T. **ABCF ATPases involved in protein synthesis, ribosome assembly and antibiotic resistance: structural and functional diversification across the tree of life.** *PNAS* 2021; 118(5): e2007873118. https://doi.org/10.1073/pnas.2007873118

---

## Problem

Two closely related bacterial strains (subspecies or species) have recombined to produce a hybrid. The donor strain has a well-annotated public reference genome; the hybrid does not. DNA-seq reads from the hybrid mapped against the donor reference reveal SNPs — nucleotide positions where the hybrid differs from the donor.

Using the donor reference as-is for RNA-seq mapping of the hybrid introduces systematic mismatches at those positions, reducing mapping precision and quantification accuracy. The corrected (hybrid-adjusted) reference eliminates this problem.

---

## Aim

Produce a corrected reference FASTA for a single bacterial chromosome by substituting donor nucleotides with hybrid-specific SNPs at their exact positions.

---

## Prerequisites

1. **Reference FASTA** — the donor genome from a public database (e.g. NCBI), single chromosome
2. **Filtered VCF** — SNPs identified between the donor and the hybrid

### How to obtain the VCF

```bash
# 1. Loose mapping of hybrid DNA-seq reads against the donor reference
#    (-O 1 -B 1: reduced gap-open and mismatch penalties for closely related strains)
bwa mem -O 1 -B 1 reference.fa hybrid_reads_R1.fq hybrid_reads_R2.fq \
    | samtools sort -o hybrid_vs_donor.bam
samtools index hybrid_vs_donor.bam

# 2. SNP calling with GATK
gatk HaplotypeCaller -R reference.fa -I hybrid_vs_donor.bam -O raw.vcf

# 3. SNP filtering (GATK hard filters, then custom filtering)
#    Keep only high-quality, single-nucleotide variants
gatk SelectVariants -R reference.fa -V raw.vcf --select-type-to-include SNP \
    -O snps_only.vcf
gatk VariantFiltration -R reference.fa -V snps_only.vcf \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
    --filter-name "basic_filter" -O filtered.vcf
```

---

## Scripts

Two implementations of the same algorithm — choose whichever fits your environment.

| Script | Language | Requires |
|--------|----------|----------|
| `genome_adjustments.sh` | Bash + gawk | gawk (standard on Linux; `brew install gawk` on macOS) |
| `genome_adjustments.py` | Python 3 | Python ≥ 3.6, no external libraries |

Both accept either a VCF file or a plain two-column TSV (position, alt allele) and produce an optional substitution log.

---

## Usage

### Bash / gawk

```bash
chmod +x genome_adjustments.sh

# From a VCF:
./genome_adjustments.sh \
    -r reference.fa \
    -v filtered.vcf \
    -o hybrid.fa

# From a two-column TSV (POS<tab>ALT, no header):
./genome_adjustments.sh \
    -r reference.fa \
    -s snps.tsv \
    -o hybrid.fa \
    -l substitutions.log
```

### Python

```bash
# From a VCF:
python genome_adjustments.py \
    -r reference.fa \
    -v filtered.vcf \
    -o hybrid.fa

# From a two-column TSV:
python genome_adjustments.py \
    -r reference.fa \
    -s snps.tsv \
    -o hybrid.fa \
    --log substitutions.log
```

---

## Inputs and outputs

**Reference FASTA** (`-r`): standard FASTA, single chromosome, any line width.

**VCF** (`-v`): standard VCF with at minimum columns POS (col 2) and ALT (col 5). Indels and multi-allelic sites are silently skipped — only single-nucleotide substitutions are applied.

**SNP TSV** (`-s`): alternative to VCF. Two columns, no header:
```
142     A
1073    C
4521    G
```

**Output FASTA** (`-o`): the adjusted hybrid reference, 60 bp per line, same header as input.

**Substitution log** (`-l` / `--log`): tab-separated, three columns — position, original nucleotide, substituted nucleotide. Useful for verifying the result.
```
position    ref_nt  alt_nt
142         T       A
1073        G       C
```

---

## How it works

Both scripts use the same algorithm:

1. Read all SNPs into a hash table (`position → alt_nucleotide`)
2. Load the reference sequence as a single string
3. Walk through SNP positions in sorted order, emitting unchanged chunks of sequence between substitution sites
4. Write the result as a 60-bp-per-line FASTA

This is O(genome\_size + n\_snps × log(n\_snps)) — a single pass through the genome regardless of SNP count, with no per-SNP file scanning.

The original bash implementation loaded the genome into a shell array (one element per nucleotide) and called `sed` once per SNP to look up each position — O(genome\_size × n\_snps). For a 4 Mb genome with 1,000 SNPs that is ~4 billion operations; the rewritten version performs ~5 million.

---

## Notes

- **Single chromosome only.** For multi-chromosome genomes, split the FASTA by chromosome and run the script once per chromosome.
- **VCF positions are 1-based.** The scripts handle the 1-to-0 index conversion internally.
- **No e-value cutoff was applied during SNP calling** in the original study — all GATK-filtered SNPs were retained. Evaluate your own VCF filtering stringency before applying.
- The FASTA header of the donor reference is preserved in the output. Rename it afterwards if needed (e.g. `sed -i 's/^>.*/\>hybrid_chromosome/' hybrid.fa`).
