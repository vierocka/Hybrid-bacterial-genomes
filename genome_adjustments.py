#!/usr/bin/env python3
"""
genome_adjustments.py

Apply point mutations (SNPs) from a filtered VCF to a single-chromosome
reference FASTA to produce an adjusted hybrid genome.

Usage:
    python genome_adjustments.py -r reference.fa -v filtered.vcf -o hybrid.fa
    python genome_adjustments.py -r reference.fa -s snps.tsv     -o hybrid.fa [--log check.log]

Arguments:
    -r / --reference   Reference FASTA (single chromosome, any line width)
    -v / --vcf         Filtered VCF (uses columns POS and ALT; indels skipped)
    -s / --snps        Alternative: two-column TSV (no header) — POS <tab> ALT
    -o / --output      Output FASTA
    --log              Optional log file: POS, REF_NT, ALT_NT for every substitution

Requires Python 3.6+, no external libraries.
"""

import argparse
import re
import sys


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("-r", "--reference", required=True,
                   help="Reference FASTA (single chromosome)")
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("-v", "--vcf",  help="Filtered VCF file")
    src.add_argument("-s", "--snps", help="Two-column TSV: POS<tab>ALT")
    p.add_argument("-o", "--output", required=True, help="Output FASTA")
    p.add_argument("--log", help="Write substitution log here")
    return p.parse_args()


def read_fasta(path):
    """Return (header, sequence) for a single-chromosome FASTA."""
    header, parts = "", []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                header = line
            else:
                parts.append(line)
    if not header:
        sys.exit(f"ERROR: no FASTA header found in {path}")
    return header, "".join(parts)


def read_snps_vcf(path):
    """Read SNPs from a VCF file. Returns dict {1-based pos: alt_nt}."""
    snps = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 5:
                continue
            pos, alt = int(fields[1]), fields[4].upper()
            if re.fullmatch(r"[ACGT]", alt):   # skip indels and multi-allelic
                snps[pos] = alt
    return snps


def read_snps_tsv(path):
    """Read SNPs from a two-column TSV (POS<tab>ALT). Returns dict {pos: alt}."""
    snps = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip().split("\t")
            pos, alt = int(fields[0]), fields[1].upper()
            if re.fullmatch(r"[ACGT]", alt):
                snps[pos] = alt
    return snps


def apply_snps(seq, snps, log_path=None):
    """
    Apply SNP substitutions to the sequence string.
    VCF positions are 1-based; Python strings are 0-based.
    Returns the modified sequence.
    """
    seq = bytearray(seq.encode())   # mutable, O(1) position access

    log_lines = ["position\tref_nt\talt_nt"] if log_path else None

    for pos in sorted(snps):
        idx = pos - 1
        if idx < 0 or idx >= len(seq):
            print(f"WARNING: position {pos} is out of range — skipped",
                  file=sys.stderr)
            continue
        ref_nt = chr(seq[idx])
        seq[idx] = ord(snps[pos])
        if log_lines is not None:
            log_lines.append(f"{pos}\t{ref_nt}\t{snps[pos]}")

    if log_path and log_lines:
        with open(log_path, "w") as f:
            f.write("\n".join(log_lines) + "\n")

    return seq.decode()


def write_fasta(path, header, seq, line_width=60):
    with open(path, "w") as f:
        f.write(header + "\n")
        for i in range(0, len(seq), line_width):
            f.write(seq[i : i + line_width] + "\n")


def main():
    args = parse_args()

    header, seq = read_fasta(args.reference)
    snps = read_snps_vcf(args.vcf) if args.vcf else read_snps_tsv(args.snps)

    print(f"Reference:     {len(seq):>12,} bp   ({args.reference})",  file=sys.stderr)
    print(f"SNPs to apply: {len(snps):>12,}",                         file=sys.stderr)

    new_seq = apply_snps(seq, snps, log_path=args.log)
    write_fasta(args.output, header, new_seq)

    print(f"Output written: {args.output}", file=sys.stderr)
    if args.log:
        print(f"Substitution log: {args.log}", file=sys.stderr)


if __name__ == "__main__":
    main()
