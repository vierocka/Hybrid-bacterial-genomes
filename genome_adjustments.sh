#!/usr/bin/env bash
# genome_adjustments.sh
# Apply point mutations (SNPs) from a filtered VCF to a single-chromosome
# reference FASTA to produce an adjusted hybrid genome.
#
# Usage:
#   genome_adjustments.sh -r reference.fa -v filtered.vcf  -o hybrid.fa
#   genome_adjustments.sh -r reference.fa -s snps.tsv      -o hybrid.fa [-l check.log]
#
# -r  Reference FASTA (single chromosome, any line width)
# -v  Filtered VCF file  (uses columns 2=POS, 5=ALT; skips indels automatically)
# -s  Alternative SNP input: two-column TSV, no header — POS <tab> ALT
# -o  Output FASTA
# -l  Optional log file: lists every substitution applied (POS, REF_NT, ALT_NT)
#
# Requires: gawk (GNU awk), available by default on Linux; on macOS install via brew.

set -euo pipefail

usage() {
    grep "^#" "$0" | grep -v "^#!/" | sed 's/^# \{0,1\}//'
    exit 1
}

REF="" VCF="" SNPS="" OUT="" LOG=""
while [[ $# -gt 0 ]]; do
    case $1 in
        -r) REF="$2";  shift 2 ;;
        -v) VCF="$2";  shift 2 ;;
        -s) SNPS="$2"; shift 2 ;;
        -o) OUT="$2";  shift 2 ;;
        -l) LOG="$2";  shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

[[ -z "$REF" || -z "$OUT" ]]          && { echo "ERROR: -r and -o are required."; usage; }
[[ -z "$VCF" && -z "$SNPS" ]]         && { echo "ERROR: -v or -s is required.";   usage; }
[[ -n "$VCF"  && ! -f "$VCF"  ]]      && { echo "ERROR: VCF not found: $VCF";     exit 1; }
[[ -n "$SNPS" && ! -f "$SNPS" ]]      && { echo "ERROR: SNP file not found: $SNPS"; exit 1; }
[[ ! -f "$REF" ]]                      && { echo "ERROR: Reference not found: $REF"; exit 1; }

# Build a sorted (by position) SNP table: POS <tab> ALT, one SNP per line.
# VCF: skip header lines, keep only single-nucleotide ALTs (no indels).
# TSV: skip comment lines.
tmp_snp=$(mktemp)
trap "rm -f '$tmp_snp'" EXIT

if [[ -n "$VCF" ]]; then
    awk '!/^#/ && $5 ~ /^[ACGTacgt]$/ { print $2"\t"$5 }' "$VCF" \
        | sort -k1,1n > "$tmp_snp"
else
    awk '!/^#/ && $2 ~ /^[ACGTacgt]$/ { print $1"\t"$2 }' "$SNPS" \
        | sort -k1,1n > "$tmp_snp"
fi

n_snps=$(wc -l < "$tmp_snp")
echo "Reference: $REF"
echo "SNPs to apply: $n_snps"

# Apply substitutions with gawk.
# Strategy: read SNPs into an array, read the genome into one string, then
# walk the sorted SNP positions and emit chunks of unchanged sequence between
# them — O(genome_size + n_snps*log(n_snps)), no per-SNP file scanning.
gawk -v log_file="$LOG" '
FNR == NR {
    # Pass 1: load SNPs
    pos = $1 + 0
    snp[pos] = $2
    next
}
# Pass 2: read reference FASTA
/^>/ { header = $0; next }
      { seq = seq $0  }

END {
    # Sort SNP positions numerically
    n = 0
    for (p in snp) pos_arr[++n] = p + 0
    n = asort(pos_arr)   # gawk built-in numeric sort

    # Write log header
    if (log_file != "") print "position\tref_nt\talt_nt" > log_file

    # Walk sorted positions, emitting chunks of unchanged sequence between SNPs
    print header
    out = ""
    prev = 1
    for (i = 1; i <= n; i++) {
        p   = pos_arr[i]
        ref_nt = substr(seq, p, 1)
        out = out substr(seq, prev, p - prev) snp[p]
        if (log_file != "")
            print p"\t"ref_nt"\t"snp[p] >> log_file
        prev = p + 1
    }
    out = out substr(seq, prev)

    # Print with 60-char line wrapping
    len = length(out)
    for (i = 1; i <= len; i += 60)
        print substr(out, i, 60)
}
' "$tmp_snp" "$REF" > "$OUT"

echo "Output written: $OUT"
[[ -n "$LOG" ]] && echo "Substitution log: $LOG"
