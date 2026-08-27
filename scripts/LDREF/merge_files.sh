#!/bin/bash
# Build the merged genome-wide LD reference that FUSION.assoc_test_trans.R needs.
#
# The repo ships 1000G.EUR.{1..22}.{bed,bim,fam} but not the merged fileset, because
# FUSION's trans association test calls read_plink() on a SINGLE bfile prefix -- per
# chromosome files cannot be used directly. Run this once after cloning:
#
#     bash LDREF/merge_files.sh /path/to/plink            # PLINK 1.9
#     PLINK1=/path/to/plink bash LDREF/merge_files.sh
#
# Then point config.sh's ld_ref at LDREF/1000G.EUR.merged (that is already the default).
#
# --keep-allele-order matters: without it PLINK 1.9 re-codes A1 to the minor allele, which
# would flip the sign of every TWAS Z whose weight was trained against the other coding.
set -euo pipefail

# Run from this script's own directory: the --merge-list below names the per-chromosome
# filesets by bare prefix, so a caller in any other CWD would fail to find them.
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

plink=${1:-${PLINK1:-plink}}
command -v "$plink" >/dev/null 2>&1 || [[ -x "$plink" ]] || {
    echo "ERROR: PLINK 1.9 not found at '${plink}'." >&2
    echo "       Pass it as the first argument or set PLINK1." >&2
    exit 1
}

"$plink" \
    --bfile 1000G.EUR.1 \
    --merge-list merge_list.txt \
    --make-bed --keep-allele-order \
    --out 1000G.EUR.merged

echo "merged LD reference written to $(pwd)/1000G.EUR.merged"
