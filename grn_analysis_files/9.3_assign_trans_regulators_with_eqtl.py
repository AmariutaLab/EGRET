# --- 9.3_assign_trans_regulators_with_eqtl.py ------------------------------
# GRN step 2. Replaces the nearest-gene guess with the gene each trans SNP is actually
# an eQTL for in GTEx, falling back to nearest gene where there is no eQTL.
#
#   in : grn_analysis/trans_connections/; eqtl_file (GTEx eQTLs with an rsid column)
#   out: grn_analysis/trans_connections_by_eqtl/{tissue}_trans_connections_by_eqtl.txt
#   run: 9.1_run_GRN_analysis.sh
# --------------------------------------------------------------------------
import os
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--tissue', required=True, help='tissue for analysis')
parser.add_argument('--output_dir', required=True, help='base output directory')
parser.add_argument('--eqtl_file_path', required=True,
                    help='file containg eqtl sumstats for specific tissue')
args = parser.parse_args()

tissue = args.tissue
output_dir = args.output_dir
eqtl_file_path = args.eqtl_file_path

trans_connections = pd.read_csv(
    f"{output_dir}/grn_analysis/trans_connections/{tissue}_trans_connections.txt",
    sep='\t'
)
eqtl_sumstats = pd.read_csv(
    f"{eqtl_file_path}",
    sep='\t'
)

def require_columns(frame, required, label):
    missing = sorted(set(required) - set(frame.columns))
    if missing:
        raise ValueError(f"{label} is missing required columns: {', '.join(missing)}")

require_columns(
    trans_connections,
    ['target_gene', 'snp', 'nearest_gene_id'],
    'trans connections'
)
require_columns(eqtl_sumstats, ['rsid', 'gene_id', 'pval_nominal'], 'eQTL summary statistics')

merged = trans_connections.merge(
    eqtl_sumstats[['rsid', 'gene_id', 'pval_nominal']],
    left_on='snp',
    right_on='rsid',
    how='left'
)

merged.rename(columns={'gene_id': 'egene'}, inplace=True)

# Separate matched and unmatched rsids
matched = merged[merged['rsid'].notna()]
unmatched = merged[merged['rsid'].isna()]

# Keep the most significant eGene association per target-gene/SNP pair.
matched = (
    matched
    .sort_values(
        ['target_gene', 'rsid', 'pval_nominal']
    )
    .groupby(
        ['target_gene', 'rsid'],
        group_keys=False
    )
    .head(1)
)

# Recombine
merged = pd.concat(
    [matched, unmatched],
    ignore_index=True
)

# Fill missing egene values
merged.loc[
    merged['egene'].isna(),
    'egene'
] = merged['nearest_gene_id']

out_dir = f"{output_dir}/grn_analysis/trans_connections_by_eqtl"
os.makedirs(out_dir, exist_ok=True)
merged.to_csv(f"{out_dir}/{tissue}_trans_connections_by_eqtl.txt", sep="\t", index=False)

print(
    f"Wrote {len(merged)} connections: "
    f"{merged['target_gene'].nunique()} targets, "
    f"{merged['egene'].nunique()} regulators, "
    f"{merged['snp'].nunique()} SNPs"
)
