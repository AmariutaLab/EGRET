import os
import argparse
import pandas as pd
from convert_variant_id_to_rsid import convert_variant_id_to_rsid

parser = argparse.ArgumentParser()
parser.add_argument('--tissue', required=True, help='tissue for analysis')
parser.add_argument('--output_dir', required=True, help='base output directory')
parser.add_argument('--eqtl_dir', required=True,
                    help='directory containing GTEx v8 eQTL signif_variant_gene_pairs files')
parser.add_argument('--bim_file', required=True,
                    help='genotype .bim file used to map GTEx variant IDs to rsids')
args = parser.parse_args()

tissue = args.tissue
output_dir = args.output_dir
eqtl_dir = args.eqtl_dir
bim_file = args.bim_file

print(tissue)
convert_variant_id_to_rsid(tissue, eqtl_dir, bim_file)

trans_connections = pd.read_csv(
    f"{output_dir}/grn_analysis/trans_connections/{tissue}_trans_connections.txt",
    sep='\t'
)
eqtl_sumstats = pd.read_csv(
    f"{eqtl_dir}/{tissue}_with_rsid.txt",
    sep='\t'
)

merged = trans_connections.merge(
    eqtl_sumstats[['rsid', 'gene_id']],
    left_on='snp',
    right_on='rsid',
    how='left'
)

merged.rename(columns={'gene_id': 'egene'}, inplace=True)
merged.loc[merged['egene'].isna(), "egene"] = merged['nearest_gene_id']

out_dir = f"{output_dir}/grn_analysis/trans_connections_by_eqtl"
os.makedirs(out_dir, exist_ok=True)
merged.to_csv(f"{out_dir}/{tissue}_trans_connections_by_eqtl.txt", sep="\t", index=False)

print(len(merged['target_gene'].unique()))
print(len(merged['egene'].unique()))
print(len(merged['snp'].unique()))
