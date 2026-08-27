import os
import argparse
import pandas as pd
import pyreadr
from scipy.spatial import KDTree

parser = argparse.ArgumentParser()
parser.add_argument('--tissue', required=True, help='tissue for analysis')
parser.add_argument('--output_dir', required=True, help='base output directory')
parser.add_argument('--project', required=True,
                    help='model config subdir (e.g. EGRET or cis_MatrixeQTL_GBAT_transPCO_FDR_0.1)')
parser.add_argument('--gene_info_file_path', default='../data/GTEx_V8.txt.gz',
                    help='path to gene info file')
parser.add_argument('--pval_threshold', type=float, default=0.01,
                    help='p_trans_part threshold for selecting trans-component genes')
args = parser.parse_args()

tissue = args.tissue
output_dir = args.output_dir
project = args.project
gene_info_file_path = args.gene_info_file_path
pval_threshold = args.pval_threshold

sumstats_path = f"{output_dir}/results_sumstats/{tissue}/{project}.txt"
sumstats = pd.read_csv(sumstats_path, header=0, sep='\t')

column_mapping = {
    'r2_gw': 'gw',
    'r2_cis_part': 'cis',
    'r2_trans_part': 'trans'
}
print('here')
sumstats['r2_gw'] = sumstats['r2_gw'].fillna(0)
sumstats['r2_cis_part'] = sumstats['r2_cis_part'].fillna(0)
sumstats['r2_trans_part'] = sumstats['r2_trans_part'].fillna(0)

sumstats['best_model'] = sumstats[['r2_gw', 'r2_cis_part', 'r2_trans_part']].idxmax(axis=1).map(column_mapping)

sumstats['trans_component'] = (
    (sumstats['p_trans_part'] < pval_threshold) &
    ((sumstats['best_model'] == 'gw') | (sumstats['best_model'] == 'trans'))
)

tissue_genes_path = f"{output_dir}/expression_files/{tissue}_expression.txt.gz"
tissue_genes = pd.read_csv(tissue_genes_path, header=0, sep='\t')

gene_annot = pd.read_csv(gene_info_file_path, sep='\t')
has_name = 'name' in gene_annot.columns
annot_cols = ['geneId', '#chrom', 'chromStart', 'chromEnd']
if has_name:
    annot_cols.insert(1, 'name')
gene_annot = gene_annot[annot_cols]
gene_annot['#chrom'] = gene_annot['#chrom'].str.replace('chr', '').astype(str)
gene_annot = gene_annot[gene_annot['geneId'].isin(tissue_genes['gene_id'])]

gene_kdtrees = {}
for chrom in gene_annot['#chrom'].unique():
    positions = gene_annot[gene_annot['#chrom'] == chrom][['chromStart']].values
    gene_kdtrees[chrom] = KDTree(positions)

results = []

for _, row in sumstats[sumstats['trans_component']].iterrows():
    gene = row['gene']
    print(gene)
    best_model = row['best_model']
    gene_chr = gene_annot.loc[gene_annot['geneId'] == gene, '#chrom'].values[0]

    rdat_path = f"{output_dir}/xtune_fusion_models/{tissue}/{project}/{gene}.wgt.RDat"
    try:
        rdat = pyreadr.read_r(rdat_path)
    except Exception as e:
        print(f"Failed to read {rdat_path}: {e}")
        continue

    model = rdat['cv.performance'].idxmax(axis=1).iloc[0]
    model_weights = rdat['wgt.matrix'][model]
    model_weights = model_weights[model_weights != 0]

    snp_info = rdat['snps']
    snp_info.index = rdat['wgt.matrix'].index  # Align index with weights

    snp_chr = snp_info.loc[model_weights.index, 'V1'].astype(str)
    trans_weights = model_weights[snp_chr != gene_chr]

    for snp_id, weight in trans_weights.items():
        snp_row = snp_info.loc[snp_id]
        snp_chr = str(snp_row['V1'])
        snp_pos = snp_row['V4']

        if snp_chr in gene_kdtrees:
            dist, idx = gene_kdtrees[snp_chr].query([[snp_pos]])
            nearest_gene = gene_annot[(gene_annot['#chrom'] == snp_chr)].iloc[idx[0]]

            result = {
                'target_gene': gene,
                'snp': snp_id,
                'snp_chr': snp_chr,
                'snp_pos': snp_pos,
                'nearest_gene_id': nearest_gene['geneId'],
                'distance': dist[0],
                'weight': weight
            }
            if has_name:
                result['nearest_gene'] = nearest_gene['name']
            results.append(result)

trans_sources = pd.DataFrame(results)
print(trans_sources.head())

out_dir = f"{output_dir}/grn_analysis/trans_connections"
os.makedirs(out_dir, exist_ok=True)
output_path = f"{out_dir}/{tissue}_trans_connections.txt"
trans_sources.to_csv(output_path, sep="\t", index=False)
