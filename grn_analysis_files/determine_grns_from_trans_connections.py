import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform

parser = argparse.ArgumentParser()
parser.add_argument('--tissue', required=True, help='tissue for analysis')
parser.add_argument('--output_dir', required=True, help='base output directory')
args = parser.parse_args()

tissue = args.tissue
output_dir = args.output_dir

DISTANCE_THRESHOLD = 0.9
MIN_CLUSTER_SIZE = 3
MIN_REGULATOR_TARGETS = 2

trans_connections = pd.read_csv(
    f"{output_dir}/grn_analysis/trans_connections_by_eqtl/{tissue}_trans_connections_by_eqtl.txt",
    sep="\t", header=0
)

trans_connections.loc[trans_connections['egene'].isna(), "egene"] = trans_connections['nearest_gene_id']
trans_connections = trans_connections[trans_connections['egene'].notna()]

regulators = trans_connections["egene"].unique()
targets = trans_connections["target_gene"].unique()

trans_connection_matrix = pd.DataFrame(
    0, index=regulators, columns=targets, dtype=int
)

for _, row in trans_connections.iterrows():
    trans_connection_matrix.loc[row["egene"], row["target_gene"]] = 1

G = trans_connection_matrix.T @ trans_connection_matrix

norm_sim = G.div(np.sqrt(np.outer(np.diag(G), np.diag(G))), axis=0)
dist_matrix = 1 - norm_sim.fillna(0)
np.fill_diagonal(dist_matrix.values, 0.0)

condensed_dist = squareform(dist_matrix, checks=False)
Z = linkage(condensed_dist, method="average")

output_dir_clusters = f"{output_dir}/grn_analysis/coregulation_grns/{tissue}"
os.makedirs(output_dir_clusters, exist_ok=True)

regulator_output_dir = f"{output_dir}/grn_analysis/coregulation_grn_regulators/{tissue}"
os.makedirs(regulator_output_dir, exist_ok=True)

plt.figure(figsize=(6, 4))
dendrogram(Z, labels=G.index, leaf_rotation=90)
plt.title("Hierarchical clustering of genes (based on shared regulators)")
plt.savefig(f"{output_dir_clusters}/grn_dendrogram.png", format='png')
plt.close()

clusters = fcluster(Z, t=DISTANCE_THRESHOLD, criterion='distance')
gene_clusters = pd.DataFrame({"Gene": G.index, "Cluster": clusters})
cluster_sizes = gene_clusters.groupby('Cluster').size()
multi_gene_clusters = cluster_sizes[cluster_sizes > MIN_CLUSTER_SIZE].index

for cluster_id in multi_gene_clusters:
    filtered_genes = gene_clusters[gene_clusters['Cluster'] == cluster_id]['Gene']

    submatrix = trans_connection_matrix.loc[:, filtered_genes]

    regulator_counts = submatrix.sum(axis=1).sort_values(ascending=False)
    top_regulators = regulator_counts[regulator_counts >= MIN_REGULATOR_TARGETS]

    out_path = os.path.join(output_dir_clusters, f"coregulation_cluster_{cluster_id}.txt")
    filtered_genes.to_csv(out_path, sep='\t', index=False, header=False)

    reg_path = os.path.join(regulator_output_dir, f"coregulation_cluster_{cluster_id}_regulators.txt")
    top_regulators.to_csv(reg_path, sep='\t', header=['num_targets'])

    print(f"{tissue} | Cluster {cluster_id}: {len(filtered_genes)} genes, {len(top_regulators)} key regulators")
