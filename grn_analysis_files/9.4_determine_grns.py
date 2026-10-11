# --- 9.4_determine_grns.py -------------------------------------------------
# GRN step 3. Builds a regulator x target incidence matrix, cosine-normalises it, and
# hierarchically clusters targets into co-regulation networks.
#
#   in : grn_analysis/trans_connections_by_eqtl/
#   out: grn_analysis/coregulation_grns/{tissue}/ and coregulation_grn_regulators/{tissue}/
#   run: 9.1_run_GRN_analysis.sh
# --------------------------------------------------------------------------
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
# Clustering knobs, deliberately CLI-only rather than config.sh entries: they tune the shape
# of the dendrogram, not the identity of the run. Override when exploring granularity.
parser.add_argument('--distance_threshold', type=float, default=0.9,
                    help='fcluster cophenetic-distance cutoff (default: 0.9)')
parser.add_argument('--min_cluster_size', type=int, default=4,
                    help='minimum target genes for a cluster to be reported (default: 4)')
parser.add_argument('--min_regulator_targets', type=int, default=2,
                    help='minimum targets for a regulator to be reported (default: 2)')
args = parser.parse_args()
tissue = args.tissue
output_dir = args.output_dir

DISTANCE_THRESHOLD = args.distance_threshold
MIN_CLUSTER_SIZE = args.min_cluster_size
MIN_REGULATOR_TARGETS = args.min_regulator_targets

trans_connections = pd.read_csv(
    f"{output_dir}/grn_analysis/trans_connections_by_eqtl/{tissue}_trans_connections_by_eqtl.txt",
    sep="\t", header=0
)

required_columns = {'target_gene', 'egene', 'nearest_gene_id'}
missing_columns = sorted(required_columns - set(trans_connections.columns))
if missing_columns:
    raise ValueError(
        "trans connections are missing required columns: " + ", ".join(missing_columns)
    )

output_dir_clusters = f"{output_dir}/grn_analysis/coregulation_grns/{tissue}"
os.makedirs(output_dir_clusters, exist_ok=True)
regulator_output_dir = f"{output_dir}/grn_analysis/coregulation_grn_regulators/{tissue}"
os.makedirs(regulator_output_dir, exist_ok=True)

# A rerun may produce fewer clusters. Remove only files owned by this stage so stale modules
# cannot leak into the downstream TWAS summary.
for directory, prefix in (
    (output_dir_clusters, "coregulation_cluster_"),
    (regulator_output_dir, "coregulation_cluster_"),
):
    for name in os.listdir(directory):
        if name.startswith(prefix) and name.endswith(".txt"):
            os.remove(os.path.join(directory, name))
old_dendrogram = os.path.join(output_dir_clusters, "grn_dendrogram.png")
if os.path.exists(old_dendrogram):
    os.remove(old_dendrogram)

trans_connections.loc[trans_connections['egene'].isna(), "egene"] = trans_connections['nearest_gene_id']
trans_connections = trans_connections[trans_connections['egene'].notna()]

regulators = trans_connections["egene"].unique()
targets = trans_connections["target_gene"].unique()

if len(regulators) == 0 or len(targets) < 2:
    print(
        f"{tissue}: no clusterable network "
        f"({len(regulators)} regulators, {len(targets)} targets); wrote no modules"
    )
    raise SystemExit(0)

trans_connection_matrix = pd.DataFrame(
    0, index=regulators, columns=targets, dtype=int
)

for _, row in trans_connections.iterrows():
    trans_connection_matrix.loc[row["egene"], row["target_gene"]] = 1

G = trans_connection_matrix.T @ trans_connection_matrix

norm_sim = G.div(np.sqrt(np.outer(np.diag(G), np.diag(G))), axis=0)
dist_matrix = 1 - norm_sim.fillna(0)
dist_values = dist_matrix.to_numpy(copy=True)
np.fill_diagonal(dist_values, 0.0)

condensed_dist = squareform(dist_values, checks=False)
Z = linkage(condensed_dist, method="average")

plt.figure(figsize=(6, 4))
dendrogram(Z, labels=G.index, leaf_rotation=90)
plt.title("Hierarchical clustering of genes (based on shared regulators)")
plt.savefig(f"{output_dir_clusters}/grn_dendrogram.png", format='png')
plt.close()

clusters = fcluster(Z, t=DISTANCE_THRESHOLD, criterion='distance')
gene_clusters = pd.DataFrame({"Gene": G.index, "Cluster": clusters})
cluster_sizes = gene_clusters.groupby('Cluster').size()
multi_gene_clusters = cluster_sizes[cluster_sizes >= MIN_CLUSTER_SIZE].index

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
