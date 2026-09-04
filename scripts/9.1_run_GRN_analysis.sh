#!/bin/bash
#SBATCH --job-name="GRN_analysis"
#SBATCH --output="batch_submissions/GRN_analysis.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32G
#SBATCH --account=csd832
#SBATCH --export=ALL
#SBATCH -t 08:00:00
# --- 9.1_run_GRN_analysis.sh -----------------------------------------------
# Runs the four-step GRN chain: identify trans regulators, assign them via GTEx eQTLs,
# cluster into co-regulation networks, then summarise TWAS Z per module.
# Step 4 is skipped, not failed, when --gwas_sumstat_dir is absent.
#
#   in : results_sumstats/{tissue}/{egret_subdir}.txt; eqtl_file; gene_info
#   out: grn_analysis/{trans_connections,trans_connections_by_eqtl,coregulation_grns}/
#   run: submitted by 9.0_submit_GRN_analysis.sh, or by hand
# --------------------------------------------------------------------------
# The headers above are not decoration: this cluster's bank_limit plugin rejects any job
# submitted without ntasks, so without them `sbatch 9.1_run_GRN_analysis.sh` fails outright.
set -euo pipefail
# Locate workflow_lib.sh; see its header for why BASH_SOURCE alone is not enough.
scripts_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd || true)
[[ -r "${scripts_dir}/workflow_lib.sh" ]] || scripts_dir=${EGRET_SCRIPTS_DIR:-${SLURM_SUBMIT_DIR:-}}
[[ -r "${scripts_dir}/workflow_lib.sh" ]] || { echo "ERROR: cannot locate workflow_lib.sh (tried BASH_SOURCE, EGRET_SCRIPTS_DIR, SLURM_SUBMIT_DIR); export EGRET_SCRIPTS_DIR to the EGRET scripts directory" >&2; exit 1; }
export EGRET_SCRIPTS_DIR="$scripts_dir"
source "${scripts_dir}/workflow_lib.sh"
USAGE="Usage: 9.1_run_GRN_analysis.sh --tissue NAME --output_dir DIR --egret_subdir NAME --gene_info FILE --eqtl_file FILE --conda_env NAME [--gwas_sumstat_dir DIR] [--gwas_traits LIST] [--grn_scripts_dir DIR] [--pval_threshold P]"
CLI_FLAGS=(--tissue --output_dir --egret_subdir --gene_info --eqtl_file --gwas_sumstat_dir --gwas_traits --grn_scripts_dir --conda_env --pval_threshold)
CLI_VARS=(tissue output_dir egret_subdir gene_info eqtl_file gwas_sumstat_dir gwas_traits grn_scripts_dir conda_env pval_threshold)
CLI_REQUIRED=(tissue output_dir egret_subdir gene_info eqtl_file conda_env)
grn_scripts_dir="$(cd "${scripts_dir}/.." && pwd)/grn_analysis_files"
gwas_sumstat_dir=""; gwas_traits=""; pval_threshold="0.01"
parse_cli "$@"
activate_workflow_env "$conda_env"

# This stage reads the cis/trans-partitioned sumstats emitted by the model-analysis stage.
require_file "${output_dir}/results_sumstats/${tissue}/${egret_subdir}.txt" "EGRET analysis summary"
require_file "$gene_info" "gene info file"
require_file "$eqtl_file" "eQTL summary statistics"

python3 "${grn_scripts_dir}/9.2_identify_trans_regulators.py" \
    --tissue "$tissue" \
    --output_dir "$output_dir" \
    --project "$egret_subdir" \
    --gene_info "$gene_info" \
    --pval_threshold "$pval_threshold"

python3 "${grn_scripts_dir}/9.3_assign_trans_regulators_with_eqtl.py" \
    --tissue "$tissue" \
    --output_dir "$output_dir" \
    --eqtl_file_path "$eqtl_file"

python3 "${grn_scripts_dir}/9.4_determine_grns.py" \
    --tissue "$tissue" \
    --output_dir "$output_dir"

# Step 4 compares TWAS Z between the EGRET and cis arms, so unlike the steps above it needs
# BOTH arms' results. Skipped rather than failed without --gwas_sumstat_dir, so the chain
# still runs end to end on a tree that has models but no TWAS.
if [[ -z "$gwas_sumstat_dir" ]]; then
    echo "Skipping module Z-score summary: --gwas_sumstat_dir not supplied"
else
    IFS=',' read -ra traits <<< "$gwas_traits"
    ((${#traits[@]})) || cli_die "--gwas_traits must list at least one trait when --gwas_sumstat_dir is supplied"
    for gwas_file in "${traits[@]}"; do
        trait=${gwas_file%.sumstats.gz}
        require_file "${output_dir}/TWAS_association_results/${tissue}/${egret_subdir}/${trait}.dat" "EGRET TWAS result for $trait"
        require_file "${output_dir}/TWAS_association_results/${tissue}/cis/${trait}.dat" "cis TWAS result for $trait"
    done
    # 9.5 defaults to a single grn_analysis/module_z_score_summary.txt. This chain runs one
    # tissue at a time, so without a per-tissue name each run silently overwrites the last.
    Rscript "${grn_scripts_dir}/9.5_average_z_score_module.R" \
        --tissues "$tissue" \
        --output_dir "$output_dir" \
        --egret_subdir "$egret_subdir" \
        --gene_info "$gene_info" \
        --gwas_sumstat_dir "$gwas_sumstat_dir" \
        --gwas_traits "$gwas_traits" \
        --pval_threshold "$pval_threshold" \
        --out "${output_dir}/grn_analysis/module_z_score_summary_${tissue}.txt"
fi
