#!/bin/bash
# --- 2.0_run_cis_models.sh -------------------------------------------------
# Stage `cis`. Builds the per-gene cis windows, then splits the gene list into
# chunks and submits one 2.2 job per chunk.
#
#   in : gene_info, genotype_files/, individuals, covariates
#   out: plink_results/{tissue}/cis/{gene}.*; FUSION/{tissue}/cis/{gene}.wgt.RDat (via 2.2)
#   run: run_EGRET_workflow.sh, stage `cis`
# --------------------------------------------------------------------------
set -euo pipefail
# Locate workflow_lib.sh; see its header for why BASH_SOURCE alone is not enough.
scripts_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd || true)
[[ -r "${scripts_dir}/workflow_lib.sh" ]] || scripts_dir=${EGRET_SCRIPTS_DIR:-${SLURM_SUBMIT_DIR:-}}
[[ -r "${scripts_dir}/workflow_lib.sh" ]] || { echo "ERROR: cannot locate workflow_lib.sh (tried BASH_SOURCE, EGRET_SCRIPTS_DIR, SLURM_SUBMIT_DIR); export EGRET_SCRIPTS_DIR to the EGRET scripts directory" >&2; exit 1; }
export EGRET_SCRIPTS_DIR="$scripts_dir"
source "${scripts_dir}/workflow_lib.sh"
USAGE="Usage: 2.0_run_cis_models.sh --tissue NAME --output_dir DIR --plink_path PATH --gene_info FILE --individuals FILE --covariates FILE --plink1_path PATH --gcta_path PATH --gemma_path PATH --conda_env NAME --slurm_account NAME --slurm_partition NAME --slurm_log_dir DIR --job_id_file FILE [options] [--cis_window_size BP]"
CLI_FLAGS=(--tissue --output_dir --plink_path --gene_info --genotype_prefix --individuals --covariates --plink1_path --gcta_path --gemma_path --chunk_size --conda_env --slurm_account --slurm_partition --slurm_log_dir --job_id_file --cleanup_cis_working --cis_window_size --gene_list --hsq_set)
CLI_VARS=(tissue output_dir plink_path gene_info genotype_prefix individuals covariates plink1_path gcta_path gemma_path chunk_size conda_env slurm_account slurm_partition slurm_log_dir job_id_file cleanup_cis_working cis_window_size gene_list hsq_set)
CLI_REQUIRED=(tissue output_dir plink_path gene_info individuals covariates plink1_path gcta_path gemma_path conda_env slurm_account slurm_partition slurm_log_dir job_id_file)
genotype_prefix=GTEX_v8_genotypes_pruned; chunk_size=2000; cleanup_cis_working=FALSE
cis_window_size=1000000; gene_list=""; hsq_set=""
parse_cli "$@"; activate_workflow_env "$conda_env"
Rscript "${scripts_dir}/2.1_setup_cis_windows.R" --gene_info "$gene_info" --tissue "$tissue" \
    --plink_path "$plink_path" --output_dir "$output_dir" --genotype_prefix "$genotype_prefix" \
    --cis_window_size "$cis_window_size" ${gene_list:+--gene_list "$gene_list"}
cis_manifest="${output_dir}/plink_results/${tissue}/cis/cis_genes.txt"
require_file "$cis_manifest" "cis genotype manifest"
mapfile -t genes < "$cis_manifest"
for gene_file in "${genes[@]}"; do
    gene=${gene_file%.bed}
    for ext in bed bim fam; do
        require_file "${output_dir}/plink_results/${tissue}/cis/${gene}.${ext}" "cis genotype ${gene}.${ext}"
    done
done
num_chunks=$(((${#genes[@]} + chunk_size - 1) / chunk_size))
job_ids=()
for ((i=0; i<num_chunks; i++)); do
    start_index=$((i * chunk_size))
    jid=$(submit_with_retry --parsable --job-name="${tissue}_cis${i}" --account="$slurm_account" \
        --partition="$slurm_partition" --output="${slurm_log_dir}/cis_fusion.%j.%N.out" \
        "${scripts_dir}/2.2_compute_cis_weights_chunk.sh" --tissue "$tissue" --output_dir "$output_dir" \
        --chunk_size "$chunk_size" --start_index "$start_index" --plink_path "$plink_path" \
        --individuals "$individuals" --covariates "$covariates" --plink1_path "$plink1_path" \
        --gcta_path "$gcta_path" --gemma_path "$gemma_path" --conda_env "$conda_env" \
        --cleanup_cis_working "$cleanup_cis_working" ${hsq_set:+--hsq_set "$hsq_set"})
    [[ -n "$jid" ]] || cli_die "no job ID returned for cis chunk $i"
    job_ids+=("$jid")
done
((${#job_ids[@]})) || cli_die "no cis model jobs were submitted"
write_job_ids "$job_id_file" "${job_ids[@]}"
