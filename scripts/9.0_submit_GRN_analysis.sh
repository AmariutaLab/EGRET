#!/bin/bash
#SBATCH --job-name="GRN_driver"
#SBATCH --output="batch_submissions/GRN_driver.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=1G
#SBATCH --account=csd832
#SBATCH --export=ALL
#SBATCH -t 00:30:00
# --- 9.0_submit_GRN_analysis.sh --------------------------------------------
# Stage `grn`, dispatcher. The TWAS drivers finish as soon as they have submitted
# their per-trait jobs, so depending on them is not enough; this reads their child
# job IDs and submits 9.1 afterok on every one.
#
#   in : the TWAS stages' job-ID files
#   out: submits 9.1_run_GRN_analysis.sh
#   run: run_EGRET_workflow.sh, stage `grn` (opt-in)
# --------------------------------------------------------------------------
set -euo pipefail

scripts_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd || true)
[[ -r "${scripts_dir}/workflow_lib.sh" ]] || scripts_dir=${EGRET_SCRIPTS_DIR:-${SLURM_SUBMIT_DIR:-}}
[[ -r "${scripts_dir}/workflow_lib.sh" ]] || { echo "ERROR: cannot locate workflow_lib.sh" >&2; exit 1; }
export EGRET_SCRIPTS_DIR="$scripts_dir"
source "${scripts_dir}/workflow_lib.sh"

USAGE="Usage: 9.0_submit_GRN_analysis.sh --tissue NAME --output_dir DIR --egret_subdir NAME --gene_info FILE --eqtl_file FILE --gwas_sumstat_dir DIR --gwas_traits LIST --conda_env NAME --slurm_account NAME --slurm_partition NAME --slurm_log_dir DIR --job_id_file FILE [--twas_job_id_file FILE] [--cis_twas_job_id_file FILE] [--pval_threshold P]"
CLI_FLAGS=(--tissue --output_dir --egret_subdir --gene_info --eqtl_file --gwas_sumstat_dir --gwas_traits --conda_env --slurm_account --slurm_partition --slurm_log_dir --job_id_file --twas_job_id_file --cis_twas_job_id_file --pval_threshold)
CLI_VARS=(tissue output_dir egret_subdir gene_info eqtl_file gwas_sumstat_dir gwas_traits conda_env slurm_account slurm_partition slurm_log_dir job_id_file twas_job_id_file cis_twas_job_id_file pval_threshold)
CLI_REQUIRED=(tissue output_dir egret_subdir gene_info eqtl_file gwas_sumstat_dir gwas_traits conda_env slurm_account slurm_partition slurm_log_dir job_id_file)
twas_job_id_file=""; cis_twas_job_id_file=""; pval_threshold="0.01"
parse_cli "$@"

child_ids=()
for id_file in "$twas_job_id_file" "$cis_twas_job_id_file"; do
    [[ -n "$id_file" ]] || continue
    ids=$(read_job_ids "$id_file")
    IFS=: read -ra file_ids <<< "$ids"
    child_ids+=("${file_ids[@]}")
done

sbatch_args=(--parsable --job-name="${tissue}_grn" --account="$slurm_account"
    --partition="$slurm_partition" --output="${slurm_log_dir}/GRN_analysis.%j.%N.out")
if ((${#child_ids[@]})); then
    sbatch_args+=(--dependency="afterok:$(join_by_colon "${child_ids[@]}")")
fi

jid=$(submit_with_retry "${sbatch_args[@]}" "${scripts_dir}/9.1_run_GRN_analysis.sh" \
    --tissue "$tissue" --output_dir "$output_dir" --egret_subdir "$egret_subdir" \
    --gene_info "$gene_info" --eqtl_file "$eqtl_file" \
    --gwas_sumstat_dir "$gwas_sumstat_dir" --gwas_traits "$gwas_traits" \
    --conda_env "$conda_env" --pval_threshold "$pval_threshold")
[[ -n "$jid" ]] || cli_die "no job ID returned for GRN analysis"
write_job_ids "$job_id_file" "$jid"
