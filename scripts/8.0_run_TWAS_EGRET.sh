#!/bin/bash
#SBATCH --job-name="TWAS_driver"
#SBATCH --output="batch_submissions/TWAS_driver.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8G
#SBATCH --account=csd832
#SBATCH --export=ALL
#SBATCH -t 04:00:00
# --- 8.0_run_TWAS_EGRET.sh -------------------------------------------------
# Stage `twas`. Builds the EGRET .pos file, then submits one 8.4 job per GWAS trait.
#
#   in : results_sumstats/{tissue}/{egret_subdir}.txt; gwas_sumstat_dir
#   out: pos_files/{tissue}/{egret_subdir}.pos; TWAS_association_results/{tissue}/{egret_subdir}/{trait}.dat
#   run: run_EGRET_workflow.sh, stage `twas`
# --------------------------------------------------------------------------
set -euo pipefail
# Locate workflow_lib.sh; see its header for why BASH_SOURCE alone is not enough.
scripts_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd || true)
[[ -r "${scripts_dir}/workflow_lib.sh" ]] || scripts_dir=${EGRET_SCRIPTS_DIR:-${SLURM_SUBMIT_DIR:-}}
[[ -r "${scripts_dir}/workflow_lib.sh" ]] || { echo "ERROR: cannot locate workflow_lib.sh (tried BASH_SOURCE, EGRET_SCRIPTS_DIR, SLURM_SUBMIT_DIR); export EGRET_SCRIPTS_DIR to the EGRET scripts directory" >&2; exit 1; }
export EGRET_SCRIPTS_DIR="$scripts_dir"
source "${scripts_dir}/workflow_lib.sh"
USAGE="Usage: 8.0_run_TWAS_EGRET.sh --tissue NAME --output_dir DIR --gene_info FILE --gwas_sumstat_dir DIR --ld_ref PREFIX --egret_subdir NAME --gwas_traits CSV --conda_env NAME --slurm_account NAME --slurm_partition NAME --slurm_log_dir DIR --job_id_file FILE"
CLI_FLAGS=(--tissue --output_dir --gene_info --gwas_sumstat_dir --ld_ref --egret_subdir --gwas_traits --conda_env --slurm_account --slurm_partition --slurm_log_dir --job_id_file)
CLI_VARS=(tissue output_dir gene_info gwas_sumstat_dir ld_ref egret_subdir gwas_traits conda_env slurm_account slurm_partition slurm_log_dir job_id_file)
CLI_REQUIRED=(tissue output_dir gene_info gwas_sumstat_dir ld_ref egret_subdir gwas_traits conda_env slurm_account slurm_partition slurm_log_dir job_id_file)
parse_cli "$@"; activate_workflow_env "$conda_env"
require_file "${output_dir}/results_sumstats/${tissue}/${egret_subdir}.txt" "EGRET analysis summary"
Rscript "${scripts_dir}/8.2_create_EGRET_pos.R" --tissue "$tissue" --output_dir "$output_dir" \
    --gene_info "$gene_info" --egret_subdir "$egret_subdir"
wgt="${output_dir}/pos_files/${tissue}/${egret_subdir}.pos"
wgtdir="${output_dir}/EGRET_models/${tissue}/${egret_subdir}"
out_dir="${output_dir}/TWAS_association_results/${tissue}/${egret_subdir}"
IFS=',' read -ra traits <<< "$gwas_traits"
job_ids=()
for gwas_file in "${traits[@]}"; do
    trait=${gwas_file%%.sumstats.gz}
    jid=$(submit_with_retry --parsable --job-name="${tissue}_twas_${trait}" --account="$slurm_account" \
        --partition="$slurm_partition" --output="${slurm_log_dir}/TWAS_output.%j.%N.out" \
        "${scripts_dir}/8.4_run_TWAS_association.sh" --gwas_sumstats "${gwas_sumstat_dir}/${gwas_file}" \
        --trait "$trait" --weights "$wgt" --weights_dir "$wgtdir" --ld_ref "$ld_ref" \
        --output_dir "$out_dir" --conda_env "$conda_env")
    [[ -n "$jid" ]] || cli_die "no job ID returned for EGRET TWAS trait $trait"
    job_ids+=("$jid")
done
write_job_ids "$job_id_file" "${job_ids[@]}"
