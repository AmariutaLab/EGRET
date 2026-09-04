#!/bin/bash
# --- 5.0_run_transPCO_stage.sh ---------------------------------------------
# Stage `transpco`. Multivariate trans mapping over co-expression modules: builds
# WGCNA modules and per-chromosome genotypes, then chains the per-module MatrixeQTL,
# PCO association and FDR jobs through SLURM dependencies.
#
#   in : expression_files/{tissue}_expression.txt.gz; covariates; genotype_files/
#   out: transPCO/{tissue}/fold_{i}/{transpco_bed_dir}/{gene}.txt
#   run: run_EGRET_workflow.sh, stage `transpco`
# --------------------------------------------------------------------------
set -euo pipefail
# Locate workflow_lib.sh; see its header for why BASH_SOURCE alone is not enough.
scripts_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd || true)
[[ -r "${scripts_dir}/workflow_lib.sh" ]] || scripts_dir=${EGRET_SCRIPTS_DIR:-${SLURM_SUBMIT_DIR:-}}
[[ -r "${scripts_dir}/workflow_lib.sh" ]] || { echo "ERROR: cannot locate workflow_lib.sh (tried BASH_SOURCE, EGRET_SCRIPTS_DIR, SLURM_SUBMIT_DIR); export EGRET_SCRIPTS_DIR to the EGRET scripts directory" >&2; exit 1; }
export EGRET_SCRIPTS_DIR="$scripts_dir"
source "${scripts_dir}/workflow_lib.sh"
USAGE="Usage: 5.0_run_transPCO_stage.sh --tissue NAME --folds N --fdr VALUE --covariate_columns LIST --gene_info FILE --plink_path PATH --output_dir DIR --conda_env NAME --slurm_account NAME --slurm_partition NAME --slurm_log_dir DIR --matrix_job_id_file FILE --pco_job_id_file FILE --job_id_file FILE [options]"
CLI_FLAGS=(--tissue --folds --fdr --covariate_columns --gene_info --plink_path --output_dir --genotype_prefix --conda_env --slurm_account --slurm_partition --slurm_log_dir --matrix_job_id_file --pco_job_id_file --job_id_file --bed_dir --cleanup_transpco_raw)
CLI_VARS=(tissue folds fdr covariate_columns gene_info plink_path output_dir genotype_prefix conda_env slurm_account slurm_partition slurm_log_dir matrix_job_id_file pco_job_id_file job_id_file bed_dir cleanup_transpco_raw)
CLI_REQUIRED=(tissue folds fdr covariate_columns gene_info plink_path output_dir conda_env slurm_account slurm_partition slurm_log_dir matrix_job_id_file pco_job_id_file job_id_file)
genotype_prefix=GTEX_v8_genotypes_pruned; cleanup_transpco_raw=FALSE
# Where 5.8 writes the per-gene BEDs. Passed straight through to 5.7, which substitutes the
# historical bed_files_FDR_{fdr} when this is empty. The caller must be able to set it: 6.2
# reads the BEDs back from the config's transpco_bed_dir.
bed_dir=""
parse_cli "$@"; activate_workflow_env "$conda_env"
Rscript "${scripts_dir}/5.1_run_WGCNA_clustering.R" \
    --expression "${output_dir}/expression_files/${tissue}_expression.txt.gz" --output_dir "$output_dir" \
    --individuals "${output_dir}/fold_0_info/${tissue}/train_individuals.txt" \
    --covariates "${output_dir}/covariate_files/${tissue}_covariates.txt.gz" --tissue "$tissue" \
    --folds "$folds" --gene_info "$gene_info" --covariate_columns_for_coexpression "$covariate_columns"
Rscript "${scripts_dir}/5.2_initialize_chr_genotypes.R" --tissue "$tissue" --output_dir "$output_dir" \
    --folds "$folds" --genotype_prefix "$genotype_prefix" --plink_path "$plink_path"
module_names=()
for ((fold=0; fold<=folds; fold++)); do
    module_dir="${output_dir}/transPCO/${tissue}/fold_${fold}/modules"
    [[ -d "$module_dir" ]] || cli_die "transPCO module directory is missing: $module_dir"
    for module_file in "$module_dir"/*; do
        [[ -s "$module_file" ]] || continue
        module_names+=("$(basename "$module_file")")
    done
done
((${#module_names[@]})) || cli_die "no non-empty transPCO module files found"
mapfile -t modules < <(printf '%s\n' "${module_names[@]}" | sort -u)
matrix_ids=()
for module in "${modules[@]}"; do
    jid=$(submit_with_retry --parsable --job-name="${tissue}_meqtl_${module%.*}" --account="$slurm_account" \
        --partition="$slurm_partition" --output="${slurm_log_dir}/transPCO_matrixeQTL.%j.%N.out" \
        "${scripts_dir}/5.3_run_MatrixeQTL_job.sh" --tissue "$tissue" --module "$module" \
        --module_dir modules --output_dir "$output_dir" --folds "$folds" --plink_path "$plink_path" \
        --genotype_prefix "$genotype_prefix" --conda_env "$conda_env")
    [[ -n "$jid" ]] || cli_die "no job ID returned for transPCO MatrixeQTL module $module"
    matrix_ids+=("$jid")
done
matrix_dep=$(join_by_colon "${matrix_ids[@]}")
write_job_ids "$matrix_job_id_file" "${matrix_ids[@]}"
pco_ids=()
for module in "${modules[@]}"; do
    jid=$(submit_with_retry --parsable --job-name="${tissue}_pco_${module%.*}" --account="$slurm_account" \
        --partition="$slurm_partition" --dependency="afterok:${matrix_dep}" \
        --output="${slurm_log_dir}/transPCO_PCO_assoc.%j.%N.out" \
        "${scripts_dir}/5.5_run_PCO_association_job.sh" --tissue "$tissue" --module "$module" \
        --module_dir modules --output_dir "$output_dir" --folds "$folds" --gene_info "$gene_info" \
        --genotype_prefix "$genotype_prefix" --conda_env "$conda_env")
    [[ -n "$jid" ]] || cli_die "no job ID returned for transPCO PCO module $module"
    pco_ids+=("$jid")
done
((${#pco_ids[@]})) || cli_die "no transPCO PCO jobs were submitted"
pco_dep=$(join_by_colon "${pco_ids[@]}")
write_job_ids "$pco_job_id_file" "${pco_ids[@]}"
final_jid=$(submit_with_retry --parsable --job-name="${tissue}_transpco_results" --account="$slurm_account" \
    --partition="$slurm_partition" --dependency="afterok:${pco_dep}" \
    --output="${slurm_log_dir}/transPCO_results.%j.%N.out" \
    "${scripts_dir}/5.7_run_results_analysis_job.sh" --tissue "$tissue" --fdr "$fdr" \
    --output_dir "$output_dir" --folds "$folds" --expected_module_count "${#modules[@]}" \
    --genotype_prefix "$genotype_prefix" --conda_env "$conda_env" \
    --bed_dir "$bed_dir" \
    --cleanup_transpco_raw "$cleanup_transpco_raw")
[[ -n "$final_jid" ]] || cli_die "no job ID returned for transPCO results analysis"
write_job_ids "$job_id_file" "$final_jid"
