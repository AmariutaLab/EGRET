#!/bin/bash
# --- 6.0_train_EGRET_models.sh ---------------------------------------------
# Stage `models`. Builds the per-gene EGRET BED files for every fold, then chunks
# the gene list and submits one 6.3 training job per chunk.
#
#   in : the three trans-eQTL bed trees; plink_results/{tissue}/cis/
#   out: EGRET_models/{tissue}/{egret_subdir}/{gene}.wgt.RDat
#   run: run_EGRET_workflow.sh, stage `models`
# --------------------------------------------------------------------------
set -euo pipefail
# Locate workflow_lib.sh; see its header for why BASH_SOURCE alone is not enough.
scripts_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd || true)
[[ -r "${scripts_dir}/workflow_lib.sh" ]] || scripts_dir=${EGRET_SCRIPTS_DIR:-${SLURM_SUBMIT_DIR:-}}
[[ -r "${scripts_dir}/workflow_lib.sh" ]] || { echo "ERROR: cannot locate workflow_lib.sh (tried BASH_SOURCE, EGRET_SCRIPTS_DIR, SLURM_SUBMIT_DIR); export EGRET_SCRIPTS_DIR to the EGRET scripts directory" >&2; exit 1; }
export EGRET_SCRIPTS_DIR="$scripts_dir"
source "${scripts_dir}/workflow_lib.sh"
USAGE="Usage: 6.0_train_EGRET_models.sh --tissue NAME --output_dir DIR --expression FILE --gene_info FILE --plink_path PATH --gemma_path PATH --bed_job_id_file FILE --job_id_file FILE --conda_env NAME --slurm_account NAME --slurm_partition NAME --slurm_log_dir DIR [options]"
CLI_FLAGS=(--tissue --folds --plink_path --output_dir --matrixeqtl_bed_dir --gbat_bed_dir --transpco_bed_dir --exclude_crossmap --crossmap_dir --genotype_prefix --egret_subdir --fusion_models --bed_models --gemma_path --chunk_size --gene_info --expression --conda_env --slurm_account --slurm_partition --slurm_log_dir --bed_job_id_file --job_id_file --upstream_job_ids --cleanup_egret_plink --cleanup_egret_working --cleanup_egret_z_matrices --cleanup_bed_lists --cleanup_fold_expression)
CLI_VARS=(tissue folds plink_path output_dir matrixeqtl_bed_dir gbat_bed_dir transpco_bed_dir exclude_crossmap crossmap_dir genotype_prefix egret_subdir fusion_models bed_models gemma_path chunk_size gene_info expression conda_env slurm_account slurm_partition slurm_log_dir bed_job_id_file job_id_file upstream_job_ids cleanup_egret_plink cleanup_egret_working cleanup_egret_z_matrices cleanup_bed_lists cleanup_fold_expression)
CLI_REQUIRED=(tissue output_dir plink_path gemma_path gene_info expression conda_env slurm_account slurm_partition slurm_log_dir bed_job_id_file job_id_file)
folds=5; matrixeqtl_bed_dir=results_FDR_0.1; gbat_bed_dir=results_FDR_0.1
transpco_bed_dir=bed_files_FDR_0.1; exclude_crossmap=TRUE; crossmap_dir=background_mismatches
genotype_prefix=GTEX_v8_genotypes_pruned; egret_subdir=EGRET; fusion_models=xtune,lasso,enet,blup
bed_models=Cis,MatrixeQTL,GBAT,transPCO; chunk_size=500; upstream_job_ids=""
cleanup_egret_plink=FALSE; cleanup_egret_working=FALSE; cleanup_egret_z_matrices=FALSE
cleanup_bed_lists=FALSE; cleanup_fold_expression=FALSE
parse_cli "$@"
require_file "$expression" "regressed expression"
bed_ids=()
for ((fold=0; fold<=folds; fold++)); do
    args=(--parsable --job-name="${tissue}_bed${fold}" --account="$slurm_account" \
        --partition="$slurm_partition" --output="${slurm_log_dir}/make_EGRET_bed.%j.%N.out")
    [[ -n "$upstream_job_ids" ]] && args+=(--dependency="afterok:${upstream_job_ids}")
    jid=$(submit_with_retry "${args[@]}" "${scripts_dir}/6.1_run_make_EGRET_beds.sh" --tissue "$tissue" \
        --output_dir "$output_dir" --matrixeqtl_bed_dir "$matrixeqtl_bed_dir" \
        --gbat_bed_dir "$gbat_bed_dir" --transpco_bed_dir "$transpco_bed_dir" \
        --egret_subdir "$egret_subdir" --exclude_crossmap "$exclude_crossmap" \
        --crossmap_dir "$crossmap_dir" --genotype_prefix "$genotype_prefix" \
        --plink_path "$plink_path" --expression "$expression" --fold "$fold" \
        --gene_info "$gene_info" --bed_models "$bed_models" --conda_env "$conda_env" \
        --cleanup_bed_lists "$cleanup_bed_lists" --cleanup_fold_expression "$cleanup_fold_expression")
    [[ -n "$jid" ]] || cli_die "no job ID returned for bed fold $fold"
    bed_ids+=("$jid")
done
bed_dep=$(join_by_colon "${bed_ids[@]}")
num_genes=$(zcat "$expression" | tail -n +2 | wc -l)
write_job_ids "$bed_job_id_file" "${bed_ids[@]}"
((num_genes > 0)) || cli_die "expression file contains no genes: $expression"
num_chunks=$(((num_genes + chunk_size - 1) / chunk_size)); model_ids=()
for ((i=0; i<num_chunks; i++)); do
    start_index=$((i * chunk_size))
    jid=$(submit_with_retry --parsable --job-name="${tissue}_mdl${i}" --account="$slurm_account" \
        --partition="$slurm_partition" --dependency="afterok:${bed_dep}" \
        --output="${slurm_log_dir}/gtex_xtune.%j.%N.out" "${scripts_dir}/6.3_train_EGRET_chunk.sh" \
        --tissue "$tissue" --chunk_size "$chunk_size" --start_index "$start_index" \
        --output_dir "$output_dir" --egret_subdir "$egret_subdir" --folds "$folds" \
        --plink_path "$plink_path" --gemma_path "$gemma_path" --expression "$expression" \
        --fusion_models "$fusion_models" --gene_info "$gene_info" --conda_env "$conda_env" \
        --cleanup_egret_plink "$cleanup_egret_plink" --cleanup_egret_working "$cleanup_egret_working" \
        --cleanup_egret_z_matrices "$cleanup_egret_z_matrices")
    [[ -n "$jid" ]] || cli_die "no job ID returned for model chunk $i"
    model_ids+=("$jid")
done
write_job_ids "$job_id_file" "${model_ids[@]}"
