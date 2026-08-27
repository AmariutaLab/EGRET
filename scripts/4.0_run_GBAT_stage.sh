#!/bin/bash
#SBATCH --job-name="GBAT"
#SBATCH --output="batch_submissions/GBAT.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32G
#SBATCH --account=csd832
#SBATCH --export=ALL
#SBATCH -t 24:00:00
# --- 4.0_run_GBAT_stage.sh -------------------------------------------------
# Stage `gbat`. Polygenic trans-association: predicts cis expression from the cis
# FUSION weights, tests those predictions against every other gene, FDR-filters.
#
#   in : FUSION/{tissue}/cis/*.wgt.RDat; genotype panel
#   out: GBAT/{tissue}/fold_{i}/{gbat_bed_dir}/{gene}.txt
#   run: run_EGRET_workflow.sh, stage `gbat` (depends on stage `cis`)
# --------------------------------------------------------------------------
set -euo pipefail
# Locate workflow_lib.sh; see its header for why BASH_SOURCE alone is not enough.
scripts_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd || true)
[[ -r "${scripts_dir}/workflow_lib.sh" ]] || scripts_dir=${EGRET_SCRIPTS_DIR:-${SLURM_SUBMIT_DIR:-}}
[[ -r "${scripts_dir}/workflow_lib.sh" ]] || { echo "ERROR: cannot locate workflow_lib.sh (tried BASH_SOURCE, EGRET_SCRIPTS_DIR, SLURM_SUBMIT_DIR); export EGRET_SCRIPTS_DIR to the EGRET scripts directory" >&2; exit 1; }
export EGRET_SCRIPTS_DIR="$scripts_dir"
source "${scripts_dir}/workflow_lib.sh"
USAGE="Usage: 4.0_run_GBAT_stage.sh --tissue NAME --fdr VALUE --folds N --gene_info FILE --cis_model_dir DIR --plink_path PATH --genotype_prefix PREFIX --output_dir DIR --conda_env NAME [options]"
CLI_FLAGS=(--tissue --fdr --folds --gene_info --cis_model_dir --plink_path --genotype_prefix --output_dir --conda_env --cleanup_cis_plink)
CLI_VARS=(tissue fdr folds gene_info cis_model_dir plink_path genotype_prefix output_dir conda_env cleanup_cis_plink)
CLI_REQUIRED=(tissue fdr folds gene_info cis_model_dir plink_path genotype_prefix output_dir conda_env)
cleanup_cis_plink=FALSE
parse_cli "$@"; activate_workflow_env "$conda_env"
Rscript "${scripts_dir}/4.1_run_GBAT.R" --tissue "$tissue" --gene_info "$gene_info" \
    --output_dir "$output_dir" --cis_model_dir "$cis_model_dir" --genotype_prefix "$genotype_prefix" --plink_path "$plink_path"
Rscript "${scripts_dir}/4.2_run_GBAT_association.R" --tissue "$tissue" --folds "$folds" --output_dir "$output_dir"
Rscript "${scripts_dir}/4.3_make_GBAT_bed_by_FDR.R" --tissue "$tissue" --gene_info "$gene_info" \
    --cis_model_dir "$cis_model_dir" --FDR "$fdr" --folds "$folds" --output_dir "$output_dir"
if [[ "$cleanup_cis_plink" == TRUE ]]; then
    cis_plink_dir="${output_dir}/plink_results/${tissue}/cis"
    [[ ! -d "$cis_plink_dir" ]] || rm -rf "$cis_plink_dir"
fi
