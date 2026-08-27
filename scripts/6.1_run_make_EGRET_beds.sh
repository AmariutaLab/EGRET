#!/bin/bash
#SBATCH --job-name="make_EGRET_bed"
#SBATCH --output="batch_submissions/make_EGRET_bed.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=12G
#SBATCH --account=csd832
#SBATCH --export=ALL
#SBATCH -t 24:00:00
# --- 6.1_run_make_EGRET_beds.sh --------------------------------------------
# Builds the EGRET BED files for one fold. Submitted once per fold by 6.0.
#
#   in : MatrixeQTL/, GBAT/, transPCO/ bed trees; cross_mapped/
#   out: bed_files/{tissue}/{egret_subdir}/fold_{i}/{gene}; plink_results/{tissue}/{egret_subdir}/
#   run: submitted by 6.0_train_EGRET_models.sh
# --------------------------------------------------------------------------
set -euo pipefail
# Locate workflow_lib.sh; see its header for why BASH_SOURCE alone is not enough.
scripts_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd || true)
[[ -r "${scripts_dir}/workflow_lib.sh" ]] || scripts_dir=${EGRET_SCRIPTS_DIR:-${SLURM_SUBMIT_DIR:-}}
[[ -r "${scripts_dir}/workflow_lib.sh" ]] || { echo "ERROR: cannot locate workflow_lib.sh (tried BASH_SOURCE, EGRET_SCRIPTS_DIR, SLURM_SUBMIT_DIR); export EGRET_SCRIPTS_DIR to the EGRET scripts directory" >&2; exit 1; }
export EGRET_SCRIPTS_DIR="$scripts_dir"
source "${scripts_dir}/workflow_lib.sh"
USAGE="Usage: 6.1_run_make_EGRET_beds.sh --tissue NAME --output_dir DIR --egret_subdir NAME --plink_path PATH --expression FILE --fold N --gene_info FILE --conda_env NAME [options]"
CLI_FLAGS=(--tissue --output_dir --matrixeqtl_bed_dir --gbat_bed_dir --transpco_bed_dir --egret_subdir --exclude_crossmap --crossmap_dir --genotype_prefix --plink_path --expression --fold --gene_info --bed_models --conda_env --cleanup_bed_lists --cleanup_fold_expression)
CLI_VARS=(tissue output_dir matrixeqtl_bed_dir gbat_bed_dir transpco_bed_dir egret_subdir exclude_crossmap crossmap_dir genotype_prefix plink_path expression fold gene_info bed_models conda_env cleanup_bed_lists cleanup_fold_expression)
CLI_REQUIRED=(tissue output_dir egret_subdir plink_path expression fold gene_info conda_env)
matrixeqtl_bed_dir=results_FDR_0.1; gbat_bed_dir=results_FDR_0.1; transpco_bed_dir=bed_files_FDR_0.1
exclude_crossmap=TRUE; crossmap_dir=background_mismatches; genotype_prefix=GTEX_v8_genotypes_pruned
bed_models=Cis,MatrixeQTL,GBAT,transPCO; cleanup_bed_lists=FALSE; cleanup_fold_expression=FALSE
parse_cli "$@"; activate_workflow_env "$conda_env"
Rscript "${scripts_dir}/6.2_make_EGRET_beds_by_fold.R" --tissue "$tissue" \
    --output_dir "$output_dir" --matrixeqtl_bed_dir "$matrixeqtl_bed_dir" \
    --gbat_bed_dir "$gbat_bed_dir" --transpco_bed_dir "$transpco_bed_dir" \
    --egret_subdir "$egret_subdir" --exclude_crossmap "$exclude_crossmap" --crossmap_dir "$crossmap_dir" \
    --genotype_prefix "$genotype_prefix" --plink_path "$plink_path" --expression "$expression" \
    --bed_models "$bed_models" --fold "$fold" --gene_info "$gene_info" \
    --cleanup_bed_lists "$cleanup_bed_lists"

# Cleanup: train_expression.txt has three consumers in concurrent branches (3.1, 4.2, 5.4),
# so no single branch may delete it. This stage waits on all three, making it the first
# safe point.
if [[ "$cleanup_fold_expression" == TRUE ]]; then
    rm -f "${output_dir}/fold_${fold}_info/${tissue}/train_expression.txt"
fi
