#!/bin/bash
#SBATCH --job-name="MatrixeQTL"
#SBATCH --output="batch_submissions/MatrixeQTL.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=64G
#SBATCH --account=csd832
#SBATCH --export=ALL
#SBATCH -t 24:00:00
# --- 3.0_run_MatrixeQTL_stage.sh -------------------------------------------
# Stage `matrixeqtl`. Runs the univariate SNP-gene scan across folds, then converts
# the significant hits into per-gene BED files for model training.
#
#   in : fold_{i}_info/{tissue}/train_{expression,genotypes}.txt
#   out: MatrixeQTL/{tissue}/fold_{i}/{matrixeqtl_bed_dir}/{gene}.txt
#   run: run_EGRET_workflow.sh, stage `matrixeqtl`
# --------------------------------------------------------------------------
set -euo pipefail
# Locate workflow_lib.sh; see its header for why BASH_SOURCE alone is not enough.
scripts_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd || true)
[[ -r "${scripts_dir}/workflow_lib.sh" ]] || scripts_dir=${EGRET_SCRIPTS_DIR:-${SLURM_SUBMIT_DIR:-}}
[[ -r "${scripts_dir}/workflow_lib.sh" ]] || { echo "ERROR: cannot locate workflow_lib.sh (tried BASH_SOURCE, EGRET_SCRIPTS_DIR, SLURM_SUBMIT_DIR); export EGRET_SCRIPTS_DIR to the EGRET scripts directory" >&2; exit 1; }
export EGRET_SCRIPTS_DIR="$scripts_dir"
source "${scripts_dir}/workflow_lib.sh"
USAGE="Usage: 3.0_run_MatrixeQTL_stage.sh --tissue NAME --fdr VALUE --folds N --gene_info FILE --genotype_prefix PREFIX --output_dir DIR --conda_env NAME [--cleanup_matrixeqtl_raw BOOL] [--cleanup_fold_genotypes BOOL] [--cis_window_size BP]"
CLI_FLAGS=(--tissue --fdr --folds --gene_info --genotype_prefix --output_dir --conda_env --cleanup_matrixeqtl_raw --cleanup_fold_genotypes --cis_window_size)
CLI_VARS=(tissue fdr folds gene_info genotype_prefix output_dir conda_env cleanup_matrixeqtl_raw cleanup_fold_genotypes cis_window_size)
CLI_REQUIRED=(tissue fdr folds gene_info genotype_prefix output_dir conda_env)
cleanup_matrixeqtl_raw=FALSE; cleanup_fold_genotypes=FALSE
cis_window_size=1000000
parse_cli "$@"; activate_workflow_env "$conda_env"
Rscript "${scripts_dir}/3.1_run_MatrixeQTL.R" --tissue "$tissue" --folds "$folds" --output_dir "$output_dir"
Rscript "${scripts_dir}/3.2_make_MatrixeQTL_bed_by_FDR.R" --tissue "$tissue" --folds "$folds" \
    --FDR "$fdr" --gene_info "$gene_info" --genotype_prefix "$genotype_prefix" --output_dir "$output_dir" \
    --cis_window_size "$cis_window_size"

# Cleanup: 3.2 above is the raw association tables' last consumer, so this stage owns
# their deletion. Note 3.1 writes them at MatrixeQTL/{tissue}/, not one level deeper.
if [[ "$cleanup_matrixeqtl_raw" == TRUE ]]; then
    find "${output_dir}/MatrixeQTL/${tissue}" -maxdepth 1 -type f \
        -name "association_results_fold_*_*.txt" -delete
fi

# Cleanup: fold_{i}_info/{tissue}/train_genotypes.txt is written by 1.3_setup_folds.R and
# read only by 3.1_run_MatrixeQTL.R above, so this stage is its last consumer.
if [[ "$cleanup_fold_genotypes" == TRUE ]]; then
    for ((fold=0; fold<=folds; fold++)); do
        rm -f "${output_dir}/fold_${fold}_info/${tissue}/train_genotypes.txt"
    done
fi
