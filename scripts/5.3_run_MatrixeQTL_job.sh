#!/bin/bash
#SBATCH --job-name="transPCO_matrixeQTL"
#SBATCH --output="batch_submissions/transPCO_matrixeQTL.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=12G
#SBATCH --account=csd832
#SBATCH --export=ALL
#SBATCH -t 12:00:00
# --- 5.3_run_MatrixeQTL_job.sh ---------------------------------------------
# Runs the per-module SNP-gene scan for one module. Submitted once per module.
#
#   in : transPCO/{tissue}/fold_{i}/modules/, genotype_files/
#   out: transPCO/{tissue}/fold_{i}/association_results/
#   run: submitted by 5.0_run_transPCO_stage.sh
# --------------------------------------------------------------------------
set -euo pipefail
# Locate workflow_lib.sh; see its header for why BASH_SOURCE alone is not enough.
scripts_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd || true)
[[ -r "${scripts_dir}/workflow_lib.sh" ]] || scripts_dir=${EGRET_SCRIPTS_DIR:-${SLURM_SUBMIT_DIR:-}}
[[ -r "${scripts_dir}/workflow_lib.sh" ]] || { echo "ERROR: cannot locate workflow_lib.sh (tried BASH_SOURCE, EGRET_SCRIPTS_DIR, SLURM_SUBMIT_DIR); export EGRET_SCRIPTS_DIR to the EGRET scripts directory" >&2; exit 1; }
export EGRET_SCRIPTS_DIR="$scripts_dir"
source "${scripts_dir}/workflow_lib.sh"
USAGE="Usage: 5.3_run_MatrixeQTL_job.sh --tissue NAME --module FILE --module_dir DIR --output_dir DIR --folds N --plink_path PATH --conda_env NAME [--genotype_prefix PREFIX]"
CLI_FLAGS=(--tissue --module --module_dir --output_dir --folds --plink_path --genotype_prefix --conda_env)
CLI_VARS=(tissue module module_dir output_dir folds plink_path genotype_prefix conda_env)
CLI_REQUIRED=(tissue module module_dir output_dir folds plink_path conda_env)
genotype_prefix=GTEX_v8_genotypes_pruned
parse_cli "$@"; activate_workflow_env "$conda_env"
Rscript "${scripts_dir}/5.4_MatrixeQTL_for_transPCO.R" --tissue "$tissue" --module "$module" \
    --module_dir "$module_dir" --output_dir "$output_dir" --folds "$folds" --plink_path "$plink_path" \
    --genotype_prefix "$genotype_prefix"
