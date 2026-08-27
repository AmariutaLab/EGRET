#!/bin/bash
# --- 1.0_setup_genotype_and_expression.sh ----------------------------------
# Stage `setup`. Runs the three preprocessing steps in order: LD-prune genotypes,
# filter and residualize expression, then assign cross-validation folds.
#
#   in : genotype_prefix, expression, covariates, individuals, gene_info
#   out: genotype_files/, expression_files/, covariate_files/, fold_{i}_info/
#   run: run_EGRET_workflow.sh, stage `setup`
# --------------------------------------------------------------------------
set -euo pipefail
# Locate workflow_lib.sh; see its header for why BASH_SOURCE alone is not enough.
scripts_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd || true)
[[ -r "${scripts_dir}/workflow_lib.sh" ]] || scripts_dir=${EGRET_SCRIPTS_DIR:-${SLURM_SUBMIT_DIR:-}}
[[ -r "${scripts_dir}/workflow_lib.sh" ]] || { echo "ERROR: cannot locate workflow_lib.sh (tried BASH_SOURCE, EGRET_SCRIPTS_DIR, SLURM_SUBMIT_DIR); export EGRET_SCRIPTS_DIR to the EGRET scripts directory" >&2; exit 1; }
export EGRET_SCRIPTS_DIR="$scripts_dir"
source "${scripts_dir}/workflow_lib.sh"
USAGE="Usage: 1.0_setup_genotype_and_expression.sh --genotype_prefix PATH --expression FILE --covariates FILE --tissue NAME --individuals FILE --plink_path PATH --gene_info FILE --output_dir DIR --conda_env NAME [options]"
CLI_FLAGS=(--genotype_prefix --expression --covariates --tissue --individuals --ld_prune_r2 --plink_path --genotype_output_prefix --folds --gene_info --output_dir --conda_env)
CLI_VARS=(genotype_prefix expression covariates tissue individuals ld_prune_r2 plink_path genotype_output_prefix folds gene_info output_dir conda_env)
CLI_REQUIRED=(genotype_prefix expression covariates tissue individuals plink_path gene_info output_dir conda_env)
ld_prune_r2=0.9; genotype_output_prefix=GTEX_v8_genotypes_pruned; folds=5
parse_cli "$@"; activate_workflow_env "$conda_env"
mkdir -p "${output_dir}/genotype_files"
Rscript "${scripts_dir}/1.1_setup_genotypes.R" --plink_path "$plink_path" --genotype_prefix "$genotype_prefix" \
    --LD_r2 "$ld_prune_r2" --LD_chunk 1 --LD_window 100 --out "${output_dir}/genotype_files/${genotype_output_prefix}"
Rscript "${scripts_dir}/1.2_setup_expression.R" --expression "$expression" --individuals "$individuals" \
    --covariates "$covariates" --tissue "$tissue" --gene_info "$gene_info" --output_dir "$output_dir"
Rscript "${scripts_dir}/1.3_setup_folds.R" --expression "${output_dir}/expression_files/${tissue}_expression_regressed.txt.gz" \
    --folds "$folds" --tissue "$tissue" --individuals "$individuals" --plink_path "$plink_path" \
    --genotype_prefix "${output_dir}/genotype_files/${genotype_output_prefix}" --output_dir "$output_dir"
