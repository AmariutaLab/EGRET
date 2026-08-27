#!/bin/bash
#SBATCH --job-name="transPCO_PCO_assoc"
#SBATCH --output="batch_submissions/transPCO_PCO_assoc.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=12G
#SBATCH --account=csd832
#SBATCH --export=ALL
#SBATCH -t 12:00:00
# --- 5.5_run_PCO_association_job.sh ----------------------------------------
# Runs the PCO multivariate test for one module. Submitted once per module, after
# that module's 5.3 job.
#
#   in : transPCO/{tissue}/fold_{i}/association_results/
#   out: transPCO/{tissue}/fold_{i}/PCO_association_results/
#   run: submitted by 5.0_run_transPCO_stage.sh
# --------------------------------------------------------------------------
set -euo pipefail
# Locate workflow_lib.sh; see its header for why BASH_SOURCE alone is not enough.
scripts_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd || true)
[[ -r "${scripts_dir}/workflow_lib.sh" ]] || scripts_dir=${EGRET_SCRIPTS_DIR:-${SLURM_SUBMIT_DIR:-}}
[[ -r "${scripts_dir}/workflow_lib.sh" ]] || { echo "ERROR: cannot locate workflow_lib.sh (tried BASH_SOURCE, EGRET_SCRIPTS_DIR, SLURM_SUBMIT_DIR); export EGRET_SCRIPTS_DIR to the EGRET scripts directory" >&2; exit 1; }
export EGRET_SCRIPTS_DIR="$scripts_dir"
source "${scripts_dir}/workflow_lib.sh"
USAGE="Usage: 5.5_run_PCO_association_job.sh --tissue NAME --module FILE --module_dir DIR --output_dir DIR --folds N --gene_info FILE --conda_env NAME [--genotype_prefix PREFIX]"
CLI_FLAGS=(--tissue --module --module_dir --output_dir --folds --gene_info --conda_env --genotype_prefix)
CLI_VARS=(tissue module module_dir output_dir folds gene_info conda_env genotype_prefix)
CLI_REQUIRED=(tissue module module_dir output_dir folds gene_info conda_env)
genotype_prefix=GTEX_v8_genotypes_pruned
parse_cli "$@"; activate_workflow_env "$conda_env"
module_name=${module%.*}; missing=0; has_any=0
for ((fold=0; fold<=folds; fold++)); do
    module_file="${output_dir}/transPCO/${tissue}/fold_${fold}/${module_dir}/${module}"
    [[ -f "$module_file" ]] || continue
    assoc_dir="${output_dir}/transPCO/${tissue}/fold_${fold}/association_results"
    pco_dir="${output_dir}/transPCO/${tissue}/fold_${fold}/PCO_association_results"
    count=$(count_files "$assoc_dir" -maxdepth 1 -type f -name "matrix_${module_name}_chr_*.txt.gz")
    done_count=$(count_files "$pco_dir" -maxdepth 1 -type f -name "${module_name}_chr_*.txt.gz")
    if ((count == 0 && done_count == 0)); then missing=1; else has_any=1; fi
done
((has_any == 1 && missing == 0)) || cli_die "missing association results for module $module_name"
Rscript "${scripts_dir}/5.6_PCO_association_by_module.R" --tissue "$tissue" \
    --module "$module" --module_dir "$module_dir" --output_dir "$output_dir" --folds "$folds" \
    --gene_info "$gene_info" --results_dir PCO_association_results --association_dir association_results \
    --scripts_dir "$scripts_dir" --genotype_prefix "$genotype_prefix"
