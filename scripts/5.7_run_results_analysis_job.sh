#!/bin/bash
#SBATCH --job-name="transPCO_results"
#SBATCH --output="batch_submissions/transPCO_results.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=12G
#SBATCH --account=csd832
#SBATCH --export=ALL
#SBATCH -t 02:00:00
# --- 5.7_run_results_analysis_job.sh ---------------------------------------
# Runs the FDR step once every module's PCO job has finished.
#
#   in : transPCO/{tissue}/fold_{i}/PCO_association_results/
#   out: transPCO/{tissue}/fold_{i}/{transpco_bed_dir}/{gene}.txt
#   run: submitted by 5.0_run_transPCO_stage.sh
# --------------------------------------------------------------------------
set -euo pipefail
# Locate workflow_lib.sh; see its header for why BASH_SOURCE alone is not enough.
scripts_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd || true)
[[ -r "${scripts_dir}/workflow_lib.sh" ]] || scripts_dir=${EGRET_SCRIPTS_DIR:-${SLURM_SUBMIT_DIR:-}}
[[ -r "${scripts_dir}/workflow_lib.sh" ]] || { echo "ERROR: cannot locate workflow_lib.sh (tried BASH_SOURCE, EGRET_SCRIPTS_DIR, SLURM_SUBMIT_DIR); export EGRET_SCRIPTS_DIR to the EGRET scripts directory" >&2; exit 1; }
export EGRET_SCRIPTS_DIR="$scripts_dir"
source "${scripts_dir}/workflow_lib.sh"
USAGE="Usage: 5.7_run_results_analysis_job.sh --tissue NAME --fdr VALUE --output_dir DIR --folds N --expected_module_count N --genotype_prefix PREFIX --conda_env NAME [options]"
CLI_FLAGS=(--tissue --fdr --output_dir --folds --expected_module_count --genotype_prefix --conda_env --cleanup_transpco_raw)
CLI_VARS=(tissue fdr output_dir folds expected_module_count genotype_prefix conda_env cleanup_transpco_raw)
CLI_REQUIRED=(tissue fdr output_dir folds expected_module_count genotype_prefix conda_env)
cleanup_transpco_raw=FALSE
parse_cli "$@"; activate_workflow_env "$conda_env"
[[ "$expected_module_count" =~ ^[1-9][0-9]*$ ]] ||
    cli_die "--expected_module_count must be a positive integer"
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
((${#modules[@]} == expected_module_count)) ||
    cli_die "module universe changed: expected $expected_module_count, found ${#modules[@]}"
missing_results=()
for ((fold=0; fold<=folds; fold++)); do
    module_dir="${output_dir}/transPCO/${tissue}/fold_${fold}/modules"
    pco_dir="${output_dir}/transPCO/${tissue}/fold_${fold}/PCO_association_results"
    [[ -d "$module_dir" ]] || { missing_results+=("fold_${fold}:modules_missing"); continue; }
    for module_file in "$module_dir"/*; do
        [[ -s "$module_file" ]] || continue
        module=$(basename "$module_file")
        module_name=${module%%.*}
        # 5.6 skips any module with fewer than 2 genes -- PCO is a multivariate test and is
        # undefined on a single gene -- so requiring a result for one is a guaranteed stage
        # failure. WGCNA emits these routinely. Mirror that rule.
        gene_count=$(grep -c '[^[:space:]]' "$module_file" || true)
        ((gene_count >= 2)) || continue
        has_result=0
        for result_file in "${pco_dir}/${module_name}_chr_"*.txt.gz; do
            [[ -s "$result_file" ]] && { has_result=1; break; }
        done
        ((has_result == 1)) || missing_results+=("fold_${fold}:${module_name}")
    done
done
((${#missing_results[@]} == 0)) ||
    cli_die "missing PCO results for ${missing_results[*]}"
Rscript "${scripts_dir}/5.8_make_transPCO_bed_by_FDR.R" --tissue "$tissue" --FDR "$fdr" \
    --output_dir "$output_dir" --folds "$folds" --bed_dir "bed_files_FDR_${fdr}" \
    --PCO_association_dir PCO_association_results --module_dir modules --genotype_prefix "$genotype_prefix"
# Cleanup: the dependency-gated end of the transPCO branch, so the first point at which the
# per-module scratch is safe to drop. genotype_files and expression_files are shared across
# every module job, which is why no per-module job may remove them.
if [[ "$cleanup_transpco_raw" == TRUE ]]; then
    for ((fold=0; fold<=folds; fold++)); do
        fold_dir="${output_dir}/transPCO/${tissue}/fold_${fold}"
        for sub in PCO_association_results association_results genotype_files expression_files; do
            [[ ! -d "${fold_dir}/${sub}" ]] || rm -rf "${fold_dir}/${sub}"
        done
    done
fi
