#!/bin/bash
#SBATCH --job-name="EGRET_workflow"
#SBATCH --output="EGRET_workflow.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=30G
#SBATCH --account=csd832
#SBATCH --export=ALL
#SBATCH -t 08:00:00
# --- run_EGRET_workflow.sh -------------------------------------------------
# EGRET entry point. Validates the configuration, then submits every enabled stage as
# a SLURM DAG wired together by job dependencies.
# Stages it does NOT run still have their outputs asserted, so a resumed run cannot
# silently train on missing inputs.
#
#   in : a config file (--config, default config.sh)
#   out: the whole output_dir tree; job IDs under output_dir/.job_ids/{run_id}/
#   run: sbatch run_EGRET_workflow.sh --config config.sh
# --------------------------------------------------------------------------

set -euo pipefail
# Locate workflow_lib.sh; see its header for why BASH_SOURCE alone is not enough.
scripts_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd || true)
[[ -r "${scripts_dir}/workflow_lib.sh" ]] || scripts_dir=${EGRET_SCRIPTS_DIR:-${SLURM_SUBMIT_DIR:-}}
[[ -r "${scripts_dir}/workflow_lib.sh" ]] || { echo "ERROR: cannot locate workflow_lib.sh (tried BASH_SOURCE, EGRET_SCRIPTS_DIR, SLURM_SUBMIT_DIR); export EGRET_SCRIPTS_DIR to the EGRET scripts directory" >&2; exit 1; }
export EGRET_SCRIPTS_DIR="$scripts_dir"
source "${scripts_dir}/workflow_lib.sh"

USAGE="Usage: run_EGRET_workflow.sh [--config FILE]\nConfiguration defaults to \${EGRET_CONFIG:-scripts/config.sh}.\nNote: this script's own SLURM log lands in the submit directory, not in slurm_log_dir --\nSLURM reads the #SBATCH directives above before config.sh is sourced, so the directive\ncannot reference the configured path. Every job it submits does use slurm_log_dir."
config_file=${EGRET_CONFIG:-${scripts_dir}/config.sh}
CLI_FLAGS=(--config); CLI_VARS=(config_file); CLI_REQUIRED=()
parse_cli "$@"

[[ -r "$config_file" ]] || cli_die "configuration file is not readable: $config_file"
source "$config_file"
grn_pval_threshold=${grn_pval_threshold:-0.01}
required_config=(genotype_prefix expression individuals covariates gene_info plink_path plink1_path
    gcta_path gemma_path tissue output_dir genotype_output_prefix egret_subdir folds ld_prune_r2 fdr chunk_size
    bed_models fusion_models matrixeqtl_bed_dir gbat_bed_dir transpco_bed_dir stages conda_env
    exclude_crossmap crossmap_dir covariate_columns_for_coexpression cleanup_cis_plink
    cleanup_cis_working cleanup_matrixeqtl_raw cleanup_transpco_raw cleanup_egret_plink
    cleanup_egret_working
    cleanup_egret_z_matrices cleanup_bed_lists cleanup_fold_genotypes cleanup_fold_expression
    slurm_account slurm_partition slurm_log_dir gwas_sumstat_dir ld_ref gwas_sumstats)
for config_name in "${required_config[@]}"; do
    [[ -n "${!config_name:-}" ]] || cli_die "missing configuration parameter: $config_name"
done
for stage_name in $stages; do
    case "$stage_name" in
        setup|cis|matrixeqtl|gbat|transpco|models|analysis|twas|twas_cis|grn) ;;
        *) cli_die "unknown stage in configuration: $stage_name" ;;
    esac
done
[[ -n "$stages" ]] || cli_die "configuration parameter stages must not be empty"
boolean_config=(exclude_crossmap cleanup_cis_plink cleanup_cis_working cleanup_matrixeqtl_raw
    cleanup_transpco_raw cleanup_egret_plink cleanup_egret_working cleanup_egret_z_matrices
    cleanup_bed_lists cleanup_fold_genotypes cleanup_fold_expression)
for config_name in "${boolean_config[@]}"; do
    [[ "${!config_name}" == TRUE || "${!config_name}" == FALSE ]] ||
        cli_die "configuration parameter $config_name must be TRUE or FALSE"
done
if [[ $(declare -p gwas_sumstats 2>/dev/null) == "declare -a"* ]]; then
    gwas_traits=$(IFS=,; echo "${gwas_sumstats[*]}")
else
    gwas_traits=$gwas_sumstats
fi

stage_enabled() { [[ " $stages " == *" $1 "* ]]; }
any_stage_enabled() {
    local stage
    for stage in "$@"; do
        stage_enabled "$stage" && return 0
    done
    return 1
}
model_uses_source() { [[ ",${bed_models// /}," == *",$1,"* ]]; }

if stage_enabled grn; then
    [[ "$grn_pval_threshold" =~ ^(0([.][0-9]+)?|1([.]0+)?)$ ]] &&
        awk -v p="$grn_pval_threshold" 'BEGIN { exit !(p > 0 && p <= 1) }' ||
        cli_die "configuration parameter grn_pval_threshold must be greater than 0 and at most 1"
fi

# required_config only proves these are non-empty; check they exist before submitting
# anything, or a typo surfaces hours later inside a chunk job. Every problem is collected
# and reported together rather than one per submission.
config_problems=()
check_path() {   # check_path VARNAME KIND
    local name=$1 kind=$2 path=${!1}
    case "$kind" in
        file) [[ -f "$path" ]] || config_problems+=("$name: no such file: $path") ;;
        exec) if [[ -x "$path" ]]; then :
              elif command -v "$path" >/dev/null 2>&1; then :
              else config_problems+=("$name: not an executable and not on PATH: $path"); fi ;;
        dir)  [[ -d "$path" ]] || config_problems+=("$name: no such directory: $path") ;;
        # A PLINK binary prefix: .bed/.bim/.fam must all be present.
        bfile) local missing=()
               for ext in bed bim fam; do [[ -f "${path}.${ext}" ]] || missing+=(".$ext"); done
               ((${#missing[@]} == 0)) || config_problems+=("$name: missing ${missing[*]} for prefix $path") ;;
    esac
}
check_match() {  # check_match LABEL DIR GLOB
    local label=$1 dir=$2 pattern=$3
    [[ -d "$dir" ]] && compgen -G "${dir}/${pattern}" >/dev/null ||
        config_problems+=("$label: no files matching ${dir}/${pattern}")
}
if stage_enabled setup; then
    check_path genotype_prefix bfile
    check_path expression file
fi
if any_stage_enabled setup cis; then
    check_path individuals file
    check_path covariates file
fi
if any_stage_enabled setup cis matrixeqtl gbat transpco models twas twas_cis grn; then
    check_path gene_info file
fi
if any_stage_enabled setup cis matrixeqtl gbat transpco models; then
    check_path plink_path exec
fi
if stage_enabled cis; then
    check_path plink1_path exec
    check_path gcta_path exec
fi
if any_stage_enabled cis models; then
    check_path gemma_path exec
fi
if any_stage_enabled twas twas_cis grn; then
    check_path gwas_sumstat_dir dir
    IFS=',' read -ra configured_gwas <<< "$gwas_traits"
    for gwas_file in "${configured_gwas[@]}"; do
        [[ -n "$gwas_file" ]] || { config_problems+=("gwas_sumstats: contains an empty trait name"); continue; }
        [[ -f "${gwas_sumstat_dir}/${gwas_file}" ]] ||
            config_problems+=("gwas_sumstats: no such file: ${gwas_sumstat_dir}/${gwas_file}")
    done
fi
if stage_enabled twas || stage_enabled twas_cis; then
    check_path ld_ref bfile
fi
if stage_enabled grn; then
    if [[ -n "${eqtl_file:-}" ]]; then
        check_path eqtl_file file
    else
        config_problems+=("eqtl_file: missing configuration parameter (required by grn)")
    fi
fi
if stage_enabled models && [[ "$exclude_crossmap" == TRUE ]]; then
    check_match "cross-mappability BEDs" "${output_dir}/cross_mapped/${crossmap_dir}" '*.bed'
fi
if ((${#config_problems[@]})); then
    printf 'ERROR: configuration refers to paths that do not exist (%s):\n' "$config_file" >&2
    printf '  - %s\n' "${config_problems[@]}" >&2
    exit 2
fi

run_id=${EGRET_RUN_ID:-${SLURM_JOB_ID:-manual_$(date +%Y%m%dT%H%M%S)_$$}}
job_id_dir="${output_dir}/.job_ids/${run_id}"
run_log_dir="${slurm_log_dir}/${run_id}"
mkdir -p "$job_id_dir" "$run_log_dir"
activate_workflow_env "$conda_env"

submit_stage() {
    local name=$1 dependency=$2 log_name=$3 script=$4
    shift 4
    local args=(--parsable --job-name="$name" --account="$slurm_account"
        --partition="$slurm_partition" --output="${run_log_dir}/${log_name}.%j.%N.out")
    [[ -n "$dependency" ]] && args+=(--dependency="afterok:${dependency}")
    submit_with_retry "${args[@]}" "$script" "$@"
}
record_single_job() {
    local file=$1 jid=$2
    [[ -n "$jid" ]] || cli_die "submit_with_retry returned an empty job ID for $file"
    write_job_ids "$file" "$jid"
}

echo "EGRET run ID: $run_id"
echo "Current-run job IDs: $job_id_dir"

if stage_enabled setup; then
    "${scripts_dir}/1.0_setup_genotype_and_expression.sh" \
        --genotype_prefix "$genotype_prefix" --expression "$expression" \
        --covariates "$covariates" --tissue "$tissue" --individuals "$individuals" \
        --ld_prune_r2 "$ld_prune_r2" --plink_path "$plink_path" \
        --genotype_output_prefix "$genotype_output_prefix" --folds "$folds" \
        --gene_info "$gene_info" --output_dir "$output_dir" --conda_env "$conda_env"
elif any_stage_enabled cis matrixeqtl gbat transpco models; then
    require_file "${output_dir}/genotype_files/${genotype_output_prefix}.bim" "setup genotype BIM"
    require_file "${output_dir}/expression_files/${tissue}_expression_regressed.txt.gz" "regressed expression"
    for ((fold=0; fold<=folds; fold++)); do
        require_file "${output_dir}/fold_${fold}_info/${tissue}/train_individuals.txt" "fold_${fold} assignments"
    done
fi

cis_ids=""
if stage_enabled cis; then
    cis_file="${job_id_dir}/cis.jids"
    "${scripts_dir}/2.0_run_cis_models.sh" --tissue "$tissue" --output_dir "$output_dir" \
        --plink_path "$plink_path" --gene_info "$gene_info" \
        --genotype_prefix "$genotype_output_prefix" --individuals "$individuals" \
        --covariates "$covariates" --plink1_path "$plink1_path" --gcta_path "$gcta_path" \
        --gemma_path "$gemma_path" --chunk_size "$chunk_size" --conda_env "$conda_env" \
        --slurm_account "$slurm_account" --slurm_partition "$slurm_partition" \
        --slurm_log_dir "$run_log_dir" --job_id_file "$cis_file" \
        --cleanup_cis_working "$cleanup_cis_working"
    cis_ids=$(read_job_ids "$cis_file")
else
    if any_stage_enabled gbat analysis twas_cis; then
        require_dir_with_match "${output_dir}/FUSION/${tissue}/cis" '*.wgt.RDat' "cis model weights"
    fi
    if stage_enabled models && model_uses_source Cis; then
        require_dir_with_match "${output_dir}/bed_files/${tissue}/cis" '*.bed' "cis model BED windows"
    fi
fi

matrix_ids=""
if stage_enabled matrixeqtl; then
    jid=$(submit_stage "${tissue}_matrixeqtl" "" MatrixeQTL "${scripts_dir}/3.0_run_MatrixeQTL_stage.sh" \
        --tissue "$tissue" --fdr "$fdr" --folds "$folds" --gene_info "$gene_info" \
        --genotype_prefix "$genotype_output_prefix" --output_dir "$output_dir" --conda_env "$conda_env" \
        --cleanup_matrixeqtl_raw "$cleanup_matrixeqtl_raw" \
        --cleanup_fold_genotypes "$cleanup_fold_genotypes")
    record_single_job "${job_id_dir}/matrixeqtl.jids" "$jid"; matrix_ids=$jid
elif stage_enabled models && model_uses_source MatrixeQTL; then
    for ((fold=0; fold<=folds; fold++)); do
        require_dir_with_match "${output_dir}/MatrixeQTL/${tissue}/fold_${fold}/${matrixeqtl_bed_dir}" '*.txt' "MatrixeQTL fold_${fold} bed source"
    done
fi

gbat_ids=""
if stage_enabled gbat; then
    jid=$(submit_stage "${tissue}_gbat" "$cis_ids" GBAT "${scripts_dir}/4.0_run_GBAT_stage.sh" \
        --tissue "$tissue" --fdr "$fdr" --folds "$folds" --gene_info "$gene_info" \
        --cis_model_dir "${output_dir}/FUSION/${tissue}/cis" --plink_path "$plink_path" \
        --genotype_prefix "$genotype_prefix" --output_dir "$output_dir" \
        --cleanup_cis_plink "$cleanup_cis_plink" --conda_env "$conda_env")
    record_single_job "${job_id_dir}/gbat.jids" "$jid"; gbat_ids=$jid
elif stage_enabled models && model_uses_source GBAT; then
    for ((fold=0; fold<=folds; fold++)); do
        require_dir_with_match "${output_dir}/GBAT/${tissue}/fold_${fold}/${gbat_bed_dir}" '*.txt' "GBAT fold_${fold} bed source"
    done
fi

transpco_ids=""
if stage_enabled transpco; then
    transpco_file="${job_id_dir}/transpco.jids"
    "${scripts_dir}/5.0_run_transPCO_stage.sh" --tissue "$tissue" --folds "$folds" --fdr "$fdr" \
        --covariate_columns "$covariate_columns_for_coexpression" --gene_info "$gene_info" \
        --plink_path "$plink_path" --output_dir "$output_dir" \
        --genotype_prefix "$genotype_output_prefix" --conda_env "$conda_env" \
        --slurm_account "$slurm_account" --slurm_partition "$slurm_partition" \
        --slurm_log_dir "$run_log_dir" \
        --matrix_job_id_file "${job_id_dir}/transpco_matrix.jids" \
        --pco_job_id_file "${job_id_dir}/transpco_pco.jids" --job_id_file "$transpco_file" \
        --bed_dir "$transpco_bed_dir" \
        --cleanup_transpco_raw "$cleanup_transpco_raw"
    transpco_ids=$(read_job_ids "$transpco_file")
elif stage_enabled models && model_uses_source transPCO; then
    for ((fold=0; fold<=folds; fold++)); do
        require_dir_with_match "${output_dir}/transPCO/${tissue}/fold_${fold}/${transpco_bed_dir}" '*.txt' "transPCO fold_${fold} bed source"
    done
fi

upstream_ids=()
# 2.0 creates the cis BED windows synchronously before submitting its chunks, so the bed
# stage need not wait on them. Analysis does, and takes cis_ids directly below.
[[ -n "$matrix_ids" ]] && upstream_ids+=("$matrix_ids")
[[ -n "$gbat_ids" ]] && upstream_ids+=("$gbat_ids")
[[ -n "$transpco_ids" ]] && upstream_ids+=("$transpco_ids")
model_ids=""
if stage_enabled models; then
    model_file="${job_id_dir}/models.jids"
    "${scripts_dir}/6.0_train_EGRET_models.sh" --tissue "$tissue" --folds "$folds" \
        --plink_path "$plink_path" --output_dir "$output_dir" \
        --matrixeqtl_bed_dir "$matrixeqtl_bed_dir" --gbat_bed_dir "$gbat_bed_dir" \
        --transpco_bed_dir "$transpco_bed_dir" --exclude_crossmap "$exclude_crossmap" \
        --crossmap_dir "$crossmap_dir" --genotype_prefix "$genotype_output_prefix" \
        --egret_subdir "$egret_subdir" --fusion_models "$fusion_models" \
        --bed_models "$bed_models" --gemma_path "$gemma_path" --chunk_size "$chunk_size" \
        --gene_info "$gene_info" --expression "${output_dir}/expression_files/${tissue}_expression_regressed.txt.gz" \
        --conda_env "$conda_env" --slurm_account "$slurm_account" \
        --slurm_partition "$slurm_partition" --slurm_log_dir "$run_log_dir" \
        --bed_job_id_file "${job_id_dir}/beds.jids" --job_id_file "$model_file" \
        --upstream_job_ids "$(join_by_colon "${upstream_ids[@]}")" \
        --cleanup_egret_plink "$cleanup_egret_plink" \
        --cleanup_egret_working "$cleanup_egret_working" \
        --cleanup_egret_z_matrices "$cleanup_egret_z_matrices" \
        --cleanup_bed_lists "$cleanup_bed_lists" \
        --cleanup_fold_expression "$cleanup_fold_expression"
    model_ids=$(read_job_ids "$model_file")
elif any_stage_enabled analysis twas; then
    require_dir_with_match "${output_dir}/EGRET_models/${tissue}/${egret_subdir}" '*.wgt.RDat' "EGRET model weights"
fi

analysis_ids=""
if stage_enabled analysis; then
    analysis_upstream_ids=()
    [[ -z "$cis_ids" ]] || analysis_upstream_ids+=("$cis_ids")
    [[ -z "$model_ids" ]] || analysis_upstream_ids+=("$model_ids")
    analysis_deps=""
    ((${#analysis_upstream_ids[@]} == 0)) || analysis_deps=$(join_by_colon "${analysis_upstream_ids[@]}")
    jid=$(submit_stage "${tissue}_analysis" "$analysis_deps" EGRET_analysis "${scripts_dir}/7.0_run_model_analysis.sh" \
        --tissue "$tissue" --egret_subdir "$egret_subdir" --output_dir "$output_dir" --conda_env "$conda_env" \
        --cleanup_egret_plink "$cleanup_egret_plink" --cleanup_egret_working "$cleanup_egret_working" \
        --cleanup_egret_z_matrices "$cleanup_egret_z_matrices" --cleanup_cis_working "$cleanup_cis_working")
    record_single_job "${job_id_dir}/analysis.jids" "$jid"; analysis_ids=$jid
else
    if any_stage_enabled twas grn; then
        require_file "${output_dir}/results_sumstats/${tissue}/${egret_subdir}.txt" "EGRET analysis summary"
    fi
    if any_stage_enabled twas_cis grn; then
        require_file "${output_dir}/results_sumstats/${tissue}/cis.txt" "cis analysis summary"
    fi
fi

twas_driver_ids=""
if stage_enabled twas; then
    twas_deps=$analysis_ids
    [[ -n "$twas_deps" ]] || twas_deps=$model_ids
    jid=$(submit_stage "${tissue}_twas_driver" "$twas_deps" TWAS_driver "${scripts_dir}/8.0_run_TWAS_EGRET.sh" \
        --tissue "$tissue" --output_dir "$output_dir" --gene_info "$gene_info" \
        --gwas_sumstat_dir "$gwas_sumstat_dir" --ld_ref "$ld_ref" \
        --egret_subdir "$egret_subdir" --gwas_traits "$gwas_traits" --conda_env "$conda_env" \
        --slurm_account "$slurm_account" --slurm_partition "$slurm_partition" --slurm_log_dir "$run_log_dir" \
        --job_id_file "${job_id_dir}/twas.jids")
    record_single_job "${job_id_dir}/twas_driver.jids" "$jid"; twas_driver_ids=$jid
elif stage_enabled grn; then
    require_dir_with_match "${output_dir}/TWAS_association_results/${tissue}/${egret_subdir}" '*.dat' "EGRET TWAS results"
fi

twas_cis_ids=""
if stage_enabled twas_cis; then
    twas_cis_deps=$analysis_ids
    [[ -n "$twas_cis_deps" ]] || twas_cis_deps=$cis_ids
    jid=$(submit_stage "${tissue}_cis_twas_driver" "$twas_cis_deps" TWAS_cis_driver "${scripts_dir}/8.1_run_TWAS_cis.sh" \
        --tissue "$tissue" --output_dir "$output_dir" --gene_info "$gene_info" \
        --gwas_sumstat_dir "$gwas_sumstat_dir" --ld_ref "$ld_ref" --gwas_traits "$gwas_traits" \
        --conda_env "$conda_env" --slurm_account "$slurm_account" \
        --slurm_partition "$slurm_partition" --slurm_log_dir "$run_log_dir" \
        --job_id_file "${job_id_dir}/twas_cis.jids")
    record_single_job "${job_id_dir}/twas_cis_driver.jids" "$jid"; twas_cis_ids=$jid
elif stage_enabled grn; then
    require_dir_with_match "${output_dir}/TWAS_association_results/${tissue}/cis" '*.dat' "cis TWAS results"
fi

# The TWAS jobs are dispatchers: they finish once their per-trait jobs are submitted, so
# depending on them is not enough. 9.0 runs after both, reads their child-ID files, and
# submits the GRN worker afterok on every child. An omitted arm supplies no child-ID file.
if stage_enabled grn; then
    # An array, not `printf | grep -v '^$'`: with both IDs empty that grep exits 1, and under
    # `set -euo pipefail` the same spelling in an assignment kills the script with no output.
    grn_dep_ids=()
    [[ -z "$analysis_ids" ]] || grn_dep_ids+=("$analysis_ids")
    [[ -z "$twas_driver_ids" ]] || grn_dep_ids+=("$twas_driver_ids")
    [[ -z "$twas_cis_ids" ]] || grn_dep_ids+=("$twas_cis_ids")
    grn_deps=""
    ((${#grn_dep_ids[@]} == 0)) || grn_deps=$(join_by_colon "${grn_dep_ids[@]}")
    twas_job_id_file=""; cis_twas_job_id_file=""
    [[ -z "$twas_driver_ids" ]] || twas_job_id_file="${job_id_dir}/twas.jids"
    [[ -z "$twas_cis_ids" ]] || cis_twas_job_id_file="${job_id_dir}/twas_cis.jids"
    jid=$(submit_stage "${tissue}_grn_driver" "$grn_deps" GRN_driver "${scripts_dir}/9.0_submit_GRN_analysis.sh" \
        --tissue "$tissue" --output_dir "$output_dir" --egret_subdir "$egret_subdir" \
        --gene_info "$gene_info" --eqtl_file "$eqtl_file" \
        --gwas_sumstat_dir "$gwas_sumstat_dir" --gwas_traits "$gwas_traits" \
        --pval_threshold "$grn_pval_threshold" \
        --conda_env "$conda_env" --slurm_account "$slurm_account" \
        --slurm_partition "$slurm_partition" --slurm_log_dir "$run_log_dir" \
        --twas_job_id_file "$twas_job_id_file" --cis_twas_job_id_file "$cis_twas_job_id_file" \
        --job_id_file "${job_id_dir}/grn.jids")
    record_single_job "${job_id_dir}/grn_driver.jids" "$jid"
fi

echo "All selected stages submitted for run $run_id"
