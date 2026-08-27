#!/bin/bash
#SBATCH --job-name="EGRET_analysis"
#SBATCH --output="batch_submissions/EGRET_analysis.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16G
#SBATCH --account=csd832
#SBATCH --export=ALL
#SBATCH -t 04:00:00
# --- 7.0_run_model_analysis.sh ---------------------------------------------
# Stage `analysis`. Summarises both arms into the sumstats the TWAS and GRN stages
# read, then runs the terminal cleanup sweep over the scratch trees.
#
#   in : EGRET_models/{tissue}/{egret_subdir}/; FUSION/{tissue}/cis/
#   out: results_sumstats/{tissue}/{egret_subdir}.txt and cis.txt
#   run: run_EGRET_workflow.sh, stage `analysis`
# --------------------------------------------------------------------------
set -euo pipefail
# Locate workflow_lib.sh; see its header for why BASH_SOURCE alone is not enough.
scripts_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd || true)
[[ -r "${scripts_dir}/workflow_lib.sh" ]] || scripts_dir=${EGRET_SCRIPTS_DIR:-${SLURM_SUBMIT_DIR:-}}
[[ -r "${scripts_dir}/workflow_lib.sh" ]] || { echo "ERROR: cannot locate workflow_lib.sh (tried BASH_SOURCE, EGRET_SCRIPTS_DIR, SLURM_SUBMIT_DIR); export EGRET_SCRIPTS_DIR to the EGRET scripts directory" >&2; exit 1; }
export EGRET_SCRIPTS_DIR="$scripts_dir"
source "${scripts_dir}/workflow_lib.sh"
USAGE="Usage: 7.0_run_model_analysis.sh --tissue NAME --egret_subdir NAME --output_dir DIR --conda_env NAME [--cleanup_egret_plink BOOL] [--cleanup_egret_working BOOL]"
CLI_FLAGS=(--tissue --egret_subdir --output_dir --conda_env --cleanup_egret_plink --cleanup_egret_working --cleanup_egret_z_matrices --cleanup_cis_working)
CLI_VARS=(tissue egret_subdir output_dir conda_env cleanup_egret_plink cleanup_egret_working cleanup_egret_z_matrices cleanup_cis_working)
CLI_REQUIRED=(tissue egret_subdir output_dir conda_env)
cleanup_egret_plink=FALSE; cleanup_egret_working=FALSE; cleanup_egret_z_matrices=FALSE; cleanup_cis_working=FALSE
parse_cli "$@"; activate_workflow_env "$conda_env"
Rscript "${scripts_dir}/7.1_summarize_EGRET_models.R" --tissue "$tissue" \
    --egret_subdir "$egret_subdir" --output_dir "$output_dir"
cis_weight_dir="${output_dir}/FUSION/${tissue}/cis"
if [[ -d "$cis_weight_dir" ]]; then
    Rscript "${scripts_dir}/7.2_summarize_cis_models.R" --weight_dir "$cis_weight_dir" \
        --output_file cis.txt --tissue "$tissue" --output_dir "$output_dir"
else
    cli_die "cis weight directory is missing: $cis_weight_dir"
fi

# Terminal cleanup sweep. The per-gene hooks in 2.2 and 6.3 are success-gated, so failed
# genes keep their scratch; this job waits on every chunk, so by now it is all dead weight.
# Each sweep is opt-in. A run that stops before this stage still leaks.
if [[ "$cleanup_egret_working" == TRUE ]]; then
    working_dir="${output_dir}/working/${tissue}/${egret_subdir}"
    [[ ! -d "$working_dir" ]] || rm -rf "$working_dir"
fi
if [[ "$cleanup_egret_plink" == TRUE ]]; then
    egret_plink_dir="${output_dir}/plink_results/${tissue}/${egret_subdir}"
    [[ ! -d "$egret_plink_dir" ]] || rm -rf "$egret_plink_dir"
fi
if [[ "$cleanup_egret_z_matrices" == TRUE ]]; then
    egret_z_dir="${output_dir}/z_matrices/${tissue}/${egret_subdir}"
    [[ ! -d "$egret_z_dir" ]] || rm -rf "$egret_z_dir"
fi
# The cis stage has no sweep of its own; its weights live in FUSION/{tissue}/cis, so both
# of these trees are pure scratch by now.
if [[ "$cleanup_cis_working" == TRUE ]]; then
    cis_working_dir="${output_dir}/working/${tissue}/cis"
    [[ ! -d "$cis_working_dir" ]] || rm -rf "$cis_working_dir"
    cis_tmp_dir="${output_dir}/tmp/${tissue}/cis"
    [[ ! -d "$cis_tmp_dir" ]] || rm -rf "$cis_tmp_dir"
fi
