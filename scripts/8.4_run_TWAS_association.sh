#!/bin/bash
#SBATCH --job-name="Running_TWAS"
#SBATCH --output="batch_submissions/TWAS_output_%j.%N.out"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --account=csd832
#SBATCH --export=ALL
#SBATCH -t 48:00:00
# --- 8.4_run_TWAS_association.sh -------------------------------------------
# Runs the FUSION TWAS for one GWAS trait. Submitted once per trait by 8.0 / 8.1.
# Skips traits whose output already exists, so a re-run resumes.
#
#   in : pos_files/{tissue}/{arm}.pos; a GWAS sumstats file; ld_ref
#   out: TWAS_association_results/{tissue}/{arm}/{trait}.dat
#   run: submitted by 8.0_run_TWAS_EGRET.sh and 8.1_run_TWAS_cis.sh
# --------------------------------------------------------------------------
set -euo pipefail
# Locate workflow_lib.sh; see its header for why BASH_SOURCE alone is not enough.
scripts_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd || true)
[[ -r "${scripts_dir}/workflow_lib.sh" ]] || scripts_dir=${EGRET_SCRIPTS_DIR:-${SLURM_SUBMIT_DIR:-}}
[[ -r "${scripts_dir}/workflow_lib.sh" ]] || { echo "ERROR: cannot locate workflow_lib.sh (tried BASH_SOURCE, EGRET_SCRIPTS_DIR, SLURM_SUBMIT_DIR); export EGRET_SCRIPTS_DIR to the EGRET scripts directory" >&2; exit 1; }
export EGRET_SCRIPTS_DIR="$scripts_dir"
source "${scripts_dir}/workflow_lib.sh"
USAGE="Usage: 8.4_run_TWAS_association.sh --gwas_sumstats FILE --trait NAME --weights FILE --weights_dir DIR --ld_ref PREFIX --output_dir DIR --conda_env NAME"
CLI_FLAGS=(--gwas_sumstats --trait --weights --weights_dir --ld_ref --output_dir --conda_env)
CLI_VARS=(gwas_sumstats trait weights weights_dir ld_ref output_dir conda_env)
CLI_REQUIRED=(gwas_sumstats trait weights weights_dir ld_ref output_dir conda_env)
parse_cli "$@"; activate_workflow_env "$conda_env"
mkdir -p "$output_dir"
result="${output_dir}/${trait}.dat"
if [[ -f "$result" && $(wc -l < "$result") -gt 1 ]]; then echo "Output already exists for $trait"; exit 0; fi
Rscript "${scripts_dir}/FUSION.assoc_test_trans.R" --chr 1 --ref_ld_chr "$ld_ref" \
    --sumstats "$gwas_sumstats" --weights "$weights" --weights_dir "$weights_dir" --out "$result"
