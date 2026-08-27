#!/bin/bash
#SBATCH --job-name="Fusion"
#SBATCH --output="batch_submissions/gtex_xtune.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=12G
#SBATCH --account=csd832
#SBATCH --export=ALL
#SBATCH -t 24:00:00
# --- 6.3_train_EGRET_chunk.sh ----------------------------------------------
# Trains EGRET models for one chunk of genes. Submitted once per chunk by 6.0.
# Per-gene failures are contained the same way as in 2.2.
#
#   in : bed_files/ and plink_results/ for {egret_subdir}; z_matrices/
#   out: EGRET_models/{tissue}/{egret_subdir}/{gene}.wgt.RDat
#   run: submitted by 6.0_train_EGRET_models.sh
# --------------------------------------------------------------------------
set -euo pipefail
ulimit -c 0
# Locate workflow_lib.sh; see its header for why BASH_SOURCE alone is not enough.
scripts_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd || true)
[[ -r "${scripts_dir}/workflow_lib.sh" ]] || scripts_dir=${EGRET_SCRIPTS_DIR:-${SLURM_SUBMIT_DIR:-}}
[[ -r "${scripts_dir}/workflow_lib.sh" ]] || { echo "ERROR: cannot locate workflow_lib.sh (tried BASH_SOURCE, EGRET_SCRIPTS_DIR, SLURM_SUBMIT_DIR); export EGRET_SCRIPTS_DIR to the EGRET scripts directory" >&2; exit 1; }
export EGRET_SCRIPTS_DIR="$scripts_dir"
source "${scripts_dir}/workflow_lib.sh"
USAGE="Usage: 6.3_train_EGRET_chunk.sh --tissue NAME --chunk_size N --start_index N --output_dir DIR --egret_subdir NAME --plink_path PATH --gemma_path PATH --expression FILE --gene_info FILE --conda_env NAME [options]"
CLI_FLAGS=(--tissue --chunk_size --start_index --output_dir --egret_subdir --folds --plink_path --gemma_path --expression --fusion_models --gene_info --conda_env --cleanup_egret_plink --cleanup_egret_working --cleanup_egret_z_matrices)
CLI_VARS=(tissue chunk_size start_index output_dir egret_subdir folds plink_path gemma_path expression fusion_models gene_info conda_env cleanup_egret_plink cleanup_egret_working cleanup_egret_z_matrices)
CLI_REQUIRED=(tissue chunk_size start_index output_dir egret_subdir plink_path gemma_path expression gene_info conda_env)
folds=5; fusion_models=xtune,lasso,enet,blup; cleanup_egret_plink=FALSE; cleanup_egret_working=FALSE
cleanup_egret_z_matrices=FALSE
parse_cli "$@"; activate_workflow_env "$conda_env"
weights="${output_dir}/EGRET_models/${tissue}/${egret_subdir}"; wd="${output_dir}/working/${tissue}/${egret_subdir}"
z_matrix_dir="${output_dir}/z_matrices/${tissue}/${egret_subdir}"; plink_dir="${output_dir}/plink_results/${tissue}/${egret_subdir}"
mkdir -p "$weights" "$wd"
individuals="${output_dir}/fold_0_info/${tissue}/train_individuals.txt"
expression_index_file=$(mktemp "${wd}/egret_expression_index.${start_index}.XXXXXX")
trap 'rm -f "$expression_index_file"' EXIT
first_expression_row=$((start_index + 2))
last_expression_row=$((start_index + chunk_size + 1))
zcat "$expression" | awk -F '\t' -v OFS='\t' -v inds_file="$individuals" \
    -v first="$first_expression_row" -v last="$last_expression_row" '
    BEGIN { while ((getline id < inds_file) > 0) keep[id] = 1 }
    NR == 1 {
        for (i = 1; i <= NF; i++) if ($i in keep) columns[++n_columns] = i
        next
    }
    NR >= first && NR <= last {
        printf "%s", $1
        for (i = 1; i <= n_columns; i++) printf "%s%s", OFS, $(columns[i])
        printf "%s", ORS
    }
' > "$expression_index_file"
genes=()
declare -A expression_by_gene
while IFS= read -r expression_row; do
    gene=${expression_row%%$'\t'*}
    genes+=("$gene")
    expression_by_gene["$gene"]=${expression_row#*$'\t'}
done < "$expression_index_file"
attempted=0; failed=0
for gene in "${genes[@]}"; do
    [[ -n "$gene" ]] || continue
    [[ ! -s "${weights}/${gene}.wgt.RDat" ]] || { echo "Skipping $gene: weight exists"; continue; }
    skip=0
    for ((fold=0; fold<=folds; fold++)); do
        [[ -s "${plink_dir}/fold_${fold}/${gene}.bed" ]] || { echo "Skipping $gene: missing fold_${fold} PLINK input"; skip=1; break; }
        # A single-variant window makes FUSION's cross-validation fold collapse a
        # one-column matrix to a vector; see 2.2_compute_cis_weights_chunk.sh for the full reasoning.
        (( $(wc -l < "${plink_dir}/fold_${fold}/${gene}.bim") >= 2 )) || { echo "Skipping $gene: <2 variants in fold_${fold}"; skip=1; break; }
    done
    ((skip == 0)) || continue
    attempted=$((attempted + 1))
    individuals_plink="${wd}/${gene}_keep.txt"; awk '{print "0\t"$1}' "$individuals" > "$individuals_plink"
    gene_expression=${expression_by_gene[$gene]}
    # Guarded like the Rscript call below. Nested inside the per-gene loop, so a failure
    # sets a flag and breaks -- `continue` would only advance the fold.
    prep_ok=1
    for ((fold=0; fold<=folds; fold++)); do
        mkdir -p "${wd}/fold_${fold}"
        if ! "$plink_path" --bfile "${plink_dir}/fold_${fold}/${gene}" --make-bed --keep "$individuals_plink" \
            --indiv-sort file "$individuals_plink" --out "${wd}/fold_${fold}/${gene}"; then
            echo "WARNING: plink failed for ${gene} fold_${fold}; skipping gene" >&2
            prep_ok=0; break
        fi
        if ! { paste --delimiters='\t' <(cut -f1-5 "${wd}/fold_${fold}/${gene}.fam") \
            <(printf '%s\n' "$gene_expression" | tr '\t' '\n') > "${wd}/fold_${fold}/${gene}.mod.fam" \
            && mv "${wd}/fold_${fold}/${gene}.mod.fam" "${wd}/fold_${fold}/${gene}.fam"; }; then
            echo "WARNING: phenotype write failed for ${gene} fold_${fold}; skipping gene" >&2
            prep_ok=0; break
        fi
    done
    if (( prep_ok == 0 )); then failed=$((failed + 1)); continue; fi
    # Per-gene failures must not abort the chunk -- see 2.2_compute_cis_weights_chunk.sh.
    if ! Rscript "${scripts_dir}/6.4_train_EGRET_model_weights.R" --gene "$gene" --working_dir "$wd" \
        --fusion_models "$fusion_models" --weights_dir "$weights" --tissue "$tissue" --plink_path "$plink_path" \
        --gemma_path "$gemma_path" --z_matrix_dir "$z_matrix_dir" --gene_info "$gene_info" \
        --folds "$folds" --output_dir "$output_dir"; then
        echo "WARNING: 6.4_train_EGRET_model_weights.R failed for ${gene}; continuing" >&2
        failed=$((failed + 1))
    fi
    if [[ -f "${weights}/${gene}.wgt.RDat" ]]; then
        # Cleanup: the per-gene plink inputs produced by the bed stage. Kept on failure so
        # a retry does not have to re-run 6.2_make_EGRET_beds_by_fold.R for this gene.
        if [[ "$cleanup_egret_plink" == TRUE ]]; then
            for ((fold=0; fold<=folds; fold++)); do
                rm -f "${plink_dir}/fold_${fold}/${gene}.bed" "${plink_dir}/fold_${fold}/${gene}.bim" \
                    "${plink_dir}/fold_${fold}/${gene}.fam" "${plink_dir}/fold_${fold}/${gene}.log"
            done
        fi
        # Cleanup: z-matrices are written by 6.2 and read exactly once, by 6.4 above.
        # Success-gated like the other hooks; 6.2 regenerates them on a retry.
        if [[ "$cleanup_egret_z_matrices" == TRUE ]]; then
            for ((fold=0; fold<=folds; fold++)); do
                rm -f "${z_matrix_dir}/fold_${fold}/${gene}_z_matrix.txt"
            done
        fi
    fi
    # Cleanup: working/ is pure scratch -- plink re-writes of files still in plink_results/,
    # plus GEMMA output. Input to nothing and rebuilt in seconds, so it is swept whether or
    # not the gene trained; success-gating this is what makes a tree grow without bound.
    if [[ "$cleanup_egret_working" == TRUE ]]; then
        for ((fold=0; fold<=folds; fold++)); do
            rm -rf "${wd}/fold_${fold}/${gene}."*
        done
        rm -f "$individuals_plink"
    fi
done
if (( attempted > 0 && failed == attempted )); then
    cli_die "all ${attempted} genes in this chunk failed"
fi
echo "EGRET chunk complete: ${attempted} genes attempted, ${failed} failed"
