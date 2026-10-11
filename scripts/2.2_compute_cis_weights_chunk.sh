#!/bin/bash
#SBATCH --job-name="Running_Fusion"
#SBATCH --output="batch_submissions/cis_fusion.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2G
#SBATCH --account=csd832
#SBATCH --export=ALL
#SBATCH -t 24:00:00
# --- 2.2_compute_cis_weights_chunk.sh --------------------------------------
# Trains cis FUSION weights for one chunk of genes. Submitted once per chunk by 2.0.
# Per-gene failures are contained: a bad gene is counted and skipped, and the chunk
# fails only if every gene in it failed.
#
#   in : plink_results/{tissue}/cis/{gene}.*; covariates
#   out: FUSION/{tissue}/cis/{gene}.wgt.RDat
#   run: submitted by 2.0_run_cis_models.sh
# --------------------------------------------------------------------------
set -euo pipefail
# Locate workflow_lib.sh; see its header for why BASH_SOURCE alone is not enough.
scripts_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd || true)
[[ -r "${scripts_dir}/workflow_lib.sh" ]] || scripts_dir=${EGRET_SCRIPTS_DIR:-${SLURM_SUBMIT_DIR:-}}
[[ -r "${scripts_dir}/workflow_lib.sh" ]] || { echo "ERROR: cannot locate workflow_lib.sh (tried BASH_SOURCE, EGRET_SCRIPTS_DIR, SLURM_SUBMIT_DIR); export EGRET_SCRIPTS_DIR to the EGRET scripts directory" >&2; exit 1; }
export EGRET_SCRIPTS_DIR="$scripts_dir"
source "${scripts_dir}/workflow_lib.sh"
USAGE="Usage: 2.2_compute_cis_weights_chunk.sh --tissue NAME --output_dir DIR --chunk_size N --start_index N --plink_path PATH --individuals FILE --covariates FILE --plink1_path PATH --gcta_path PATH --gemma_path PATH --conda_env NAME [options]"
CLI_FLAGS=(--tissue --output_dir --chunk_size --start_index --plink_path --individuals --covariates --plink1_path --gcta_path --gemma_path --conda_env --cleanup_cis_working --hsq_set)
CLI_VARS=(tissue output_dir chunk_size start_index plink_path individuals covariates plink1_path gcta_path gemma_path conda_env cleanup_cis_working hsq_set)
CLI_REQUIRED=(tissue output_dir chunk_size start_index plink_path individuals covariates plink1_path gcta_path gemma_path conda_env)
cleanup_cis_working=FALSE; hsq_set=""
parse_cli "$@"; activate_workflow_env "$conda_env"
wd="${output_dir}/working/${tissue}/cis"; tmpdir="${output_dir}/tmp/${tissue}/cis"; weights="${output_dir}/FUSION/${tissue}/cis"
mkdir -p "$wd" "$tmpdir" "$weights"
cis_manifest="${output_dir}/plink_results/${tissue}/cis/cis_genes.txt"
require_file "$cis_manifest" "cis genotype manifest"
mapfile -t all_genes < "$cis_manifest"
genes=("${all_genes[@]:start_index:chunk_size}")
(( ${#genes[@]} > 0 )) || exit 0
gefile="${output_dir}/expression_files/${tissue}_expression.txt.gz"
chunk_gene_file=$(mktemp "${wd}/cis_chunk_genes.${start_index}.XXXXXX")
expression_index_file=$(mktemp "${wd}/cis_expression_index.${start_index}.XXXXXX")
trap 'rm -f "$chunk_gene_file" "$expression_index_file"' EXIT
for gene_file in "${genes[@]}"; do
    printf '%s\n' "${gene_file%.bed}"
done > "$chunk_gene_file"
awk -F '\t' -v OFS='\t' -v inds_file="$individuals" '
    BEGIN { while ((getline id < inds_file) > 0) keep[id] = 1 }
    FNR == NR { wanted[$1] = 1; next }
    FNR == 1 {
        for (i = 1; i <= NF; i++) if ($i in keep) columns[++n_columns] = i
        next
    }
    $1 in wanted {
        printf "%s", $1
        for (i = 1; i <= n_columns; i++) printf "%s%s", OFS, $(columns[i])
        printf "%s", ORS
    }
' "$chunk_gene_file" <(zcat "$gefile") > "$expression_index_file"
declare -A expression_by_gene
while IFS= read -r expression_row; do
    expression_gene=${expression_row%%$'\t'*}
    expression_by_gene["$expression_gene"]=${expression_row#*$'\t'}
done < "$expression_index_file"
attempted=0; failed=0
for gene_file in "${genes[@]}"; do
    gene=${gene_file%.bed}
    # Resume guard: a re-run skips genes already trained. Genes FUSION dropped for low
    # heritability write no .wgt.RDat and are retried, which is cheap and correct.
    if [[ -s "${weights}/${gene}.wgt.RDat" ]]; then
        echo "Skipping $gene: weight exists"
        continue
    fi
    attempted=$((attempted + 1))
    individuals_plink="${wd}/${gene}_keep.txt"; awk '{print "0\t"$1}' "$individuals" > "$individuals_plink"
    # Guarded: under `set -e` one gene's failure would abort the chunk and drop the rest.
    if ! "$plink_path" --bfile "${output_dir}/plink_results/${tissue}/cis/${gene}" --make-bed --keep "$individuals_plink" --indiv-sort file "$individuals_plink" --out "${wd}/${gene}"; then
        echo "WARNING: plink failed for ${gene}; continuing" >&2
        failed=$((failed + 1)); continue
    fi
    # The gene list comes from a directory listing, so a gene may be absent from the
    # expression file. Without this, a blank phenotype column would be pasted into the .fam.
    [[ -v expression_by_gene[$gene] ]] || cli_die "gene $gene is absent from $gefile"
    ge_donors=${expression_by_gene[$gene]}
    # Guarded as a unit: this consumes plink's .fam, so it fails if that .fam is malformed.
    if ! { paste --delimiters='\t' <(cut -f1-5 "${wd}/${gene}.fam") <(printf '%s\n' "$ge_donors" | tr '\t' '\n') > "${wd}/${gene}.mod.fam" && mv "${wd}/${gene}.mod.fam" "${wd}/${gene}.fam"; }; then
        echo "WARNING: phenotype write failed for ${gene}; continuing" >&2
        failed=$((failed + 1)); continue
    fi
    # FUSION's cross-validation does `genos$bed[i, ] %*% pred.wgt`, which fails with
    # "non-conformable arguments" when the cis window holds a single variant: R drops the
    # one-column matrix to a vector. A one-SNP model is worthless anyway, so skip it here
    # rather than patching vendored FUSION code.
    nvar=$(wc -l < "${wd}/${gene}.bim")
    if (( nvar < 2 )); then
        echo "Skipping $gene: only ${nvar} variant(s) in the cis window"
        continue
    fi
    # Guarded: under `set -e` one gene's failure would abort the chunk and drop the rest.
    # --hsq_set skips GCTA and fixes h2 (used only as the lasso penalty), so genes whose
    # REML fails to converge or goes negative still get a model.
    if ! Rscript "${scripts_dir}/FUSION.compute_weights.R" --genotype_prefix "${wd}/${gene}" --crossval 5 \
        --fusion_models lasso,blup,enet,top1 --hsq_p 1 --tmp "${tmpdir}/${gene}" --out "${weights}/${gene}" \
        --covariates "$covariates" --gcta_path "$gcta_path" --plink_path "$plink1_path" --gemma_path "$gemma_path" \
        ${hsq_set:+--hsq_set "$hsq_set"}; then
        echo "WARNING: FUSION.compute_weights.R failed for ${gene}; continuing" >&2
        failed=$((failed + 1))
    fi
    if [[ "$cleanup_cis_working" == TRUE && -f "${weights}/${gene}.wgt.RDat" ]]; then
        rm -rf "${wd}/${gene}."* "$individuals_plink" "${tmpdir}/${gene}"
    fi
done
# A chunk in which every single gene failed is a real problem (bad inputs, missing binary);
# a chunk with some failures is normal and must not block the dependent stages.
if (( attempted > 0 && failed == attempted )); then
    cli_die "all ${attempted} genes in this chunk failed"
fi
echo "cis chunk complete: ${attempted} genes attempted, ${failed} failed"
