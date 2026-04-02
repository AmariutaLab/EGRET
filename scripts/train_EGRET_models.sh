#!/bin/bash
tissue=$1
folds=$2
models=$3
plink_path=${4:-"plink2"}
output_dir=$5
MatrixeQTL_bed_dir=${6:-"results_FDR_0.1"}
GBAT_bed_dir=${7:-"results_FDR_0.1"}
transPCO_bed_dir=${8:-"bed_files_FDR_0.1"}
exclude_crossmap=${9:-"TRUE"}
crossmap_dir=${10:-"background_mismatches"}
genotype_prefix=${11:-"GTEX_v8_genotypes_pruned"}
egret_output_subdir=${12:-"EGRET"}
fusion_models=${13:-"xtune,lasso,enet,blup"}
gemma_path=${14:-"./gemma-0.98.5-linux-static-AMD64"}
chunk_size=${15:-500}
gene_info_file_path=${16:-"../../data/GTEx_V8.txt.gz"}

expression_file_path=${output_dir}/expression_files/${tissue}_expression_regressed.txt.gz

if false ; then
bed_job_ids=()
for fold in $(seq 0 $folds); do
    jid=$(sbatch --parsable run_make_EGRET_bed_files.sh \
        "$tissue" \
        "$output_dir" \
        "$MatrixeQTL_bed_dir" \
        "$GBAT_bed_dir" \
        "$transPCO_bed_dir" \
        "$egret_output_subdir" \
        "$exclude_crossmap" \
        "$crossmap_dir" \
        "$genotype_prefix" \
        "$plink_path" \
        "$expression_file_path" \
        "$fold" \
        "$gene_info_file_path")
    bed_job_ids+=($jid)
    echo "Submitted bed file job $jid for fold $fold"
done
echo "Submitted ${#bed_job_ids[@]} bed file jobs"

# Build dependency string so model jobs wait for all bed file jobs
dep_str=$(IFS=:; echo "${bed_job_ids[*]}")
dependency="--dependency=afterok:${dep_str}"
fi

# Estimate gene count from expression file (upper bound)
num_genes=$(zcat "${output_dir}/expression_files/${tissue}_expression.txt.gz" | tail -n +2 | wc -l)
echo "Estimated gene count: $num_genes"

num_chunks=$(((num_genes + chunk_size - 1) / chunk_size))
echo "Submitting $num_chunks model chunks"

# Loop through each chunk
for ((i = 0; i < num_chunks; i++)); do
        start_index=$((i * chunk_size))
        end_index=$((start_index + chunk_size - 1))
        sbatch $dependency 14_create_EGRET_models.sh "$tissue" "$chunk_size" "$end_index" \
            "$output_dir" "$egret_output_subdir" "$folds" "$plink_path" \
            "$gemma_path" "$expression_file_path" "$fusion_models" "$gene_info_file_path"
done
