#!/bin/bash
tissue=$1
folds=$2
models=$3
genotypes_file_path=$4
expression_file_path=$5
plink_path=${6:-"plink2"}
output_dir=$7
MatrixeQTL_bed_dir=${8:-"results_FDR_0.1"}
GBAT_bed_dir=${9:-"results_FDR_0.1"}
transPCO_bed_dir=${10:-"bed_files_FDR_0.1"}
exclude_crossmap=${11:-"TRUE"}
crossmap_dir=${12:-"background_mismatches"}
genotype_prefix=${13:-"GTEX_v8_genotypes_pruned"}
egret_output_subdir=${14:-"EGRET"}
fusion_models=${15:-"xtune,lasso,enet,blup"}
gemma_path=${16:-"../full_workflow/gemma-0.98.5-linux-static-AMD64"}
chunk_size=${17:-500}
gene_info_file_path=${18:-"../../data/GTEx_V8.txt.gz"}

for fold in $(seq 0 $folds); do
    Rscript 12.4_make_xtune_bed_no_cross_mappable_by_fold.R \
        --tissue "$tissue" \
        --base_dir "$output_dir" \
        --MatrixeQTL_bed_dir "$MatrixeQTL_bed_dir" \
        --GBAT_bed_dir "$GBAT_bed_dir" \
        --transPCO_bed_dir "$transPCO_bed_dir" \
        --output_dir "$egret_output_subdir" \
        --exclude_crossmap "$exclude_crossmap" \
        --crossmap_dir "$crossmap_dir" \
        --genotype_prefix "$genotype_prefix" \
        --plink_path "$plink_path" \
        --expression_file "$expression_file_path" \
        --models Cis,MatrixeQTL,GBAT,transPCO \
        --fold "$fold"

done

plink_dir="${output_dir}/plink_results/${tissue}/${egret_output_subdir}/fold_0"
files=($(ls "$plink_dir" | grep bed))
echo ${#files[@]}

num_chunks=$(((${#files[@]} + $chunk_size - 1) / $chunk_size))
echo $num_chunks
# Loop through each chunk
for ((i = 0; i < num_chunks; i++)); do
        start_index=$((i * chunk_size))
        end_index=$((start_index + chunk_size - 1))
        sbatch 14_create_EGRET_models.sh "$tissue" "$chunk_size" "$end_index" \
            "$output_dir" "$egret_output_subdir" "$folds" "$plink_path" \
            "$gemma_path" "$expression_file_path" "$fusion_models" "$gene_info_file_path"
done
