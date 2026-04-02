#!/bin/bash
tissue=$1
output_dir=$2
gene_info_file_path=$3
gwas_sumstat_dir=$4
ld_ref=$5
egret_output_subdir=$6
gwas_sumstats_csv=$7

# Create pos files
Rscript 16.5_create_EGRET_pos.R \
    --tissue $tissue \
    --output_dir $output_dir \
    --gene_info_file_path $gene_info_file_path

# Set pos file and weights directory
wgt=${output_dir}/pos_files/${tissue}/cis_MatrixeQTL_GBAT_transPCO_FDR_0.1.pos
wgtdir=${output_dir}/xtune_fusion_models/${tissue}/cis_MatrixeQTL_GBAT_transPCO_FDR_0.1/

# Output directory for TWAS results
out_dir=${output_dir}/TWAS_association_results/${tissue}/${egret_output_subdir}

# Split comma-separated sumstats into array
IFS=',' read -ra sumstats_array <<< "$gwas_sumstats_csv"

for gwas_file in "${sumstats_array[@]}"; do
    gwas_sumstat_path=${gwas_sumstat_dir}/${gwas_file}
    trait="${gwas_file%%.sumstats.gz}"

    echo "Submitting TWAS job for ${trait}"
    sbatch 17_run_TWAS_association.sh \
        $tissue \
        $gwas_sumstat_path \
        $trait \
        $wgt \
        $wgtdir \
        $ld_ref \
        $out_dir
done
