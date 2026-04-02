#!/bin/bash
#SBATCH --job-name="make_EGRET_bed"
#SBATCH --output="batch_submissions/make_EGRET_bed.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=12G
#SBATCH --account=csd832
#SBATCH --export=ALL
#SBATCH -t 24:00:00

source ~/.bashrc
conda activate transTWAS_env

tissue=$1
output_dir=$2
MatrixeQTL_bed_dir=$3
GBAT_bed_dir=$4
transPCO_bed_dir=$5
egret_output_subdir=$6
exclude_crossmap=$7
crossmap_dir=$8
genotype_prefix=$9
plink_path=${10}
expression_file_path=${11}
fold=${12}
gene_info_file_path=${13}

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
    --fold "$fold" \
    --gene_information "$gene_info_file_path"
