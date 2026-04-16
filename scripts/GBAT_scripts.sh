#!/bin/bash

tissue=$1                   # tissue for which the analysis will be done on
FDR=$2                      # FDR threshold at which to select variants
folds=$3                    # number of folds
gene_info=$4                # gene info file containing gene id, chr, start location
cis_model_dir=$5            # directory containing pretrained cis models with lasso as a model name
plink_path=$6               # path to plink executable
genotype_file_path=$7       # path to genotype files to use for imputation of cis models
output_dir=$8               # output directory where results will be stored
scripts_dir=$9              # directory where scripts are located

if false ; then

Rscript ${scripts_dir}/5_run_GBAT_updated.R \
    --tissue $tissue \
    --gene_info $gene_info \
    --output_dir $output_dir \
    --cis_model_dir $cis_model_dir \
    --genotype_file_path $genotype_file_path \
    --plink_path $plink_path


Rscript ${scripts_dir}/5.1_run_GBAT_association.R \
    --tissue $tissue \
    --output_dir $output_dir

fi

Rscript ${scripts_dir}/5.3_make_GBAT_bed_by_FDR.R \
    --tissue $tissue \
    --gene_info $gene_info \
    --cis_model_dir $cis_model_dir \
    --FDR $FDR \
    --output_dir $output_dir