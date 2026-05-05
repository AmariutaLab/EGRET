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
cleanup_cis_plink=${10:-"FALSE"}  # delete plink_results/{tissue}/cis/ after GBAT (last consumer)


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

Rscript ${scripts_dir}/5.3_make_GBAT_bed_by_FDR.R \
    --tissue $tissue \
    --gene_info $gene_info \
    --cis_model_dir $cis_model_dir \
    --FDR $FDR \
    --output_dir $output_dir

# Hook B: GBAT is the last consumer of plink_results/{tissue}/cis/. Drop it now.
if [[ "$cleanup_cis_plink" == "TRUE" ]]; then
    cis_plink_dir="${output_dir}/plink_results/${tissue}/cis"
    if [ -d "$cis_plink_dir" ]; then
        echo "Removing $cis_plink_dir (cleanup_cis_plink=TRUE)"
        rm -rf "$cis_plink_dir"
    fi
fi