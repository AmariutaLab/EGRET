#!/bin/bash

tissue=$1                   # tissue for which the analysis will be done on
FDR=$2                      # FDR threshold at which to select variants
folds=$3                    # number of folds
gene_info=$4                # gene info file containing gene id, chr, start location
genotype_output_prefix=$5   # to access the bim file for snp location
output_dir=$6               # output directory where results will be stored
scripts_dir=$7              # directory where scripts are located


Rscript ${scripts_dir}/7_run_MatrixeQTL.R \
    --tissue ${tissue} \
    --folds ${folds} \
    --output_dir ${output_dir}


Rscript ${scripts_dir}/8_make_MatrixeQTL_bed_by_FDR.R \
    --tissue ${tissue} \
    --folds ${folds} \
    --FDR ${FDR} \
    --gene_info ${gene_info} \
    --genotype_output_prefix ${genotype_output_prefix} \
    --output_dir ${output_dir}