#!/bin/bash

genotypes_file_path=$1     # path to bed/bim/fam file containing genotype information for all individuals
expression_file_path=$2    # path to file containing expression with column containing gene_id as ENSG format and other columns labeled by their individual ids
covariates_file_path=$3    # path to file containing covariate information with rows as given by their individial ids and columns as covariates to be regressed out
tissue=$4                  # tissue for analysis
individuals_file_path=$5   # path to file with no header containing a list of individual ids for analysis
LD_prune_r2=$6             # LD r2 limit by which to prune variants by. Default is 0.9
plink_path=$7              # path to plink executable
genotype_output_prefix=$8  # prefix which to store genotypes by
folds=$9                   # number of folds for crossvalidation
gene_info_file_path=${10}  # file containing gene info such as gene id, gene name, chr, and start location
output_dir=${11}           # directory where output files will be stored

Rscript 0_setup_genotypes.R \
	--plink_path $plink_path \
	--bfile $genotypes_file_path \
	--LD_r2 $LD_prune_r2 \
	--LD_chunk 1 \
	--LD_window 100 \
	--out ${output_dir}/genotype_files/${genotype_output_prefix}

Rscript 1_setup_expression.R \
	--expression $expression_file_path \
	--individuals  $individuals_file_path \
	--covariates  $covariates_file_path \
	--tissue $tissue \
	--gene_info $gene_info_file_path \
	--output_dir $output_dir

Rscript 2_setup_folds.R \
	--expression $output_dir/expression_files/${tissue}_expression_regressed.txt.gz \
	--folds 5 \
	--tissue ${tissue} \
	--individuals $individuals_file_path \
	--plink_path $plink_path \
	--bfile $output_dir/genotype_files/${genotype_output_prefix} \
	--output_dir $output_dir
