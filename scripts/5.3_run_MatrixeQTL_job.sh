#!/bin/bash
#SBATCH --job-name="transPCO_matrixeQTL"
#SBATCH --output="batch_submissions/transPCO_matrixeQTL.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=12G
#SBATCH --account=csd832
#SBATCH --export=ALL
#SBATCH -t 12:00:00

source ~/.bashrc
conda activate transTWAS_env

tissue=$1
module=$2
module_dir=$3
output_dir=$4
folds=$5
plink_path=$6
genotype_prefix=${7:-"GTEX_v8_genotypes_pruned"}
scripts_dir=$8

Rscript ${scripts_dir}/9.10_cleaned_up_MatrixeQTL_for_transPCO.R \
    --tissue $tissue \
    --module $module \
    --module_dir $module_dir \
    --output_dir $output_dir \
    --folds $folds \
    --plink_path $plink_path \
    --genotype_prefix $genotype_prefix
