#!/bin/bash
#SBATCH --job-name="transPCO_results"
#SBATCH --output="batch_submissions/transPCO_results.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=12G
#SBATCH --account=csd832
#SBATCH --export=ALL
#SBATCH -t 02:00:00

source ~/.bashrc
conda activate transTWAS_env

tissue=$1
pval_threshold=$2
output_dir=$3
folds=$4
expected_module_count=$5
scripts_dir=$6
genotype_prefix=$7

# Verify PCO output files exist before running results analysis
total_found=0
for fold in $(seq 0 $folds); do
    pco_dir="${output_dir}/transPCO/${tissue}/fold_${fold}/PCO_association_results"
    if [ -d "$pco_dir" ]; then
        found=$(ls "${pco_dir}"/*.txt.gz 2>/dev/null | wc -l)
        total_found=$((total_found + found))
    fi
done

if [ "$total_found" -eq 0 ]; then
    echo "ERROR: No PCO output files found across any fold. Exiting."
    exit 1
fi

# Check that each fold has output for every expected module (22 chromosomes each)
expected_per_fold=$((expected_module_count * 22))
for fold in $(seq 0 $folds); do
    pco_dir="${output_dir}/transPCO/${tissue}/fold_${fold}/PCO_association_results"
    found=$(ls "${pco_dir}"/*.txt.gz 2>/dev/null | wc -l)
    if [ "$found" -lt "$expected_per_fold" ]; then
        echo "WARNING: fold_${fold} has ${found} PCO files, expected ${expected_per_fold} (${expected_module_count} modules x 22 chromosomes)"
    fi
done
echo "Found ${total_found} total PCO output files across all folds"

Rscript ${scripts_dir}/11_transPCO_results_analysis_by_FDR.R \
    --tissue $tissue \
    --FDR $pval_threshold \
    --output_dir $output_dir \
    --folds $folds \
    --PCO_association_dir PCO_association_results \
    --module_dir modules \
    --genotype_prefix $genotype_prefix
