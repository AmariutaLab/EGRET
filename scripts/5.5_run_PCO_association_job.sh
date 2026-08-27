#!/bin/bash
#SBATCH --job-name="transPCO_PCO_assoc"
#SBATCH --output="batch_submissions/transPCO_PCO_assoc.%j.%N.out"
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
gene_info=$6
scripts_dir=$7

# Check that association result files exist before running PCO
module_name="${module%.*}"  # strip file extension

missing=0
has_any=0
for fold in $(seq 0 $folds); do
    module_file="${output_dir}/transPCO/${tissue}/fold_${fold}/${module_dir}/${module}"
    if [ ! -f "$module_file" ]; then
        echo "NOTICE: Module ${module_name} not present in fold_${fold} (expected for small modules)"
        continue
    fi
    assoc_dir="${output_dir}/transPCO/${tissue}/fold_${fold}/association_results"
    count=$(ls "${assoc_dir}"/matrix_${module_name}_chr_*.txt.gz 2>/dev/null | wc -l)
    if [ "$count" -eq 0 ]; then
        echo "ERROR: No association files found for ${module_name} in fold_${fold} at ${assoc_dir}"
        missing=1
    else
        has_any=1
    fi
done

if [ "$has_any" -eq 0 ] || [ "$missing" -eq 1 ]; then
    echo "Exiting: missing association result files for module ${module_name}"
    exit 1
fi

Rscript ${scripts_dir}/9.3_transPCO_multivariate_association_by_module.R \
    --tissue $tissue \
    --module $module \
    --module_dir $module_dir \
    --output_dir $output_dir \
    --folds $folds \
    --gene_info $gene_info \
    --results_dir PCO_association_results \
    --association_dir association_results \
    --scripts_dir $scripts_dir
