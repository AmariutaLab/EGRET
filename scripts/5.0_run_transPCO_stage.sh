#!/bin/bash
tissue=$1
folds=$2
FDR=$3
covariate_columns_for_coexpression=$4
gene_info=$5
plink_path=$6
output_dir=$7
genotype_prefix=${8:-"GTEX_v8_genotypes_pruned"}
scripts_dir=${9}
slurm_log_dir=${10:-""}
cleanup_transpco_raw=${11:-"FALSE"}


Rscript ${scripts_dir}/9.0_run_WGCNA_clustering.R \
	--expression ${output_dir}/expression_files/${tissue}_expression.txt.gz \
	--output_dir ${output_dir}/ \
	--individuals ${output_dir}/fold_0_info/${tissue}/train_individuals.txt \
	--covariates ${output_dir}/covariate_files/${tissue}_covariates.txt.gz \
	--tissue $tissue \
    --gene_info $gene_info \
	--covariate_columns_for_coexpression "$covariate_columns_for_coexpression"


Rscript ${scripts_dir}/9.1_initialize_chr_genotypes.R \
    --tissue $tissue \
    --output_dir $output_dir \
    --folds $folds \
    --genotype_prefix $genotype_prefix \
    --plink_path $plink_path



# Find the fold with the most modules (use as the canonical module list)
max_files=0
for fold in $(seq 0 $folds); do
    dir="${output_dir}/transPCO/${tissue}/fold_${fold}/modules"
    file_count=$(find "$dir" -type f 2>/dev/null | wc -l)
    if (( file_count > max_files )); then
        max_files=$file_count
        max_module_dir=$dir
    fi
done
echo "Using module list from: $max_module_dir"

job_ids=()
for file in "$max_module_dir"/*; do
    jid=$(sbatch --parsable ${slurm_log_dir:+--output=${slurm_log_dir}/transPCO_matrixeQTL.%j.%N.out} ${scripts_dir}/run_transPCO_matrixeQTL_job.sh \
        $tissue \
        $(basename $file) \
        modules \
        $output_dir \
        $folds \
        $plink_path \
        $genotype_prefix \
        $scripts_dir)
    job_ids+=($jid)
    echo "Submitted job $jid for module $(basename $file)"
done
echo "Submitted ${#job_ids[@]} module jobs"

# Build dependency string from all MatrixEQTL job IDs
if [ ${#job_ids[@]} -gt 0 ]; then
    dep_str=$(IFS=:; echo "${job_ids[*]}")
    dependency="--dependency=afterok:${dep_str}"
else
    dependency=""
fi


# Submit one PCO association job per module, dependent on all MatrixEQTL jobs
modules=$(find ${output_dir}/transPCO/${tissue}/fold_*/modules/ -type f -exec basename {} \; | sort | uniq)

pco_job_ids=()
for module in $modules
do
    jid=$(sbatch --parsable $dependency ${slurm_log_dir:+--output=${slurm_log_dir}/transPCO_PCO_assoc.%j.%N.out} ${scripts_dir}/run_transPCO_PCO_association_job.sh \
        $tissue \
        $module \
        modules \
        $output_dir \
        $folds \
        $gene_info \
        $scripts_dir)
    pco_job_ids+=($jid)
    echo "Submitted PCO job $jid for module $module"
done
echo "Submitted ${#pco_job_ids[@]} PCO association jobs"


# Count expected modules for verification
module_count=$(echo "$modules" | wc -w)

# Submit results analysis after all PCO jobs finish (afterany so it can check file counts itself)
if [ ${#pco_job_ids[@]} -gt 0 ]; then
    pco_dep_str=$(IFS=:; echo "${pco_job_ids[*]}")
    pco_dependency="--dependency=afterany:${pco_dep_str}"
else
    pco_dependency=""
fi

sbatch $pco_dependency ${slurm_log_dir:+--output=${slurm_log_dir}/transPCO_results.%j.%N.out} ${scripts_dir}/run_transPCO_results_analysis_job.sh \
    $tissue \
    $FDR \
    $output_dir \
    $folds \
    $module_count \
    $scripts_dir \
    $genotype_prefix \
    "$cleanup_transpco_raw"
