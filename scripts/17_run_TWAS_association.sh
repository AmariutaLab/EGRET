#!/bin/bash
#SBATCH --job-name="Running TWAS"
#SBATCH --output="/expanse/lustre/scratch/kbrunton/temp_project/batch_submission_output/TWAS_output/%j.%N.out"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --account=csd832
#SBATCH --export=ALL
#SBATCH -t 48:00:00

source ~/.bashrc
conda activate transTWAS_env

tissue=$1
gwas_sumstat_path=$2
trait=$3
wgt=$4
wgtdir=$5
ld_ref=$6
out_dir=$7
scripts_dir=$8

chr=1  # doesn't matter if using trans script

mkdir -p $out_dir

if [ -f "${out_dir}/${trait}.dat" ]; then
    echo "Output already exists for ${trait}, skipping"
    exit 0
fi

Rscript ${scripts_dir}/FUSION.assoc_test_trans.R \
    --chr $chr \
    --ref_ld_chr $ld_ref \
    --sumstats $gwas_sumstat_path \
    --weights $wgt \
    --weights_dir $wgtdir \
    --out ${out_dir}/${trait}.dat

echo "done with ${trait}"
