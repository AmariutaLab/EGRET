#!/bin/bash
#SBATCH --job-name="Fusion"
#SBATCH --output="batch_submissions/gtex_xtune.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH --account=csd832
#SBATCH --export=ALL
#SBATCH -t 24:00:00


tissue=$1
chunk_size=$2
end_index=$3
base_dir=${4:-"."}
model_config=${5:-"EGRET"}
folds=${6:-5}
plink_path=${7:-"plink2"}
gemma_path=${8:-"../../../gemma-0.98.5-linux-static-AMD64"}
expression_file_path=${9:-"expression_files/${tissue}_expression.txt.gz"}
fusion_models=${10:-"xtune,lasso,enet,blup"}
gene_info_file_path=${11:-"../../data/GTEx_V8.txt.gz"}

if [ -n "$tissue" ]; then
    echo $tissue
else
    echo "no tissue defined"
    return
fi

home_dir=/expanse/lustre/projects/ddp412/kbrunton

#Define output directions
weights="${base_dir}/xtune_fusion_models/${tissue}/${model_config}_no_cis"
mkdir -p "$weights"
wd="${base_dir}/working/${tissue}/${model_config}_no_cis"
mkdir -p "$wd"
output="${base_dir}/xtune_fusion_results/${tissue}/${model_config}_no_cis"
mkdir -p "$output"

z_matrix_dir="${base_dir}/z_matrices/${tissue}/${model_config}_no_cis"

plink_dir="${base_dir}/plink_results/${tissue}/${model_config}"

for gene in $(ls ${plink_dir}/fold_0/ | grep .bed | head -n $end_index | tail -n $chunk_size)
#for gene in $(find transPCO/Whole_Blood/fold_1/bed_files/ transPCO/Whole_Blood/fold_2/bed_files/ transPCO/Whole_Blood/fold_3/bed_files/ transPCO/Whole_Blood/fold_4/bed_files/ transPCO/Whole_Blood/fold_5/bed_files/ -type f -exec basename {} \; | sort | uniq)
do
gene=${gene::-4}

#Extracting gene expression information
gefile="$expression_file_path"
individuals="${base_dir}/fold_0_info/${tissue}/train_individuals.txt"

alldonors=$(zcat $gefile | head -n 1)
colind=$(echo $alldonors | sed 's| |\n|g' | nl | grep -f $individuals | awk 'BEGIN {ORS=","} {print $1}') #donor list always the same
colind2=${colind%,}

for i in $(seq 0 $folds)
do
	mkdir -p $wd/fold_$i
	${plink_path} --bfile $plink_dir/fold_$i/$gene --make-bed --keep $individuals --indiv-sort f $individuals --out $wd/fold_$i/${gene}
	rm $wd/fold_$i/${gene}.log

	rowid=$(zcat $gefile | awk 'NR > 1 {print $1}' | nl | grep ${gene} | awk '{print $1 + 1}')
	ge_EURdonors=$(zcat $gefile | head -n $rowid | tail -n 1 | cut -f $colind2)
	paste --delimiters='\t' <(cut -f1-5 $wd/fold_$i/${gene}.fam) <(echo $ge_EURdonors | sed 's/ /\n/g') > $wd/fold_$i/${gene}.mod.fam
	mv $wd/fold_$i/${gene}.mod.fam $wd/fold_$i/${gene}.fam

done
echo "running xtune fusion"
Rscript xtune_fusion_cis_trans.R \
    --gene $gene \
    --working_dir "$wd" \
    --models "$fusion_models" \
    --output_dir "$output" \
    --weights_dir "$weights" \
    --tissue "$tissue" \
    --PATH_plink "$plink_path" \
    --PATH_gemma "$gemma_path" \
    --covar "${base_dir}/covariate_files/${tissue}_covariates.txt.gz" \
    --z_matrix_dir "$z_matrix_dir" \
    --gene_info_file "$gene_info_file_path" \
    --folds "$folds"
done
