#!/bin/bash
#SBATCH --job-name="Running Fusion"
#SBATCH --output="batch_submissions/cis_fusion.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2G
#SBATCH --account=csd832
#SBATCH --export=ALL
#SBATCH -t 24:00:00


tissue=$1
output_dir=$2
chunk_size=$3
end_index=$4
plink_path=$5
individuals_file_path=$6
covariates_file_path=$7
plink1_path=$8
gcta_path=$9
gemma_path=${10}
scripts_dir=${11}
cleanup_cis_working=${12:-"FALSE"}

if [ -n "$tissue" ]; then
    echo $tissue
else
    echo "no tissue defined"
    exit 1
fi

#Define output directions
wd=${output_dir}/working/$tissue/cis
mkdir -p $wd
tmpdir=${output_dir}/tmp/$tissue/cis
mkdir -p $tmpdir
weights=${output_dir}/FUSION/$tissue/cis
mkdir -p $weights

for gene in $(ls ${output_dir}/plink_results/${tissue}/cis | grep .bed | head -n $end_index | tail -n $chunk_size)
do
	gene=$(basename $gene .bed)
	gefile=${output_dir}/expression_files/${tissue}_expression.txt.gz
	individuals=${individuals_file_path}
	individuals_plink=${wd}/${gene}_keep.txt
	awk '{print "0\t"$1}' $individuals > $individuals_plink
	covar=${covariates_file_path}

	alldonors=$(zcat $gefile | head -n 1)
	colind=$(echo $alldonors | sed 's| |\n|g' | nl | grep -f $individuals | awk 'BEGIN {ORS=","} {print $1}' | sed 's/,$//')

	${plink_path} --bfile ${output_dir}/plink_results/${tissue}/cis/$gene --make-bed --keep $individuals_plink --indiv-sort file $individuals_plink --out $wd/${gene}
	rm $wd/${gene}.log

	rowid=$(zcat $gefile | nl | grep ${gene} | awk '{print $1}')
	ge_EURdonors=$(zcat $gefile | head -n $rowid | tail -n 1 | cut -f $colind)
	paste --delimiters='\t' <(cut -f1-5 $wd/${gene}.fam) <(echo $ge_EURdonors | sed 's/ /\n/g') > $wd/${gene}.mod.fam
	mv $wd/${gene}.mod.fam $wd/${gene}.fam

	TMP=$tmpdir/${gene}
	OUT=$weights/${gene}

	Rscript ${scripts_dir}/FUSION.compute_weights.R --bfile $wd/${gene} --crossval 5 --models lasso,blup,enet,top1 --hsq_p 1 --tmp $TMP --out $OUT --covar $covar --PATH_gcta ${gcta_path} --PATH_plink ${plink1_path} --PATH_gemma ${gemma_path} --noclean FALSE

	# Hook A: drop per-gene FUSION working/tmp files once the cis weight is written
	if [[ "$cleanup_cis_working" == "TRUE" && -f "$OUT.wgt.RDat" ]]; then
		rm -rf $wd/${gene}.* $tmpdir/${gene}*
	fi
done
