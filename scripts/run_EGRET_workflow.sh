#!/bin/bash
#SBATCH --job-name="Run_Script"
#SBATCH --output="TEST_SETUP.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=30G
#SBATCH --account=csd832
#SBATCH --export=ALL
#SBATCH -t 08:00:00

source ~/.bashrc
conda activate transTWAS_env

home_dir="/expanse/lustre/projects/ddp412/kbrunton/EGRET"
scripts_dir="${home_dir}/scripts"
tissue="Whole_Blood"

genotypes_file_path="/expanse/lustre/projects/ddp412/tamariutabartell/gtex/GTEx_v8_genotype_EUR_HM3_exclude_dups.allchr.reorder"
expression_file_path="/expanse/lustre/projects/ddp412/tamariutabartell/gtex/${tissue}.v8.normalized_expression.bed.gz"
individuals_file_path="/expanse/lustre/projects/ddp412/kbrunton/TWAS_across_tissues/individuals_per_tissue/${tissue}_individuals.txt"
covariates_file_path="/expanse/lustre/projects/ddp412/tamariutabartell/o2_files/gtex_covars/Covar_all_${tissue}.txt"
LD_prune_r2="0.9"
plink_path="${home_dir}/plink2"
plink1_path="/expanse/lustre/projects/ddp412/kbrunton/plink"
gcta_path="/expanse/lustre/projects/ddp412/kbrunton/fusion_twas-master/gcta_nr_robust"
genotype_output_prefix="GTEX_v8_genotypes_pruned"
folds="5"
gene_info_file_path="/expanse/lustre/projects/ddp412/kbrunton/data/GTEx_V8.txt.gz"
FDR="0.1"
covariate_columns_for_coexpression="PC1 PC2 PC3 PC4 PC5 InferredCov1 InferredCov2 InferredCov3 InferredCov4 InferredCov5 InferredCov6 InferredCov7 InferredCov8 InferredCov9 InferredCov10 pcr     platform sex"
models="lasso,enet,blup,xtune"
output_dir="test_output/"
MatrixeQTL_bed_dir="results_FDR_${FDR}"
GBAT_bed_dir="results_FDR_${FDR}"
transPCO_bed_dir="bed_files_FDR_${FDR}"
egret_output_subdir="EGRET"
fusion_models="xtune,lasso,enet,blup"
gemma_path="${scripts_dir}/gemma-0.98.5-linux-static-AMD64"
chunk_size=500
gwas_sumstat_dir="/expanse/lustre/projects/ddp412/kbrunton/TCSC/sumstats"
ld_ref="/expanse/lustre/projects/ddp412/kbrunton/fusion_twas-master/LDREF/1000G.EUR.merged"
crossmap_file="cross_mappability.tsv.gz"  # only used if regenerating cross-mappable beds (e.g. non-hg38)

gwas_sumstats=("PASS_IBD_deLange2017.sumstats.gz" "PASS_Rheumatoid_Arthritis.sumstats.gz" "UKB_460K.blood_PLATELET_COUNT.sumstats.gz" "UKB_460K.blood_RBC_DISTRIB_WIDTH.sumstats.gz" "UKB_460K.blood_RED_COUNT.sumstats.gz" "UKB_460K.blood_WHITE_COUNT.sumstats.gz" "UKB_460K.disease_ALLERGY_ECZEMA_DIAGNOSED.sumstats.gz" "UKB_460K.disease_HYPOTHYROIDISM_SELF_REP.sumstats.gz" "PASS_AtrialFibrillation_Nielsen2018.sumstats.gz" "PASS_HDL.sumstats.gz" "PASS_IschemicStroke_Malik2018.sumstats.gz" "PASS_LDL.sumstats.gz" "UKB_460K.biochemistry_Cholesterol.sumstats.gz" "UKB_460K.body_BMIz.sumstats.gz" "UKB_460K.bp_DIASTOLICadjMEDz.sumstats.gz" "PASS_ADHD_Demontis2018.sumstats.gz" "PASS_Alzheimers_deRojas2021.sumstats.gz" "PASS_Autism.sumstats.gz" "PASS_Schizophrenia_Pardinas2018.sumstats.gz" "PASS_AnorexiaNervosa_Watson2019.sumstats.gz" "UKB_460K.cov_EDU_YEARS.sumstats.gz" "UKB_460K.lung_FVCzSMOKE.sumstats.gz" "PASS_Intelligence_SavageJansen2018.sumstats.gz" "UKB_460K.mental_NEUROTICISM.sumstats.gz" "UKB_460K.body_HEIGHTz.sumstats.gz" "UKB_460K.other_MORNINGPERSON.sumstats.gz" "UKB_460K.repro_MENARCHE_AGE.sumstats.gz" "UKB_460K.biochemistry_TotalProtein.sumstats.gz" "UKB_460K.body_WHRadjBMIz.sumstats.gz" "UKB_460K.lung_FEV1FVCzSMOKE.sumstats.gz" "PASS_ReactionTime_Davies2018.sumstats.gz" "PASS_Myopia_Hysi2020.sumstats.gz" "UKB_460K.biochemistry_Creatinine.sumstats.gz" "UKB_460K.bmd_HEEL_TSCOREz.sumstats.gz" "PASS_SleepDuration_Dashti2019.sumstats.gz" "UKB_460K.biochemistry_IGF1.sumstats.gz" "PASS_GeneralRiskTolerance_KarlssonLinner2019.sumstats.gz" "PASS_Eosino_Vuckovic2020.sumstats.gz" "UKB_460K.biochemistry_AspartateAminotransferase.sumstats.gz" "PASS_Insomnia_Jansen2019.sumstats.gz" "PASS_Height1.sumstats.gz" "PASS_LymP_Vuckovic2020.sumstats.gz" "UKB_460K.repro_NumberChildrenEverBorn_Pooled.sumstats.gz" "UKB_460K.biochemistry_Testosterone_Male.sumstats.gz" "PASS_BMI1.sumstats.gz" "PASS_AgeFirstBirth.sumstats.gz" "PASS_Glaucoma_Craig2020.sumstats.gz" "PASS_RTC_Vuckovic2020.sumstats.gz" "UKB_460K.body_BALDING1.sumstats.gz" "UKB_460K.biochemistry_AlkalinePhosphatase.sumstats.gz" "PASS_DrinksPerWeek_Liu2019.sumstats.gz" "UKB_460K.biochemistry_Phosphate.sumstats.gz" "PASS_MedicationUse_Wu2019.sumstats.gz" "UKB_460K.biochemistry_VitaminD.sumstats.gz" "PASS_MDD_Wray2018.sumstats.gz" "PASS_RD_Zhao2021.sumstats.gz" "PASS_MonoP_Vuckovic2020.sumstats.gz" "UKB_460K.biochemistry_TotalBilirubin.sumstats.gz" "UKB_460K.repro_MENOPAUSE_AGE.sumstats.gz" "PASS_SA_Grasby2020.sumstats.gz" "PASS_HipOA_Tachmazidou2019.sumstats.gz" "UKB_460K.pigment_SUNBURN.sumstats.gz" "PASS_NumberChildrenEverBorn.sumstats.gz" "PASS_BipolarDisorder_Ruderfer2018.sumstats.gz" "PASS_Years_of_Education1.sumstats.gz" "PASS_PancreasVol_Liu2021.sumstats.gz" "PASS_CaudateVol_Satizabal2019.sumstats.gz" "PASS_BrainstemVol_Satizabal2019.sumstats.gz" "PASS_TH_Grasby2020.sumstats.gz" "PASS_BasoP_Vuckovic2020.sumstats.gz" "PASS_Anorexia.sumstats.gz" "PASS_Ever_Smoked.sumstats.gz" "PASS_MO_Zhao2021.sumstats.gz" "PASS_ProstateCancer.sumstats.gz" "UKB_460K.cancer_PROSTATE.sumstats.gz" "PASS_AccumbensVol_Satizabal2019.sumstats.gz" "UKB_460K.cancer_BREAST.sumstats.gz" "PASS_Type_2_Diabetes.sumstats.gz" "PASS_BipolarDisorder_Ruderfer2018.sumstats.gz")

gwas_sumstats=("PASS_IBD_deLange2017.sumstats.gz")

gwas_sumstats_csv=$(IFS=','; echo "${gwas_sumstats[*]}")

if false ; then
${scripts_dir}/setup_genotype_and_expression.sh \
    $genotypes_file_path \
    $expression_file_path \
    $covariates_file_path \
    $tissue \
    $individuals_file_path \
    $LD_prune_r2 \
    $plink_path \
    $genotype_output_prefix \
    $folds \
    $gene_info_file_path \
    $output_dir \
    $scripts_dir

${scripts_dir}/run_cis_models.sh \
    $tissue \
    $output_dir \
    $plink_path \
    $gene_info_file_path \
    $genotype_output_prefix \
    $individuals_file_path \
    $covariates_file_path \
    $plink1_path \
    $gcta_path \
    $gemma_path \
    $chunk_size \
    $scripts_dir

${scripts_dir}/MatrixeQTL_scripts.sh \
    $tissue \
    $FDR \
    $folds \
    $gene_info_file_path \
    $genotype_output_prefix \
    $output_dir \
    $scripts_dir

cis_model_dir="${home_dir}/FUSION/${tissue}/cis/"
${scripts_dir}/GBAT_scripts.sh \
    $tissue \
    $FDR \
    $folds \
    $gene_info_file_path \
    $cis_model_dir \
    $plink_path \
    $genotypes_file_path \
    $output_dir \
    $scripts_dir


${scripts_dir}/transPCO_scripts.sh \
    $tissue \
    $folds \
    $FDR \
    "$covariate_columns_for_coexpression" \
    $gene_info_file_path \
    $plink_path \
    $output_dir \
    $genotype_output_prefix \
    $scripts_dir

# Optional: regenerate cross-mappable bed files. Only needed if working with a non-hg38
# genome build; hg38 cross-mappable beds are bundled in the repo.
# Rscript ${scripts_dir}/12.0_make_crossmappable_gene_beds.R \
#     --gene_info $gene_info_file_path \
#     --crossmap_file $crossmap_file \
#     --output_dir $output_dir

${scripts_dir}/train_EGRET_models.sh \
    $tissue \
    $folds \
    $models \
    $plink_path \
    $output_dir \
    $MatrixeQTL_bed_dir \
    $GBAT_bed_dir \
    $transPCO_bed_dir \
    TRUE \
    background_mismatches \
    $genotype_output_prefix \
    $egret_output_subdir \
    $fusion_models \
    $gemma_path \
    $chunk_size \
    $gene_info_file_path \
    $scripts_dir


sbatch ${scripts_dir}/run_EGRET_model_analysis.sh \
    $tissue \
    $egret_output_subdir \
    $output_dir \
    $scripts_dir

fi

${scripts_dir}/run_TWAS_scripts.sh \
    $tissue \
    $output_dir \
    $gene_info_file_path \
    $gwas_sumstat_dir \
    $ld_ref \
    $egret_output_subdir \
    $gwas_sumstats_csv \
    $scripts_dir

# Cis-only TWAS (uses cis weights from FUSION) — useful as a baseline against the full EGRET TWAS.
${scripts_dir}/run_TWAS_cis_scripts.sh \
    $tissue \
    $output_dir \
    $gene_info_file_path \
    $gwas_sumstat_dir \
    $ld_ref \
    $egret_output_subdir \
    $gwas_sumstats_csv \
    $scripts_dir
