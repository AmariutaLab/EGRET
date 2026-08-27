#!/bin/bash
# --- config.sh -------------------------------------------------------------
# Every tunable parameter for a production run: input paths, binaries, stage selection,
# FDR and fold counts, and the cleanup switches.
# Sourced by run_EGRET_workflow.sh; edit this, not the scripts.
#
#   in : nothing -- sourced, never executed
#   out: shell variables
#   run: sourced by run_EGRET_workflow.sh (override with --config or EGRET_CONFIG)
# --------------------------------------------------------------------------

# Repo anchors, so a clone at any path works unedited. EGRET_SCRIPTS_DIR is tried first
# (run_EGRET_workflow.sh exports it), which keeps a config copied outside the repo anchored
# to the real scripts/ directory; this file's own location is the fallback.
_scripts_dir="${EGRET_SCRIPTS_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
_repo_dir="$(cd "${_scripts_dir}/.." && pwd)"

# Input data
genotype_prefix="/expanse/lustre/projects/ddp412/tamariutabartell/gtex/GTEx_v8_genotype_EUR_HM3_exclude_dups.allchr.reorder"
expression="/expanse/lustre/projects/ddp412/tamariutabartell/gtex/Whole_Blood.v8.normalized_expression.bed.gz"
individuals="/expanse/lustre/projects/ddp412/kbrunton/TWAS_across_tissues/individuals_per_tissue/Whole_Blood_individuals.txt"
covariates="/expanse/lustre/projects/ddp412/tamariutabartell/o2_files/gtex_covars/Covar_all_Whole_Blood.txt"
gene_info="/expanse/lustre/projects/ddp412/kbrunton/data/GTEx_V8.txt.gz"
# GRN eQTL sumstats. Must be the *_with_rsid.txt variant, not the stock GTEx file: 9.3
# joins on an `rsid` column that only convert_variant_id_to_rsid.py adds.
eqtl_file="/expanse/lustre/projects/ddp412/kbrunton/EGRET/grn_analysis_files/GTEx_Analysis_v8_eQTL/Whole_Blood_with_rsid.txt"
# A target gene enters the GRN graph when its trans partition clears this p-value.
# Distinct from `fdr`, which governs trans-eQTL SNP selection upstream.
grn_pval_threshold="0.01"

# External programs. plink2 is gitignored, so a fresh clone must supply it; the path is
# repo-relative, so dropping the binary in place is all that is needed.
plink_path="${_repo_dir}/plink2"
plink1_path="/expanse/lustre/projects/ddp412/kbrunton/plink"
gcta_path="/expanse/lustre/projects/ddp412/kbrunton/fusion_twas-master/gcta_nr_robust"
gemma_path="${_scripts_dir}/gemma-0.98.5-linux-static-AMD64"   # bundled in the repo

# Tissue, models, and output layout
tissue="Whole_Blood"
# Must be absolute: a relative path resolves against the submitting CWD, so the same
# config would write to a different tree depending on where sbatch was run from.
output_dir="${_scripts_dir}/test_output"
genotype_output_prefix="GTEX_v8_genotypes_pruned"
egret_subdir="EGRET"
folds="5"
ld_prune_r2="0.9"
fdr="0.1"
chunk_size="500"
bed_models="Cis,MatrixeQTL,GBAT,transPCO"
fusion_models="xtune,lasso,enet,blup"
matrixeqtl_bed_dir="results_FDR_0.1"
gbat_bed_dir="results_FDR_0.1"
transpco_bed_dir="bed_files_FDR_0.1"
exclude_crossmap="TRUE"
crossmap_dir="background_mismatches"
# Consumed by the manual cross-mappability prep step, which must be run once per output_dir
# before the models stage whenever exclude_crossmap=TRUE (the BED files are not bundled):
#   Rscript 0.1_make_crossmappable_gene_beds.R --gene_info "$gene_info" \
#       --crossmap_file "$crossmap_file" --output_dir "$output_dir" \
#       --crossmap_dir "$crossmap_dir"
# Skipping this prep step is a hard error in 6.2, not a silent fall-through.
crossmap_file="cross_mappability.tsv.gz"
covariate_columns_for_coexpression="PC1 PC2 PC3 PC4 PC5 InferredCov1 InferredCov2 InferredCov3 InferredCov4 InferredCov5 InferredCov6 InferredCov7 InferredCov8 InferredCov9 InferredCov10 pcr platform sex"

# Stage selection. EGRET_STAGES lets a resumed run re-enter partway through. Skipped
# stages still have their outputs asserted, so a typo cannot silently skip real work.
# `grn` is opt-in: name it explicitly to enable it.
stages="${EGRET_STAGES:-setup cis matrixeqtl gbat transpco models analysis twas twas_cis}"
cleanup_cis_plink="TRUE"
cleanup_cis_working="TRUE"
cleanup_matrixeqtl_raw="TRUE"
cleanup_transpco_raw="TRUE"
cleanup_egret_plink="TRUE"
cleanup_egret_working="TRUE"
# Inodes, not bytes, are the binding constraint here: one tissue can exceed a 2M inode
# quota, most of it in working/ and plink_results/. Set any FALSE to keep intermediates.
cleanup_egret_z_matrices="TRUE"   # z_matrices/{t}/{sub}/fold_*/{gene}_z_matrix.txt (~117k files)
cleanup_bed_lists="TRUE"          # bed_files/{t}/{sub}/fold_*/{gene}       (~118k files)
cleanup_fold_genotypes="TRUE"     # fold_*_info/{t}/train_genotypes.txt
cleanup_fold_expression="TRUE"    # fold_*_info/{t}/train_expression.txt

# TWAS inputs
gwas_sumstat_dir="/expanse/lustre/projects/ddp412/kbrunton/TCSC/sumstats"
# The bundled per-chromosome LDREF ships with the repo; the merged fileset FUSION needs is
# built once with `bash LDREF/merge_files.sh <plink1.9>`.
ld_ref="${_scripts_dir}/LDREF/1000G.EUR.merged"
gwas_sumstats=(
  PASS_IBD_deLange2017.sumstats.gz PASS_Rheumatoid_Arthritis.sumstats.gz
  UKB_460K.blood_PLATELET_COUNT.sumstats.gz UKB_460K.blood_RBC_DISTRIB_WIDTH.sumstats.gz
  UKB_460K.blood_RED_COUNT.sumstats.gz UKB_460K.blood_WHITE_COUNT.sumstats.gz
  UKB_460K.disease_ALLERGY_ECZEMA_DIAGNOSED.sumstats.gz UKB_460K.disease_HYPOTHYROIDISM_SELF_REP.sumstats.gz
  PASS_AtrialFibrillation_Nielsen2018.sumstats.gz PASS_HDL.sumstats.gz PASS_IschemicStroke_Malik2018.sumstats.gz
  PASS_LDL.sumstats.gz UKB_460K.biochemistry_Cholesterol.sumstats.gz UKB_460K.body_BMIz.sumstats.gz
  UKB_460K.bp_DIASTOLICadjMEDz.sumstats.gz PASS_ADHD_Demontis2018.sumstats.gz
  PASS_Alzheimers_deRojas2021.sumstats.gz PASS_Autism.sumstats.gz PASS_Schizophrenia_Pardinas2018.sumstats.gz
  PASS_AnorexiaNervosa_Watson2019.sumstats.gz UKB_460K.cov_EDU_YEARS.sumstats.gz
  UKB_460K.lung_FVCzSMOKE.sumstats.gz PASS_Intelligence_SavageJansen2018.sumstats.gz
  UKB_460K.mental_NEUROTICISM.sumstats.gz UKB_460K.body_HEIGHTz.sumstats.gz
  UKB_460K.other_MORNINGPERSON.sumstats.gz UKB_460K.repro_MENARCHE_AGE.sumstats.gz
  UKB_460K.biochemistry_TotalProtein.sumstats.gz UKB_460K.body_WHRadjBMIz.sumstats.gz
  UKB_460K.lung_FEV1FVCzSMOKE.sumstats.gz PASS_ReactionTime_Davies2018.sumstats.gz
  PASS_Myopia_Hysi2020.sumstats.gz UKB_460K.biochemistry_Creatinine.sumstats.gz
  UKB_460K.bmd_HEEL_TSCOREz.sumstats.gz PASS_SleepDuration_Dashti2019.sumstats.gz
  UKB_460K.biochemistry_IGF1.sumstats.gz PASS_GeneralRiskTolerance_KarlssonLinner2019.sumstats.gz
  PASS_Eosino_Vuckovic2020.sumstats.gz UKB_460K.biochemistry_AspartateAminotransferase.sumstats.gz
  PASS_Insomnia_Jansen2019.sumstats.gz PASS_Height1.sumstats.gz PASS_LymP_Vuckovic2020.sumstats.gz
  UKB_460K.repro_NumberChildrenEverBorn_Pooled.sumstats.gz UKB_460K.biochemistry_Testosterone_Male.sumstats.gz
  PASS_BMI1.sumstats.gz PASS_AgeFirstBirth.sumstats.gz PASS_Glaucoma_Craig2020.sumstats.gz
  PASS_RTC_Vuckovic2020.sumstats.gz UKB_460K.body_BALDING1.sumstats.gz
  UKB_460K.biochemistry_AlkalinePhosphatase.sumstats.gz PASS_DrinksPerWeek_Liu2019.sumstats.gz
  UKB_460K.biochemistry_Phosphate.sumstats.gz PASS_MedicationUse_Wu2019.sumstats.gz
  UKB_460K.biochemistry_VitaminD.sumstats.gz PASS_MDD_Wray2018.sumstats.gz PASS_RD_Zhao2021.sumstats.gz
  PASS_MonoP_Vuckovic2020.sumstats.gz UKB_460K.biochemistry_TotalBilirubin.sumstats.gz
  UKB_460K.repro_MENOPAUSE_AGE.sumstats.gz PASS_SA_Grasby2020.sumstats.gz
  PASS_HipOA_Tachmazidou2019.sumstats.gz UKB_460K.pigment_SUNBURN.sumstats.gz
  PASS_NumberChildrenEverBorn.sumstats.gz PASS_BipolarDisorder_Ruderfer2018.sumstats.gz
  PASS_Years_of_Education1.sumstats.gz PASS_PancreasVol_Liu2021.sumstats.gz
  PASS_CaudateVol_Satizabal2019.sumstats.gz PASS_BrainstemVol_Satizabal2019.sumstats.gz
  PASS_TH_Grasby2020.sumstats.gz PASS_BasoP_Vuckovic2020.sumstats.gz PASS_Anorexia.sumstats.gz
  PASS_Ever_Smoked.sumstats.gz PASS_MO_Zhao2021.sumstats.gz PASS_ProstateCancer.sumstats.gz
  UKB_460K.cancer_PROSTATE.sumstats.gz PASS_AccumbensVol_Satizabal2019.sumstats.gz
  UKB_460K.cancer_BREAST.sumstats.gz PASS_Type_2_Diabetes.sumstats.gz
)

# Runtime and scheduler
conda_env="transTWAS_env"
slurm_account="csd832"
slurm_partition="shared"
slurm_log_dir="${SCRATCH:-/expanse/lustre/scratch/kbrunton/temp_project}/EGRET_slurm_logs"
