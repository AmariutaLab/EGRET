#!/bin/bash
# --- config_validation.sh --------------------------------------------------
# Reduced-scale twin of config.sh for the end-to-end validation run: fewer genes, fewer
# SNPs, two GWAS traits. Everything else matches production.
#
#   in : nothing -- sourced, never executed
#   out: shell variables
#   run: sbatch run_EGRET_workflow.sh --config config_validation.sh
# --------------------------------------------------------------------------
# Reduced configuration for an end-to-end validation run.
#
#   ./validate_egret.sh prepare
#   sbatch run_EGRET_workflow.sh --config config_validation.sh
#   ./validate_egret.sh check
#
# Sized to exercise every stage and every cleanup hook without a full-scale run: ~1.8k
# genes (enough for WGCNA to build >=10-gene modules), ~80k SNPs spread over all 22
# chromosomes (so no per-chromosome loop hits an empty one), and 2 GWAS traits.
#
# Every cleanup_* flag below is TRUE so the run bounds its own disk use. Report the
# footprint with `./validate_egret.sh usage`; remove it with `rm -rf scripts/validation`.

# Input data -- produced by ./validate_egret.sh prepare. Anchors resolve as in config.sh.
_scripts_dir="${EGRET_SCRIPTS_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
_repo_dir="$(cd "${_scripts_dir}/.." && pwd)"
_val_dir="${_scripts_dir}/validation"
genotype_prefix="${_val_dir}/inputs/validation_genotypes"
expression="${_val_dir}/inputs/validation_expression.bed.gz"
individuals="/expanse/lustre/projects/ddp412/kbrunton/TWAS_across_tissues/individuals_per_tissue/Whole_Blood_individuals.txt"
covariates="/expanse/lustre/projects/ddp412/tamariutabartell/o2_files/gtex_covars/Covar_all_Whole_Blood.txt"
gene_info="/expanse/lustre/projects/ddp412/kbrunton/data/GTEx_V8.txt.gz"
# GRN eQTL sumstats. Must be the *_with_rsid.txt variant -- 9.3 joins on an `rsid` column.
eqtl_file="/expanse/lustre/projects/ddp412/kbrunton/EGRET/grn_analysis_files/GTEx_Analysis_v8_eQTL/Whole_Blood_with_rsid.txt"
# A target gene enters the GRN graph when its trans partition clears this p-value.
grn_pval_threshold="0.01"

# External programs
plink_path="${_repo_dir}/plink2"
plink1_path="/expanse/lustre/projects/ddp412/kbrunton/plink"
gcta_path="/expanse/lustre/projects/ddp412/kbrunton/fusion_twas-master/gcta_nr_robust"
gemma_path="${_scripts_dir}/gemma-0.98.5-linux-static-AMD64"

# Tissue, models, and output layout
tissue="Whole_Blood"
output_dir="${_val_dir}/output"
genotype_output_prefix="GTEX_v8_genotypes_pruned"
egret_subdir="EGRET"
# folds=5 deliberately, NOT 2. 1.3_setup_folds.R trains each fold on (folds-1)/folds of the
# cohort, while 2.2_compute_cis_weights_chunk.sh hands cis-FUSION a hardcoded --crossval 5.
# Only folds=5 gives the two arms equal training size; anything lower shrinks EGRET's
# training set while cis stays at 80%, depressing EGRET's R2 across every model class and
# biasing any EGRET-vs-cis comparison against EGRET. Lower it only for a plumbing smoke
# test, never for a comparison.
folds="5"
ld_prune_r2="0.9"
fdr="0.1"
chunk_size="300"
bed_models="Cis,MatrixeQTL,GBAT,transPCO"
fusion_models="xtune,lasso,enet,blup"
matrixeqtl_bed_dir="results_FDR_0.1"
gbat_bed_dir="results_FDR_0.1"
transpco_bed_dir="bed_files_FDR_0.1"
exclude_crossmap="FALSE"
crossmap_dir="background_mismatches"
crossmap_file="cross_mappability.tsv.gz"
covariate_columns_for_coexpression="PC1 PC2 PC3 PC4 PC5 InferredCov1 InferredCov2 InferredCov3 InferredCov4 InferredCov5 InferredCov6 InferredCov7 InferredCov8 InferredCov9 InferredCov10 pcr platform sex"

# Stage selection. Skipped stages still have their outputs asserted, so a typo in this
# list cannot silently skip real work.
stages="${EGRET_STAGES:-setup cis matrixeqtl gbat transpco models analysis twas twas_cis}"
cleanup_cis_plink="TRUE"
cleanup_cis_working="TRUE"
cleanup_matrixeqtl_raw="TRUE"
cleanup_transpco_raw="TRUE"
cleanup_egret_plink="TRUE"
cleanup_egret_working="TRUE"
cleanup_egret_z_matrices="TRUE"
cleanup_bed_lists="TRUE"
cleanup_fold_genotypes="TRUE"
cleanup_fold_expression="TRUE"

# TWAS inputs
gwas_sumstat_dir="/expanse/lustre/projects/ddp412/kbrunton/TCSC/sumstats"
ld_ref="${_scripts_dir}/LDREF/1000G.EUR.merged"
gwas_sumstats=(
  UKB_460K.blood_WHITE_COUNT.sumstats.gz
  PASS_Rheumatoid_Arthritis.sumstats.gz
)

# Runtime and scheduler
conda_env="egret_env"
slurm_account="csd832"
slurm_partition="shared"
slurm_log_dir="${_val_dir}/slurm_logs"
