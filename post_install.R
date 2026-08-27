# post_install.R — run after creating the conda environment
# Usage: conda activate egret_env && Rscript post_install.R
#
# Installs the R packages that are not distributed through the channels environment.yml
# configures. Run from the repository root so the local plink2R source is found.

library(devtools)

# xtune — regularized regression with external information
devtools::install_github("JingxuanH/xtune", build_vignettes = FALSE)

# ACAT — aggregated Cauchy association test. Optional: no script in the pipeline currently
# consumes it (the only ACAT caller, ModifiedPCOMerged_acat_version.R, is not part of the
# tracked workflow). Kept because swapping 9.3's PCO merge to the ACAT variant is a
# one-line change and the install is cheap.
devtools::install_github("yaowuliu/ACAT")

# plink2R — read PLINK binary files. Required by 1.3_setup_folds.R, 4.1_run_GBAT.R,
# 5.2_initialize_chr_genotypes.R, 5.4_MatrixeQTL_for_transPCO.R, and
# 6.4_train_EGRET_model_weights.R. plink2R-master/ is gitignored, so a fresh clone does not
# have it: fall back to GitHub, and fail loudly if that is unavailable too. This must never
# be a message()-and-continue -- silently skipping the install leaves an environment that
# looks complete and then dies five stages into a submitted run.
if (dir.exists("plink2R-master/plink2R/")) {
  install.packages("plink2R-master/plink2R/", repos = NULL)
} else {
  message("plink2R-master/ not found locally; installing from GitHub (gabraham/plink2R)")
  try(devtools::install_github("gabraham/plink2R", subdir = "plink2R"), silent = FALSE)
}

if (!requireNamespace("plink2R", quietly = TRUE)) {
  stop("plink2R is not installed. Download https://github.com/gabraham/plink2R into ",
       "plink2R-master/ at the repository root and re-run this script.")
}

message("post_install complete: ",
        paste(Filter(function(p) requireNamespace(p, quietly = TRUE),
                     c("xtune", "ACAT", "plink2R")), collapse = ", "))
