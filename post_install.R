# post_install.R — run after creating the conda environment
# Usage: conda activate egret_env && Rscript post_install.R

library(devtools)

# xtune — regularized regression with external information
devtools::install_github("JingxuanH/xtune", build_vignettes = FALSE)

# ACAT — aggregated Cauchy association test
devtools::install_github("yaowuliu/ACAT")

# plink2R — read PLINK binary files
if (dir.exists("plink2R-master/plink2R/")) {
  install.packages("plink2R-master/plink2R/", repos = NULL)
} else {
  message("plink2R-master/ not found. Download from https://github.com/gabraham/plink2R")
}
