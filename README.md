# EGRET (Estimating Genome-wide Regulatory Effects on the Transcriptome)

EGRET is a multivariate linear model designed to identify genome-wide loci that are predictive of gene expression levels. EGRET integrates predictions from existing trans-eQTL mapping approaches (Matrix eQTL, GBAT, and trans-PCO) and determines the optimal weighted combination of regulatory variants that best explains gene expression. Finally, EGRET utilizes a genome-wide summary statistics-based TWAS to identify novel gene–disease associations.

EGRET produces per-gene prediction models (with cross-validated R²) and, optionally, gene–trait association results from a summary-statistics-based TWAS.

## Table of Contents

- [Pipeline overview](#pipeline-overview)
- [Quick start with example data](#quick-start-with-example-data)
- [Installation](#installation)
- [Running the full pipeline](#running-the-full-pipeline)
- [Setup](#setup)
- [Running Matrix eQTL, GBAT, and _trans_-PCO](#running-matrix-eqtl-gbat-and-trans-pco)
- [Training EGRET models](#training-egret-models)
- [Running EGRET-TWAS](#running-egret-twas)
- [License](#license)
- [Contact](#contact)

## Pipeline overview

EGRET runs in four stages:

1. **Preprocessing** — LD-prune genotypes, filter and regress expression, define cross-validation folds and cis-windows.
2. **Trans-eQTL mapping** — three approaches run in parallel: Matrix eQTL (univariate SNP–gene association), GBAT (polygenic trans-association using cis models), and trans-PCO (WGCNA module-based multivariate association).
3. **Model training** — integrate the trans-eQTL results with cis weights and train per-gene prediction models (lasso, elastic net, BLUP, xtune).
4. **TWAS** — run FUSION-style gene–trait association against GWAS summary statistics.

## Quick start with example data

> **Coming soon.** A small example dataset and an end-to-end tutorial walking through the four pipeline stages will be added in a follow-up release.

## Installation

```
git clone https://github.com/AmariutaLab/EGRET.git
cd EGRET
```

### Setting up the environment

Create the conda environment and install the post-install dependencies (`xtune` and `plink2R`, installed from GitHub):

```
conda env create -f environment.yml
conda activate egret_env
Rscript post_install.R
```

`environment.yml` and `post_install.R` are authoritative for all R/Python package requirements.

### External binaries

EGRET ships with the binaries it needs:

- `plink2` — bundled at the repo root.
- `gemma-0.98.5-linux-static-AMD64` — bundled under `scripts/`.

If you want to use a different PLINK build, see https://www.cog-genomics.org/plink/ and pass the path to your binary via the `plink_path` argument below.

## Running the full pipeline

The end-to-end driver is [`scripts/run_EGRET_workflow.sh`](scripts/run_EGRET_workflow.sh). Edit the parameter block at the top (tissue name, file paths, FDR, number of folds, model types, `output_dir`, etc.) and submit it via SLURM. It chains together every per-stage script documented below.

The per-stage instructions in the rest of this README are for running each stage independently.

## Setup

To setup expression and genotype files in the format that EGRET requires, we have created a wrapper script to carry out preprocessing of genotype and expression data. All outputs are stored in the `output_dir`.

| Argument                 | Description                                                                                                   |
| ------------------------ | ------------------------------------------------------------------------------------------------------------- |
| `genotypes_file_path`    | Path to PLINK `.bed/.bim/.fam` files containing genotype data.                                                |
| `expression_file_path`   | Gene expression matrix (rows = genes, columns = individuals). Must include a `gene_id` column in ENSG format. |
| `covariates_file_path`   | Covariate file (rows = individuals, columns = covariates to regress out).                                     |
| `tissue`                 | Name of tissue for analysis. Used for labeling and output organization.                                       |
| `individuals_file_path`  | File containing list of individual IDs (one per line, no header).                                             |
| `LD_prune_r2`            | LD pruning threshold (default: `0.9`).                                                                        |
| `plink_path`             | Path to PLINK2 executable.                                                                                    |
| `genotype_output_prefix` | Prefix for pruned genotype output files.                                                                      |
| `folds`                  | Number of cross-validation folds (e.g., 5).                                                                   |
| `gene_info_file_path`    | File containing gene metadata (gene ID, name, chromosome, start position).                                    |
| `FDR`                    | False discovery rate (FDR) that will be used to select significant variants.                                  |
| `output_dir`             | Directory to store all processed outputs.                                                                     |

```
./setup_genotype_and_expression.sh \
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
    $output_dir
```

## Running Matrix eQTL, GBAT, and _trans_-PCO

In EGRET, we implement Matrix eQTL, GBAT, and _trans_-PCO methods to identify _trans_-variants with potential to enhance gene expression models. Each one is run separately and the outputs of each is used in the final training of the genome-wide model.

### Running Matrix eQTL

We utilize the R package, Matrix eQTL, to conduct pairwise association between all genetic variants and gene expression. As input this script takes the tissue you wish to carry this analysis on. It is assumed that the setup scripts were run so it looks for their respective outputs.

```
./MatrixeQTL_scripts.sh \
    $tissue \
    $FDR \
    $output_dir
```

### Running GBAT

GBAT utilizes previously trained cis-expression models to find genome-wide regulators. This code is designed to take as input models fitted by FUSION.
* If FUSION models have not already been generated, model weights can be downloaded from http://gusevlab.org/projects/fusion/

```
./GBAT_scripts.sh \
    $FDR \
    $folds \
    $gene_info \
    $output_dir
```

### Running _trans_-PCO

transPCO finds the association between a module of genes and a variant.

```
./transPCO_scripts.sh \
    $tissue \
    $folds \
    $FDR \
    $num_PCs \
    $output_dir
```

## Training EGRET models

This script takes the input from the training of the models above and trains the genome-wide EGRET model and outputs the weights and model R². We recommend using at least Matrix eQTL and GBAT.

In addition this script will filter cross-mappable reads from alignment errors. As default we supply bed files containing cross-mappable beds for each gene in the hg38 genome build. If you want to filter cross-mappable regions in a different genome build, you can run the script `12.0_make_crossmappable_gene_beds.R` with gene–gene mappability scores provided by https://doi.org/10.12688/f1000research.17145.2.

```
./train_EGRET_models.sh \
    $tissue \
    $folds \
    $models \
    $output_dir
```

## Running EGRET-TWAS

To conduct gene–trait association, we have provided the following script which takes as input GWAS summary statistics and the models trained by EGRET to produce gene–trait associations.

```
./run_TWAS_scripts.sh \
    $tissue \
    $trait \
    $gwas_sumstat_dir \
    $output_dir
```

## License

EGRET is released under the MIT License. See [`LICENSE`](LICENSE).

## Contact

Questions, bug reports, and feature requests: please open an issue at [AmariutaLab/EGRET](https://github.com/AmariutaLab/EGRET/issues).
