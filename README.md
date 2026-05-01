# EGRET (Estimating Genome-wide Regulatory Effects on the Transcriptome)

EGRET is a multivariate linear model designed to identify genome-wide loci that are predictive of gene expression levels. EGRET integrates predictions from existing trans-eQTL mapping approaches (Matrix eQTL, GBAT, and trans-PCO) and determines the optimal weighted combination of regulatory variants that best explains gene expression. Finally, EGRET utilizes a genome-wide summary statistics-based TWAS to identify novel gene–disease associations.

EGRET produces per-gene prediction models (with cross-validated R²) and, optionally, gene–trait association results from a summary-statistics-based TWAS.

## Table of Contents

- [Pipeline overview](#pipeline-overview)
- [Quick start with example data](#quick-start-with-example-data)
- [Installation](#installation)
- [Running the full pipeline](#running-the-full-pipeline)
- [Setup](#setup)
- [Building cis models](#building-cis-models)
- [Running Matrix eQTL, GBAT, and _trans_-PCO](#running-matrix-eqtl-gbat-and-trans-pco)
- [Training EGRET models](#training-egret-models)
- [Running EGRET-TWAS](#running-egret-twas)
- [License](#license)
- [Contact](#contact)

## Pipeline overview

EGRET runs in five stages:

1. **Preprocessing** — LD-prune genotypes, filter and regress expression, define cross-validation folds.
2. **Cis model training** — define cis-windows (±500 kb) per gene and fit FUSION cis weights. Required input for GBAT and used as the cis component of the final EGRET model.
3. **Trans-eQTL mapping** — three approaches run in parallel: Matrix eQTL (univariate SNP–gene association), GBAT (polygenic trans-association using cis models), and trans-PCO (WGCNA module-based multivariate association).
4. **EGRET model training** — integrate the trans-eQTL results with cis weights and train per-gene prediction models (lasso, elastic net, BLUP, xtune).
5. **TWAS** — run FUSION-style gene–trait association against GWAS summary statistics. A cis-only TWAS path is also provided as a baseline.

## Quick start with example data

> **Coming soon.** A small example dataset and an end-to-end tutorial walking through the five pipeline stages will be added in a follow-up release.

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

EGRET ships with two binaries:

- `plink2` — bundled at the repo root.
- `gemma-0.98.5-linux-static-AMD64` — bundled under `scripts/`.

Two additional binaries are required by the cis model training stage and are not bundled:

- **PLINK 1.9** — used by the FUSION cis-weight pipeline. Download from https://www.cog-genomics.org/plink/1.9/ and pass the path via `plink1_path`.
- **GCTA (`gcta_nr_robust`)** — used for BLUP weight estimation in FUSION. Available from the FUSION distribution (https://github.com/gusevlab/fusion_twas) and passed via `gcta_path`.

If you want to use a different PLINK 2 build, see https://www.cog-genomics.org/plink/ and pass the path via `plink_path`.

## Running the full pipeline

The end-to-end driver is [`scripts/run_EGRET_workflow.sh`](scripts/run_EGRET_workflow.sh). Edit the parameter block at the top (tissue name, file paths, FDR, number of folds, model types, `output_dir`, etc.) and submit it via SLURM. It chains together every per-stage script documented below.

The per-stage instructions in the rest of this README are for running each stage independently.

## Setup

The setup wrapper preprocesses genotype and expression data into the format EGRET requires. All outputs are stored in `output_dir`.

| Argument                 | Description                                                                                                   |
| ------------------------ | ------------------------------------------------------------------------------------------------------------- |
| `genotypes_file_path`    | Path to PLINK `.bed/.bim/.fam` files containing genotype data.                                                |
| `expression_file_path`   | Gene expression matrix (rows = genes, columns = individuals). Must include a `gene_id` column in ENSG format. |
| `covariates_file_path`   | Covariate file (rows = individuals, columns = covariates to regress out).                                     |
| `tissue`                 | Name of tissue for analysis. Used for labeling and output organization.                                       |
| `individuals_file_path`  | File containing list of individual IDs (one per line, no header).                                             |
| `LD_prune_r2`            | LD pruning threshold (default: `0.9`).                                                                        |
| `plink_path`             | Path to PLINK 2 executable.                                                                                   |
| `genotype_output_prefix` | Prefix for pruned genotype output files.                                                                      |
| `folds`                  | Number of cross-validation folds (e.g., 5).                                                                   |
| `gene_info_file_path`    | File containing gene metadata (gene ID, name, chromosome, start position).                                    |
| `output_dir`             | Directory to store all processed outputs.                                                                     |
| `scripts_dir`            | Path to the EGRET `scripts/` directory.                                                                       |

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
    $output_dir \
    $scripts_dir
```

## Building cis models

After setup, define cis-windows per gene and fit FUSION cis weights. This stage is required before GBAT (which uses cis models as input) and produces the cis component used by the final EGRET model.

Cis-window plink files are generated per gene, then FUSION model fitting is submitted as chunked SLURM jobs (`chunk_size` genes per job; default `2000`).

```
./run_cis_models.sh \
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
```

## Running Matrix eQTL, GBAT, and _trans_-PCO

In EGRET, we implement Matrix eQTL, GBAT, and _trans_-PCO methods to identify _trans_-variants with potential to enhance gene expression models. Each one is run separately and the outputs of each is used in the final training of the genome-wide model.

### Running Matrix eQTL

We utilize the R package Matrix eQTL to conduct pairwise association between all genetic variants and gene expression. It is assumed that the setup scripts were run so that this script can locate their outputs.

```
./MatrixeQTL_scripts.sh \
    $tissue \
    $FDR \
    $folds \
    $gene_info_file_path \
    $genotype_output_prefix \
    $output_dir \
    $scripts_dir
```

### Running GBAT

GBAT uses previously trained cis-expression models to find genome-wide regulators. It expects FUSION-format cis models (produced by [Building cis models](#building-cis-models)).
* If FUSION models were trained elsewhere, weights can also be downloaded from http://gusevlab.org/projects/fusion/

```
./GBAT_scripts.sh \
    $tissue \
    $FDR \
    $folds \
    $gene_info_file_path \
    $cis_model_dir \
    $plink_path \
    $genotypes_file_path \
    $output_dir \
    $scripts_dir
```

`cis_model_dir` is typically `${output_dir}/FUSION/${tissue}/cis/` from the previous stage.

### Running _trans_-PCO

trans-PCO finds the association between a module of co-expressed genes and a variant. `covariate_columns_for_coexpression` is a space-separated list of column names from the covariate file to regress out before computing co-expression modules.

```
./transPCO_scripts.sh \
    $tissue \
    $folds \
    $FDR \
    "$covariate_columns_for_coexpression" \
    $gene_info_file_path \
    $plink_path \
    $output_dir \
    $genotype_output_prefix \
    $scripts_dir
```

## Training EGRET models

This script takes the outputs of the trans-eQTL stages above plus the cis weights and trains genome-wide EGRET models, producing per-gene weights and cross-validated R². We recommend using at least Matrix eQTL and GBAT.

The script also filters cross-mappable reads from alignment errors. By default we supply bed files containing cross-mappable regions for each gene in the hg38 genome build. If you are working in a different genome build, regenerate them with `12.0_make_crossmappable_gene_beds.R` using gene–gene mappability scores from https://doi.org/10.12688/f1000research.17145.2:

```
Rscript 12.0_make_crossmappable_gene_beds.R \
    --gene_info $gene_info_file_path \
    --crossmap_file $crossmap_file \
    --output_dir $output_dir
```

Then train the EGRET models:

```
./train_EGRET_models.sh \
    $tissue \
    $folds \
    $models \
    $plink_path \
    $output_dir \
    $MatrixeQTL_bed_dir \
    $GBAT_bed_dir \
    $transPCO_bed_dir \
    $exclude_crossmap \
    $crossmap_dir \
    $genotype_output_prefix \
    $egret_output_subdir \
    $fusion_models \
    $gemma_path \
    $chunk_size \
    $gene_info_file_path \
    $scripts_dir
```

`exclude_crossmap` is `TRUE`/`FALSE` for whether to filter cross-mappable variants. `egret_output_subdir` is the name of the subdirectory under `xtune_fusion_models/${tissue}/` where the trained models will be written.

## Running EGRET-TWAS

To conduct gene–trait association, we provide a wrapper that takes GWAS summary statistics and the EGRET-trained models and produces gene–trait associations. `gwas_sumstats_csv` is a comma-separated list of summary statistics filenames in `gwas_sumstat_dir`.

```
./run_TWAS_scripts.sh \
    $tissue \
    $output_dir \
    $gene_info_file_path \
    $gwas_sumstat_dir \
    $ld_ref \
    $egret_output_subdir \
    $gwas_sumstats_csv \
    $scripts_dir
```

### Cis-only TWAS (baseline)

EGRET also provides a cis-only TWAS path that uses only the FUSION cis weights as a baseline for comparison against the full EGRET (cis + trans) TWAS results. The arguments are the same.

```
./run_TWAS_cis_scripts.sh \
    $tissue \
    $output_dir \
    $gene_info_file_path \
    $gwas_sumstat_dir \
    $ld_ref \
    $egret_output_subdir \
    $gwas_sumstats_csv \
    $scripts_dir
```

## License

EGRET is released under the MIT License. See [`LICENSE`](LICENSE).

## Contact

Questions, bug reports, and feature requests: please open an issue at [AmariutaLab/EGRET](https://github.com/AmariutaLab/EGRET/issues).
