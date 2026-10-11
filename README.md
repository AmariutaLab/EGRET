# EGRET (Estimating Genome-wide Regulatory Effects on the Transcriptome)

EGRET is a SLURM workflow for building tissue-specific gene expression models from cis and trans regulatory signals. It trains EGRET and cis-only FUSION models, compares their performance, runs TWAS, and can optionally build gene regulatory networks (GRNs).

## 1. Project overview

The workflow moves from cohort-level genotype and expression data to trained models and association results:

```text
setup -> cis models -> GBAT ---------+
      -> Matrix eQTL ----------------+-> EGRET models -> model analysis -> TWAS
      -> trans-PCO ------------------+                                  -> optional GRN
```

Pipeline files are numbered in execution order. For example, `5.0` coordinates the trans-PCO stage and `5.1` through `5.8` perform its individual tasks.

| Path | Contents |
|---|---|
| `environment.yml` | Conda environment for the R and Python dependencies |
| `post_install.R` | Installs xtune, ACAT, and plink2R |
| `scripts/config.sh` | Main configuration for paths, stages, model settings, and cleanup |
| `scripts/run_EGRET_workflow.sh` | Top-level workflow driver and SLURM dependency manager |
| `scripts/1.*` to `scripts/9.*` | Setup, model training, analysis, TWAS, and GRN stages |
| `scripts/validate_egret.sh` | Builds and checks a reduced validation run |
| `scripts/tests/` | Scheduler-free workflow and numerical tests |
| `scripts/LDREF/` | Per-chromosome 1000 Genomes EUR LD reference |
| `grn_analysis_files/` | GRN scripts and the eQTL rsID conversion utility |

The repository also includes GEMMA, the FUSION helper scripts used by the workflow, `qfc.so`, and related statistical routines. GRN analysis is optional and is not included in the default stage list.

## 2. Set up the environment

From the repository root:

```bash
conda env create -f environment.yml
conda activate egret_env
Rscript post_install.R
```

For an existing environment:

```bash
conda env update -n egret_env -f environment.yml
```

`post_install.R` installs the packages Conda does not supply: `xtune`, `ACAT`, and `plink2R`. Run it from the repository root.

`plink2R` is not on CRAN. `post_install.R` installs it from GitHub, or from a local `plink2R-master/plink2R/` directory if one exists. To supply that local copy:

```bash
wget https://github.com/gabraham/plink2R/archive/master.zip
unzip master.zip
```

The script stops if `plink2R` is still missing, because several stages require it.

The workflow assumes a SLURM cluster. Each job activates the environment named by `conda_env`, which defaults to `egret_env`, so it must be visible from compute nodes.

## 3. Required software and data

### Software

The following items are included:

| Software or resource | Location |
|---|---|
| GEMMA 0.98.5 | `scripts/gemma-0.98.5-linux-static-AMD64` |
| `gcta_nr_robust` | `scripts/gcta_nr_robust` |
| FUSION R helpers | `scripts/FUSION.*.R` |
| Davies method shared object | `scripts/qfc.so` |
| 1000 Genomes EUR LD panels | `scripts/LDREF/1000G.EUR.{1..22}.*` |

You must provide:

| Requirement | Config key | Used for |
|---|---|---|
| [PLINK 2](https://www.cog-genomics.org/plink/2.0/) | `plink_path` | Genotype preparation and extraction |
| [PLINK 1.9](https://www.cog-genomics.org/plink/1.9/) | `plink1_path` | Cis-model training and merging the LD reference |
| Conda and SLURM | n/a | Environment management and job submission |

PLINK 2 and `plink2R-master/` are intentionally ignored by Git, so a fresh clone will not contain them. After installing PLINK 1.9, merge the chromosome-specific LD files once:

```bash
bash scripts/LDREF/merge_files.sh /path/to/plink
```

This creates `scripts/LDREF/1000G.EUR.merged.{bed,bim,fam}`, the default `ld_ref` used by FUSION TWAS.

### Data

Set the corresponding paths in `scripts/config.sh`:

| Input | Config key | When required |
|---|---|---|
| PLINK genotype trio | `genotype_prefix` | All model-building stages |
| Normalized expression | `expression` | Setup and model training |
| Sample list | `individuals` | Setup and cis models |
| Covariates | `covariates` | Setup, cis models, and trans-PCO |
| Gene annotation | `gene_info` | All analysis branches |
| GWAS summary statistics | `gwas_sumstat_dir`, `gwas_sumstats` | TWAS and GRN |
| Cross-mappability table | n/a — passed to `0.1_make_crossmappable_gene_beds.R` | Models when `exclude_crossmap="TRUE"` |
| Tissue-matched eQTL table | `eqtl_file` | GRN only |

GTEx data are access-controlled and are not included. The validation utility creates a small test dataset by subsetting the inputs configured in `config.sh`; it is not a standalone public dataset.

The main format requirements are:

- Genotypes must be a PLINK 1 `.bed/.bim/.fam` trio. Chromosomes in the BIM must be numeric and variant IDs must be rsIDs. IIDs must match the sample list; the workflow expects FID to be `0`.
- Expression must be a gzipped, tab-delimited file with a `gene_id` column and sample columns matching the sample list.
- The sample list is headerless with one IID per line. Fold assignment follows file order, so order the samples deliberately.
- Covariates must be tab-delimited with a header, `ID` in column 2, and all columns named by `covariate_columns_for_coexpression`.
- Gene annotation must contain `geneId`, `name`, `geneType`, `#chrom`, `chromStart`, and `chromEnd`. Annotation chromosomes use the `chr` prefix.
- GWAS files must use LDSC format and end in `.sumstats.gz`.
- The GRN eQTL table must contain `rsid`, `gene_id`, and `pval_nominal`.

Stock GTEx eQTL files do not include the required `rsid` column. Add it with:

```bash
python grn_analysis_files/convert_variant_id_to_rsid.py \
  --tissue Whole_Blood \
  --eqtl_file /path/to/Whole_Blood.v8.signif_variant_gene_pairs.txt.gz \
  --bim_file /path/to/genotypes.bim \
  --out_dir /path/to/converted_eqtl
```

## 4. Run the analysis pipeline

### Configure the run

Edit `scripts/config.sh` and, at minimum, review:

- input, binary, and output paths;
- tissue name, folds, FDR, and model lists;
- GWAS traits;
- SLURM account, partition, and log directory;
- selected `stages` and cleanup settings.

Use an absolute `output_dir`. The driver checks required paths and upstream outputs before it submits work.

If cross-mappability filtering is enabled, build its BED files before the `models` stage:

```bash
Rscript scripts/0.1_make_crossmappable_gene_beds.R \
  --gene_info /path/to/gene_annotation.txt.gz \
  --crossmap_file /path/to/cross_mappability.tsv.gz \
  --output_dir /absolute/output/path \
  --crossmap_dir background_mismatches
```

The output path and cross-map directory must match `config.sh`. Set `exclude_crossmap="FALSE"` if this filter is not being used.

### Run a validation analysis

Run this before a full cohort analysis:

```bash
cd scripts
./validate_egret.sh prepare
sbatch --account=my_account --partition=my_partition \
  run_EGRET_workflow.sh --config config_validation.sh
# Wait for the submitted jobs to finish.
./validate_egret.sh check
./validate_egret.sh usage
```

`check` prints a PASS or FAIL for each expected output and exits nonzero if a check fails.

### Submit the full workflow

From the repository root:

```bash
sbatch --account=my_account --partition=my_partition scripts/run_EGRET_workflow.sh
```

The default configuration runs every stage except `grn`. Job IDs are recorded under `output_dir/.job_ids/<run-id>/`, and child-job logs go to `slurm_log_dir/<run-id>/`.

To run another tissue or output directory without editing `config.sh`:

```bash
EGRET_TISSUE=Liver EGRET_OUTPUT_DIR=/absolute/output/Liver \
  sbatch --account=my_account --partition=my_partition scripts/run_EGRET_workflow.sh
```

`expression`, `individuals`, `covariates` and `eqtl_file` are built from the tissue name, so files for that tissue must exist at the configured paths.

The stages are:

| Order | Stage name | What it does |
|---|---|---|
| 1 | `setup` | Prunes genotypes, regresses expression, and creates cross-validation folds |
| 2 | `cis` | Defines cis windows and trains cis-only FUSION models |
| 3 | `matrixeqtl` | Finds trans-eQTLs with Matrix eQTL and applies FDR |
| 4 | `gbat` | Tests associations using cis-imputed expression and applies FDR |
| 5 | `transpco` | Builds coexpression modules and finds trans associations with PCO |
| 6 | `models` | Combines selected variant sources and trains EGRET models |
| 7 | `analysis` | Summarizes EGRET and cis model performance |
| 8 | `twas`, `twas_cis` | Runs TWAS for EGRET and cis-only models |
| 9 | `grn` | Builds the optional GRN analysis from model and TWAS results |

Matrix eQTL and trans-PCO require setup outputs; GBAT also requires the cis models. The driver schedules these branches with the appropriate dependencies.

To run or resume only part of the workflow, override the stage list:

```bash
EGRET_STAGES="models analysis twas twas_cis" \
  sbatch --account=my_account --partition=my_partition scripts/run_EGRET_workflow.sh
```

Valid names are `setup cis matrixeqtl gbat transpco models analysis twas twas_cis grn`. When a stage is skipped, the driver checks any outputs needed by the selected downstream stages. If you submit stages separately, wait for one group to finish before submitting the next.

To include GRN analysis in a complete run:

```bash
EGRET_STAGES="setup cis matrixeqtl gbat transpco models analysis twas twas_cis grn" \
  sbatch --account=my_account --partition=my_partition scripts/run_EGRET_workflow.sh
```

Each stage orchestrator supports `--help` for direct development or recovery work.

## 5. Clean up the project directory

Large EGRET runs create many small intermediate files. The cleanup switches in `scripts/config.sh` default to `TRUE` and remove files only after their final pipeline consumer succeeds.

| Switch | Intermediate files removed |
|---|---|
| `cleanup_cis_plink` | Per-gene cis PLINK files after GBAT |
| `cleanup_cis_working` | Cis working and temporary files after model analysis |
| `cleanup_matrixeqtl_raw` | Raw Matrix eQTL association tables after FDR |
| `cleanup_transpco_raw` | Per-module trans-PCO association, genotype, and expression scratch |
| `cleanup_egret_plink` | Per-gene EGRET PLINK files |
| `cleanup_egret_working` | EGRET model-training scratch |
| `cleanup_egret_z_matrices` | Per-gene z-matrix files |
| `cleanup_bed_lists` | Temporary per-gene variant lists |
| `cleanup_fold_genotypes` | Fold-level genotype matrices after Matrix eQTL |
| `cleanup_fold_expression` | Fold-level expression matrices after EGRET BED preparation |

Keep a switch set to `FALSE` only when you need that stage's intermediates for debugging. A failed gene or stage keeps its scratch files so the failure can be inspected; the `analysis` stage performs a final sweep of the cis and EGRET working directories.

Final model weights, model summaries, FDR-filtered variant sets, TWAS results, GRN results, job IDs, and SLURM logs are retained. Archive those outputs before manually removing an `output_dir`.

For a validation run, inspect it first with:

```bash
cd scripts
./validate_egret.sh usage
```

After validation is complete, the generated `scripts/validation/` directory can be removed. Submission logs in `scripts/batch_submissions/` are also disposable once jobs have been reviewed; production jobs launched through the main driver use the configured `slurm_log_dir`.

## Tests

These checks do not require a live scheduler. The shell tests stub SLURM and Conda, while the numerical test uses R from the active environment:

```bash
bash scripts/tests/test_workflow_wiring.sh
bash scripts/tests/test_grn_workflow.sh
Rscript scripts/tests/test_numerical_equivalence.R
```

## License

EGRET's original code is available under the [MIT License](LICENSE). The bundled GEMMA binary, FUSION scripts, `LDREF/` panel, CompQuadForm routines, and trans-PCO scripts remain under their upstream licenses; see the third-party section of [LICENSE](LICENSE).
