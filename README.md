# EGRET (Estimating Genome-wide Regulatory Effects on the Transcriptome)

EGRET combines cis and trans eQTL evidence to train tissue-specific expression models and run transcriptome-wide association studies (TWAS). The workflow is designed for SLURM and uses R, PLINK/PLINK2, GCTA, GEMMA, FUSION, Matrix eQTL, GBAT, trans-PCO, and xtune.

## Quick start

After supplying the required cohort data and external binaries described below:

```bash
conda env create -f environment.yml && conda activate transTWAS_env
Rscript post_install.R
cd scripts && ./validate_egret.sh prepare
sbatch --account=my_account --partition=my_partition run_EGRET_workflow.sh --config config_validation.sh
./validate_egret.sh check
```

The validation run uses a reduced subset of your configured inputs. See
[Installation](#installation), [Input data formats](#input-data-formats), and
[Start with the validation run](#start-with-the-validation-run) before launching a
production analysis.

## Pipeline at a glance

Numbered filenames follow the workflow: the leading digit is the pipeline stage,
`N.0` is that stage's orchestrator, and later `N.x` files are its workers.

| Config stage | Orchestrator | Workers and purpose |
|---|---|---|
| Manual preparation | — | `0.1_make_crossmappable_gene_beds.R` builds optional cross-mappability BEDs |
| `setup` | `1.0_setup_genotype_and_expression.sh` | `1.1`–`1.3`: prune genotypes, regress expression, create folds |
| `cis` | `2.0_run_cis_models.sh` | `2.1`–`2.2`: build cis windows and train FUSION weights |
| `matrixeqtl` | `3.0_run_MatrixeQTL_stage.sh` | `3.1`–`3.2`: scan trans-eQTLs and apply FDR |
| `gbat` | `4.0_run_GBAT_stage.sh` | `4.1`–`4.3`: impute cis expression, test associations, apply FDR |
| `transpco` | `5.0_run_transPCO_stage.sh` | `5.1`–`5.8`: build modules, run Matrix eQTL/PCO, apply FDR |
| `models` | `6.0_train_EGRET_models.sh` | `6.1`–`6.4`: assemble variant sets and train EGRET models |
| `analysis` | `7.0_run_model_analysis.sh` | `7.1`–`7.2`: summarize EGRET and cis models |
| `twas`, `twas_cis` | `8.0_run_TWAS_EGRET.sh`, `8.1_run_TWAS_cis.sh` | `8.2`–`8.4`: build position files and run per-trait TWAS |
| `grn` | `9.0_submit_GRN_analysis.sh` | `9.1` dispatches GRN steps `9.2`–`9.5` |

## What the repository provides, and what you must supply

Bundled in the checkout:

| Item | Path |
|---|---|
| GEMMA 0.98.5 (linux static AMD64) | `scripts/gemma-0.98.5-linux-static-AMD64` |
| Davies' method shared object | `scripts/qfc.so` |
| 1000G EUR LD reference, per chromosome | `scripts/LDREF/1000G.EUR.{1..22}.{bed,bim,fam}` |

You must supply these yourself — none are in the repository:

| Item | Config key | Where to get it |
|---|---|---|
| PLINK 2 | `plink_path` | <https://www.cog-genomics.org/plink/2.0/> |
| PLINK 1.9 | `plink1_path` | <https://www.cog-genomics.org/plink/1.9/> |
| `gcta_nr_robust` | `gcta_path` | FUSION distribution |
| `plink2R` R package source | — | <https://github.com/gabraham/plink2R> |
| Genotypes, expression, covariates, sample list | `genotype_prefix`, `expression`, `covariates`, `individuals` | your cohort (GTEx v8 is access-controlled) |
| Gene annotation | `gene_info` | GTEx v8 gene model |
| GWAS summary statistics | `gwas_sumstat_dir` | LDSC `.sumstats.gz` format |
| Cross-mappability table | `crossmap_file` | Saha & Battle 2018 |
| GTEx eQTL table (GRN stage only) | `eqtl_file` | GTEx v8 `signif_variant_gene_pairs` |

`plink_path` and `plink2R-master/` are `.gitignore`d, so they are absent from a fresh clone
even though the default config expects them at the repository root.

## Installation

```bash
conda env create -f environment.yml     # creates transTWAS_env
conda activate transTWAS_env
Rscript post_install.R                  # xtune, ACAT, plink2R
```

`post_install.R` installs `plink2R` from `plink2R-master/` if present and otherwise from
GitHub; it **stops with an error** rather than continuing if neither works, because a
silently missing `plink2R` produces an environment that looks complete and then fails five
stages into a submitted run. For an existing checkout, `conda env update -n transTWAS_env
-f environment.yml`.

Then build the merged LD reference FUSION needs. `FUSION.assoc_test_trans.R` reads a single
PLINK prefix, so the per-chromosome files cannot be used directly:

```bash
bash scripts/LDREF/merge_files.sh /path/to/plink     # PLINK 1.9
```

That writes `scripts/LDREF/1000G.EUR.merged.*`, which is where `ld_ref` already points.

## Input data formats

The top-level driver validates required paths before submission. The following schema
requirements must also hold for downstream joins and chromosome-wise processing.

**Genotypes** (`genotype_prefix`) — PLINK 1 binary trio.
- The `.bim` chromosome column must be **numeric** (`1`, not `chr1`) and cover a contiguous
  `1..max` with no empty chromosome. Four scripts loop `1:max(bim$V1)`.
- Variant IDs must be **rsIDs**. Trans-SNP lists and the GRN eQTL join both match on them,
  so GTEx-style `chr_pos_ref_alt_b38` IDs will silently fail to join.
- FID is assumed to be literal `0`; IIDs must match the `individuals` file.

**Expression** (`expression`) — gzipped, tab-delimited, with a header.
- Must contain a column named exactly `gene_id`, plus one column per sample named exactly as
  in `individuals`. A GTEx `.bed.gz` (`#chr start end gene_id GTEX-…`) works unmodified.
- Only `protein_coding`, `lincRNA`, and `antisense` genes survive, per `gene_info`. Ensembl
  version suffixes must match between the two files.

**Individuals** (`individuals`) — headerless, one IID per line. Note that
cross-validation folds are contiguous blocks of **file order**, not shuffled: order this
file deliberately.

**Covariates** (`covariates`) — tab-delimited with a header. Column 2 must be named `ID`;
columns 3 onward are treated as covariates. It must also contain every column named in
`covariate_columns_for_coexpression`.

**Gene annotation** (`gene_info`) — tab-delimited (gzip fine), with columns `geneId`,
`name`, `geneType`, `#chrom` (**`chr`-prefixed**), `chromStart`, `chromEnd`.

**GWAS summary statistics** — LDSC format; filenames must end in `.sumstats.gz`.

**eQTL table** (`eqtl_file`, GRN stage only) — must carry an `rsid` column, which stock GTEx
files do not. Generate it once:

```bash
python grn_analysis_files/convert_variant_id_to_rsid.py \
  --tissue Whole_Blood \
  --eqtl_file .../Whole_Blood.v8.signif_variant_gene_pairs.txt.gz \
  --bim_file <genotypes>.bim \
  --out_dir grn_analysis_files/GTEx_Analysis_v8_eQTL
```

The converter matches chromosome, position, and the unordered allele pair. This tolerates
ref/alt reversal without assigning the first rsID at a multi-allelic position.

## Start with the validation run

Before committing a full-scale run, exercise every stage at reduced scale. This is the
fastest way to find out whether your inputs and binaries are wired up correctly:

```bash
cd scripts
./validate_egret.sh prepare      # ~1.8k genes, ~80k SNPs -> validation/inputs
sbatch --account=my_account --partition=my_partition \
  run_EGRET_workflow.sh --config config_validation.sh
./validate_egret.sh check        # asserts the run produced what it should
./validate_egret.sh usage        # inode and byte footprint
```

`check` prints a `PASS`/`FAIL` line per assertion and exits non-zero if any failed. The
whole thing is roughly 55 jobs; the longest single job is about 11 minutes. Discard it with
`rm -rf scripts/validation`.

`prepare` reads `genotype_prefix`, `expression`, and `gene_info` out of `config.sh`, so
those three keys must stay plain double-quoted assignments at column 0.

## Sample data

`validate_egret.sh prepare` creates a small end-to-end validation input with approximately
1,800 genes and 80,000 variants, but it does so by subsetting the cohort configured in
`config.sh`. GTEx genotypes are access-controlled, so this validation input cannot be
redistributed as public sample data.

A fixed-seed simulated dataset with the same schemas is planned. Until it is available,
users must supply their own authorized genotype, expression, covariate, and annotation
inputs before running either validation or production analyses.

## Configuration and full workflow

Copy and edit [`scripts/config.sh`](scripts/config.sh). It contains parameters only: input paths, binary paths, tissue and model settings, eQTL source directories, cleanup switches, selected stages, GWAS traits, conda environment, and SLURM account/partition.

Paths that live inside the checkout are written relative to `_scripts_dir` / `_repo_dir`,
which resolve from `EGRET_SCRIPTS_DIR` when the driver sets it and otherwise from the config
file's own location — so a clone at any path works unedited, and a config copied elsewhere
still finds the bundled GEMMA and LD reference. Keep `output_dir` absolute: a bare relative
path resolves against the directory `sbatch` was invoked from, which is how one config once
produced two divergent output trees.

Before submitting anything, the driver checks that every input file, directory, and binary
needed by the selected stages actually exists, and reports all problems at once:

```text
ERROR: configuration refers to paths that do not exist (scripts/config.sh):
  - gene_info: no such file: /nope/GTEx_V8.txt.gz
  - gcta_path: not an executable and not on PATH: /nope/gcta_nr_robust
  - ld_ref: missing .bed .bim .fam for prefix /nope/1000G.EUR.merged
```

Raw genotypes and raw expression are checked only for `setup`; individuals and covariates
are also checked for `cis`. Configured GWAS trait files are checked for TWAS/GRN stages,
cross-mappability BEDs are checked for filtered model runs, and the eQTL table is checked
only when `grn` is selected.

Submit the driver from the repository root:

```bash
sbatch --account=my_account --partition=my_partition scripts/run_EGRET_workflow.sh
```

The top-level `#SBATCH` account and partition are bootstrap defaults because SLURM directives cannot expand values loaded from `config.sh`. The command-line values above override them. Every child submission receives `slurm_account`, `slurm_partition`, its log path, and its dependency explicitly.

To use a different configuration without editing the default:

```bash
EGRET_CONFIG=/path/to/my_config.sh sbatch --account=my_account --partition=my_partition scripts/run_EGRET_workflow.sh
```

For an interactive invocation, the equivalent named option is `scripts/run_EGRET_workflow.sh --config /path/to/my_config.sh`.

The dependency graph is:

```text
1.0 setup → 2.0 cis → 4.0 GBAT ┐
          → 3.0 Matrix eQTL ────┼→ 6.0 models → 7.0 analysis → 8.0 EGRET TWAS ┐
          → 5.0 trans-PCO ──────┘                              8.1 cis TWAS ├→ 9.0 GRN
```

Job IDs are stored under `output_dir/.job_ids/<run-id>/`. Orchestrators receive explicit current-run job-ID paths; IDs from earlier runs are never read. Model analysis waits directly for both current-run cis and EGRET model chunks, and both TWAS drivers use `afterok`. If a stage is omitted from `stages`, the driver validates its outputs only when a selected downstream stage consumes them.

## Cross-mappability filtering

When `exclude_crossmap="TRUE"`, the bed builder drops SNPs falling near genes that
cross-map to the target gene — sequence similarity otherwise manufactures trans signal.
**The BED files are not bundled**, and they are not generated by any stage. Build them once
per `output_dir`, before the `models` stage:

```bash
cd scripts
Rscript 0.1_make_crossmappable_gene_beds.R \
  --gene_info "$gene_info" \
  --crossmap_file "$crossmap_file" \
  --output_dir "$output_dir" \
  --crossmap_dir background_mismatches
```

`--crossmap_dir` must match `crossmap_dir` in your config. When cross-mappability
filtering is enabled, missing BED files are a hard preflight error; the bed builder repeats
the check when run directly. Set `exclude_crossmap="FALSE"` to run without this filter.

## Resuming, and running one stage at a time

`stages` selects which stages run; `EGRET_STAGES` overrides it without editing the config:

```bash
EGRET_STAGES="models analysis twas twas_cis" \
  sbatch --account=my_account --partition=my_partition scripts/run_EGRET_workflow.sh
```

Skipping a stage is safe. For every omitted stage that feeds a selected stage, the driver
asserts that the required outputs already exist — for example, `models` with `Cis` selected
requires existing cis BED windows, while `analysis` requires cis and EGRET weights. Omitted
stages unrelated to the requested work are not checked, so `EGRET_STAGES="cis"` and other
focused runs do not require outputs from later stages.

Valid stage names: `setup cis matrixeqtl gbat transpco models analysis twas twas_cis grn`.
`grn` is opt-in and must be named explicitly.

One caveat when invoking a worker directly with `sbatch`: their `#SBATCH --output`
directives default to a relative `batch_submissions/`, which is `.gitignore`d. Create it
first (`mkdir -p scripts/batch_submissions`) or pass `--output` yourself. Jobs submitted by
the driver are unaffected — it always passes an explicit log path.

The TWAS stage is idempotent: a trait whose `.dat` already has results is skipped, so
re-running to fill gaps is cheap.

## Optional GRN analysis

Add `grn` to the stage list after `analysis twas twas_cis`. For example:

```bash
EGRET_STAGES="setup cis matrixeqtl gbat transpco models analysis twas twas_cis grn" \
  sbatch --account=my_account --partition=my_partition scripts/run_EGRET_workflow.sh
```

Both `config.sh` and `config_validation.sh` honour `EGRET_STAGES`.

Set `eqtl_file` to a tissue-matched, tab-delimited eQTL table containing `rsid`,
`gene_id`, and `pval_nominal` — see [Input data formats](#input-data-formats) for how to
produce the `rsid` column, which stock GTEx files lack. `grn_pval_threshold` controls the
trans-component significance cutoff used for both network selection and module summaries;
it defaults to `0.01` when omitted. The GRN dispatcher waits for both
TWAS dispatchers and then makes the GRN worker depend on every current-run EGRET and cis
trait job. Its IDs are recorded as `grn_driver.jids` and `grn.jids`.

## Running stages directly

Shell entry points use named options, reject unknown or positional arguments, and support
`--help`. The top-level workflow supplies their complete interfaces; direct invocation is
mainly useful when developing or resuming one stage.

```bash
scripts/6.0_train_EGRET_models.sh --help
```

Common options include `--gene_info`, `--expression`, `--covariates`,
`--output_dir`, `--genotype_prefix`, `--egret_subdir`, `--bed_models`, and
`--fusion_models`. Binary paths use `--plink_path`, `--plink1_path`,
`--gcta_path`, and `--gemma_path`.

The run directory records every submission tier: stage drivers, chunk jobs, analysis,
per-trait TWAS jobs, and the optional GRN chain. This record belongs to the current run
only and is stored under `output_dir/.job_ids/<run-id>/`.

## Reproducibility and validation

EGRET model chunks are assigned from the stable gene order in the regressed-expression file using zero-based `--start_index` and `--chunk_size`. They do not enumerate a PLINK directory that cleanup jobs may shrink.

The cis-FUSION comparison arm includes FUSION's `top1` fallback, while EGRET deliberately
does not: EGRET remains a genome-wide multi-SNP method rather than falling back to one cis
variant. Consequently, direct R-squared comparisons should either report this asymmetry or
provide a sensitivity analysis excluding genes for which `top1` is the best cis model.

For every requested Matrix eQTL, GBAT, and trans-PCO source, the bed builder reports the
number of files found, SNPs listed, SNPs retained in the final per-gene PLINK output, and
genes with at least one retained SNP. A fold fails when a requested source contributes zero
SNPs after PLINK extraction and optional cross-map filtering.

Cleanup controls remove only stage-specific intermediates after their final consumer succeeds. Generated output trees are otherwise retained.

Three test suites run without a scheduler — they stub `sbatch` and `conda`:

```bash
bash scripts/tests/test_workflow_wiring.sh    # workflow wiring, dependencies, chunking (needs ripgrep)
bash scripts/tests/test_grn_workflow.sh       # GRN dispatcher and the three GRN scripts
Rscript scripts/tests/test_numerical_equivalence.R # numerical equivalence of the optimizations
```

## Porting to another cluster

The pipeline assumes SLURM and conda. There is no `module load` anywhere; batch workers
bootstrap conda themselves by sourcing `$(conda info --base)/etc/profile.d/conda.sh`,
falling back to `~/.bashrc`.

- **Account and partition.** `slurm_account` and `slurm_partition` in the config are passed
  explicitly to every child submission. The top-level `#SBATCH` directives are only
  bootstrap defaults — SLURM reads them before `config.sh` is sourced — so override them on
  the command line, as in the `sbatch --account=… --partition=…` examples above.
- **Per-job resources are not parameterized.** Each worker hardcodes its own `--mem` and
  `-t` (from 1 G/30 min for the GRN dispatcher up to 64 G for Matrix eQTL and `-t 48:00:00`
  for TWAS association). If your cluster caps walltime below 48 h, edit those directives.
- **The parent job does real work.** `setup`, cis-window setup, WGCNA clustering, and
  per-chromosome genotype preparation run *inside* `run_EGRET_workflow.sh` itself, which is
  capped at 30 G / 8 h. At full scale on a new cluster this is the first thing likely to hit
  a wall.
- **Disk is usually the binding constraint, and it is inodes rather than bytes.** One
  tissue's `output_dir` reaches roughly 2 M files. Every `cleanup_*` switch defaults to
  `TRUE` for that reason; set individual ones to `FALSE` only when debugging that stage.

## Performance notes

Expression residualization shares one covariate decomposition across genes. Chunked cis
and EGRET workers decompress expression once per chunk, while trans-PCO workers reuse
fold-level genotype and expression inputs. Disk use is usually constrained by inode count;
the cleanup switches in `config.sh` remove intermediates after their last consumer.

## License

EGRET's original code is available under the [MIT License](LICENSE). Bundled third-party
software, reference data, and adapted statistical routines—including `FUSION.*.R`,
`davies.R`, `liu.R`, `liumod.R`, and `Modified*.R`—remain subject to their
respective upstream terms and are not relicensed by the EGRET MIT license.
