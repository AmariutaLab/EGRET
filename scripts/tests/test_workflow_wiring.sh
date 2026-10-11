#!/bin/bash
# --- test_workflow_wiring.sh -----------------------------------------------
# Asserts the workflow wires up correctly without running any science: stage selection,
# SLURM dependency edges, gene chunking, and flag pass-through. Uses a temporary tree.
#
#   in : nothing -- self-contained, builds its own fixtures
#   out: PASS/FAIL on stdout
#   run: bash scripts/tests/test_workflow_wiring.sh   (needs ripgrep)
# --------------------------------------------------------------------------
set -euo pipefail
scripts_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
test_root=$(mktemp -d)
trap 'rm -rf "$test_root"' EXIT
export TEST_ROOT=$test_root

fail() { echo "FAIL: $*" >&2; exit 1; }
assert_grep() { rg -q -- "$1" "$2" || fail "missing '$1' in $2"; }

mkdir -p "$test_root/bin" "$test_root/out/genotype_files" \
    "$test_root/out/expression_files" "$test_root/out/fold_0_info/T/train_individuals.txt.dummy" \
    "$test_root/out/FUSION/T/cis" "$test_root/out/bed_files/T/cis" \
    "$test_root/out/MatrixeQTL/T/fold_0/meqtl" "$test_root/out/GBAT/T/fold_0/gbat" \
    "$test_root/out/transPCO/T/fold_0/trans"
rmdir "$test_root/out/fold_0_info/T/train_individuals.txt.dummy"
printf '1\trs1\t0\t100\tA\tG\n' > "$test_root/out/genotype_files/geno.bim"
printf 's1\n' > "$test_root/out/fold_0_info/T/train_individuals.txt"
printf 'weight\n' > "$test_root/out/FUSION/T/cis/g.wgt.RDat"
printf '1\t100\t200\tg\t0\t+\n' > "$test_root/out/bed_files/T/cis/g.bed"
printf 'rs1\n' > "$test_root/out/MatrixeQTL/T/fold_0/meqtl/g.txt"
printf 'rs1\n' > "$test_root/out/GBAT/T/fold_0/gbat/g.txt"
printf 'rs1\n' > "$test_root/out/transPCO/T/fold_0/trans/g.txt"
printf 'gene_id\ts1\ng1\t1\ng2\t2\ng3\t3\n' | gzip -c > "$test_root/out/expression_files/T_expression_regressed.txt.gz"

mkdir -p "$test_root/fake_conda/etc/profile.d"
cat > "$test_root/bin/conda" <<'CONDA_BIN'
#!/bin/bash
if [[ "${1:-}" == info && "${2:-}" == --base ]]; then
    echo "$TEST_ROOT/fake_conda"
fi
CONDA_BIN
cat > "$test_root/fake_conda/etc/profile.d/conda.sh" <<'CONDA_HOOK'
conda() { return 0; }
CONDA_HOOK
printf '1000\n' > "$test_root/sbatch.state"
cat > "$test_root/bin/sbatch" <<'STUB'
#!/bin/bash
set -euo pipefail
n=$(cat "$SBATCH_STATE")
n=$((n + 1))
printf '%s\n' "$n" > "$SBATCH_STATE"
printf '%s\t%s\n' "$n" "$*" >> "$SBATCH_LOG"
printf '%s\n' "$n"
STUB
chmod +x "$test_root/bin/conda" "$test_root/bin/sbatch"
export SBATCH_STATE="$test_root/sbatch.state" SBATCH_LOG="$test_root/sbatch.log"

cat > "$test_root/config.sh" <<'CONFIG'
genotype_prefix="$TEST_ROOT/raw_geno"
expression="$TEST_ROOT/raw_expression"
individuals="$TEST_ROOT/individuals"
covariates="$TEST_ROOT/covariates"
gene_info="$TEST_ROOT/gene_info"
plink_path="$TEST_ROOT/bin/plink"
plink1_path="$TEST_ROOT/bin/plink1"
gcta_path="$TEST_ROOT/bin/gcta"
gemma_path="$TEST_ROOT/bin/gemma"
tissue="T"
output_dir="$TEST_ROOT/out"
genotype_output_prefix="geno"
egret_subdir="custom"
folds="0"
ld_prune_r2="0.9"
fdr="0.1"
chunk_size="2"
bed_models="Cis,MatrixeQTL,GBAT,transPCO"
fusion_models="xtune,lasso"
matrixeqtl_bed_dir="meqtl"
gbat_bed_dir="gbat"
transpco_bed_dir="trans"
exclude_crossmap="${EGRET_EXCLUDE_CROSSMAP:-FALSE}"
crossmap_dir="none"
crossmap_file="unused"
covariate_columns_for_coexpression="PC1"
stages="${EGRET_STAGES:-models analysis twas twas_cis}"
cleanup_cis_plink="FALSE"
cleanup_cis_working="FALSE"
cleanup_matrixeqtl_raw="FALSE"
cleanup_transpco_raw="FALSE"
cleanup_egret_plink="FALSE"
cleanup_egret_working="FALSE"
cleanup_egret_z_matrices="FALSE"
cleanup_bed_lists="FALSE"
cleanup_fold_genotypes="FALSE"
cleanup_fold_expression="FALSE"
gwas_sumstat_dir="$TEST_ROOT/gwas"
ld_ref="$TEST_ROOT/ldref"
gwas_sumstats="trait.sumstats.gz"
conda_env="fake"
slurm_account="acct"
slurm_partition="part"
slurm_log_dir="$TEST_ROOT/logs"
CONFIG
# The selected stages need gene_info, PLINK2, GEMMA, GWAS, and LD inputs. The configured
# setup/cis-only paths (raw inputs, individuals, covariates, PLINK1, and GCTA) deliberately
# do not exist: a regression to unconditional preflight validation will fail this run.
mkdir -p "$test_root/bin" "$test_root/gwas" "$test_root/logs"
touch "$test_root/gene_info"
touch "$test_root/gwas/trait.sumstats.gz"
for stub in plink gemma; do
    printf '#!/bin/bash\nexit 0\n' > "$test_root/bin/$stub"; chmod +x "$test_root/bin/$stub"
done
for ext in bed bim fam; do touch "$test_root/ldref.$ext"; done

# Selected-stage preflight catches trait typos and missing cross-map resources before any
# jobs are submitted. Neither setting is required by an unrelated stage.
rm "$test_root/gwas/trait.sumstats.gz"
if PATH="$test_root/bin:$PATH" EGRET_CONFIG="$test_root/config.sh" \
    EGRET_STAGES="twas" EGRET_RUN_ID=missingGwas \
    "$scripts_dir/run_EGRET_workflow.sh" > "$test_root/missing_gwas.out" 2>&1; then
    fail "missing configured GWAS file unexpectedly passed preflight"
fi
assert_grep 'gwas_sumstats: no such file:.*trait.sumstats.gz' "$test_root/missing_gwas.out"
touch "$test_root/gwas/trait.sumstats.gz"
if PATH="$test_root/bin:$PATH" EGRET_CONFIG="$test_root/config.sh" \
    EGRET_STAGES="models" EGRET_EXCLUDE_CROSSMAP=TRUE EGRET_RUN_ID=missingCrossmap \
    "$scripts_dir/run_EGRET_workflow.sh" > "$test_root/missing_crossmap.out" 2>&1; then
    fail "missing cross-mappability BEDs unexpectedly passed preflight"
fi
assert_grep 'cross-mappability BEDs: no files matching' "$test_root/missing_crossmap.out"

mkdir -p "$test_root/out/.job_ids/stale"; printf '999\n' > "$test_root/out/.job_ids/stale/models.jids"
PATH="$test_root/bin:$PATH" EGRET_CONFIG="$test_root/config.sh" \
    EGRET_STAGES="models analysis twas twas_cis" EGRET_RUN_ID=runA \
    "$scripts_dir/run_EGRET_workflow.sh" > "$test_root/workflow.out"
[[ $(cat "$test_root/out/.job_ids/runA/models.jids") == 1002:1003 ]] || fail "wrong model IDs"
[[ $(cat "$test_root/out/.job_ids/runA/beds.jids") == 1001 ]] || fail "wrong bed IDs"
[[ $(cat "$test_root/out/.job_ids/runA/analysis.jids") == 1004 ]] || fail "wrong analysis ID"
assert_grep '1002.*--dependency=afterok:1001.*6.3_train_EGRET_chunk.sh' "$test_root/sbatch.log"
assert_grep '1004.*--dependency=afterok:1002:1003.*7.0_run_model_analysis.sh' "$test_root/sbatch.log"
assert_grep '1005.*--dependency=afterok:1004.*8.0_run_TWAS_EGRET.sh' "$test_root/sbatch.log"
assert_grep '1006.*--dependency=afterok:1004.*8.1_run_TWAS_cis.sh' "$test_root/sbatch.log"
assert_grep '8.0_run_TWAS_EGRET.sh.*--job_id_file.*twas.jids' "$test_root/sbatch.log"
assert_grep '8.1_run_TWAS_cis.sh.*--job_id_file.*twas_cis.jids' "$test_root/sbatch.log"
! rg -q '999' "$test_root/sbatch.log" || fail "stale job ID was consumed"

# Analysis consumes cis weights even when EGRET model training is omitted, so it must depend
# directly on current-run cis chunks rather than receiving that ordering only through models.
# Supply the additional inputs and binaries that the newly selected cis stage really uses.
touch "$test_root/individuals" "$test_root/covariates"
for stub in plink1 gcta; do
    printf '#!/bin/bash\nexit 0\n' > "$test_root/bin/$stub"
    chmod +x "$test_root/bin/$stub"
done
mkdir -p "$test_root/out/EGRET_models/T/custom"
printf 'weight\n' > "$test_root/out/EGRET_models/T/custom/g.wgt.RDat"
cat > "$test_root/bin/Rscript" <<'CIS_RSCRIPT'
#!/bin/bash
set -euo pipefail
if [[ "${1:-}" == */2.1_setup_cis_windows.R ]]; then
    mkdir -p "$TEST_ROOT/out/plink_results/T/cis"
    printf 'g.bed\n' > "$TEST_ROOT/out/plink_results/T/cis/cis_genes.txt"
    printf 'bed\n' > "$TEST_ROOT/out/plink_results/T/cis/g.bed"
    printf '1\trs1\t0\t100\tA\tG\n' > "$TEST_ROOT/out/plink_results/T/cis/g.bim"
    printf '0\ts1\t0\t0\t0\t-9\n' > "$TEST_ROOT/out/plink_results/T/cis/g.fam"
fi
CIS_RSCRIPT
chmod +x "$test_root/bin/Rscript"
PATH="$test_root/bin:$PATH" EGRET_CONFIG="$test_root/config.sh" \
    EGRET_STAGES="cis analysis" EGRET_RUN_ID=runB \
    "$scripts_dir/run_EGRET_workflow.sh" > "$test_root/workflow_cis_analysis.out"
[[ $(cat "$test_root/out/.job_ids/runB/cis.jids") == 1007 ]] || fail "wrong cis ID"
[[ $(cat "$test_root/out/.job_ids/runB/analysis.jids") == 1008 ]] || fail "wrong cis-analysis ID"
assert_grep '1008.*--dependency=afterok:1007.*7.0_run_model_analysis.sh' "$test_root/sbatch.log"

# If analysis is intentionally omitted, TWAS must wait on the current producer directly.
# Reset the mock scheduler range so these assertions do not disturb the original ID sequence.
mkdir -p "$test_root/out/results_sumstats/T"
printf 'gene\tr2\n' > "$test_root/out/results_sumstats/T/custom.txt"
printf 'gene\tr2\n' > "$test_root/out/results_sumstats/T/cis.txt"
printf '2000\n' > "$test_root/sbatch.state"
PATH="$test_root/bin:$PATH" EGRET_CONFIG="$test_root/config.sh" \
    EGRET_STAGES="models twas" EGRET_RUN_ID=runModelsTwas \
    "$scripts_dir/run_EGRET_workflow.sh" > "$test_root/workflow_models_twas.out"
assert_grep '2004.*--dependency=afterok:2002:2003.*8.0_run_TWAS_EGRET.sh' "$test_root/sbatch.log"
PATH="$test_root/bin:$PATH" EGRET_CONFIG="$test_root/config.sh" \
    EGRET_STAGES="cis twas_cis" EGRET_RUN_ID=runCisTwas \
    "$scripts_dir/run_EGRET_workflow.sh" > "$test_root/workflow_cis_twas.out"
assert_grep '2006.*--dependency=afterok:2005.*8.1_run_TWAS_cis.sh' "$test_root/sbatch.log"
printf '1008\n' > "$test_root/sbatch.state"

# A focused analysis-only resume must not require setup/source-stage artifacts or binaries.
rm -rf "$test_root/out/genotype_files" "$test_root/out/expression_files" \
    "$test_root/out/fold_0_info" "$test_root/out/MatrixeQTL" "$test_root/out/GBAT" \
    "$test_root/out/transPCO" "$test_root/out/bed_files"
rm -f "$test_root/gene_info" "$test_root/bin/plink" "$test_root/bin/gemma" \
    "$test_root/bin/plink1" "$test_root/bin/gcta"
PATH="$test_root/bin:$PATH" EGRET_CONFIG="$test_root/config.sh" \
    EGRET_STAGES="analysis" EGRET_RUN_ID=runC \
    "$scripts_dir/run_EGRET_workflow.sh" > "$test_root/workflow_analysis_only.out"
[[ $(cat "$test_root/out/.job_ids/runC/analysis.jids") == 1009 ]] || fail "wrong analysis-only ID"
assert_grep '^1009.*7.0_run_model_analysis.sh' "$test_root/sbatch.log"
! rg -q '^1009.*--dependency' "$test_root/sbatch.log" ||
    fail "analysis-only resume received an unrelated dependency"

# Child orchestrators persist every internal job tier in the current run directory.
printf '#!/bin/bash\nexit 0\n' > "$test_root/bin/Rscript"
chmod +x "$test_root/bin/Rscript"
mkdir -p "$test_root/out/results_sumstats/T" "$test_root/out/transPCO/T/fold_0/modules"
printf 'gene\tr2\n' > "$test_root/out/results_sumstats/T/custom.txt"
printf 'gene\tr2\n' > "$test_root/out/results_sumstats/T/cis.txt"
printf 'g1\n' > "$test_root/out/transPCO/T/fold_0/modules/module_1.txt"
PATH="$test_root/bin:$PATH" "$scripts_dir/8.0_run_TWAS_EGRET.sh" --tissue T --output_dir "$test_root/out" \
    --gene_info genes --gwas_sumstat_dir "$test_root/gwas" --ld_ref ld --egret_subdir custom \
    --gwas_traits a.sumstats.gz,b.sumstats.gz --conda_env fake --slurm_account acct \
    --slurm_partition part --slurm_log_dir "$test_root/logs" \
    --job_id_file "$test_root/out/.job_ids/runA/twas.jids"
PATH="$test_root/bin:$PATH" "$scripts_dir/8.1_run_TWAS_cis.sh" --tissue T --output_dir "$test_root/out" \
    --gene_info genes --gwas_sumstat_dir "$test_root/gwas" --ld_ref ld \
    --gwas_traits a.sumstats.gz,b.sumstats.gz --conda_env fake --slurm_account acct \
    --slurm_partition part --slurm_log_dir "$test_root/logs" \
    --job_id_file "$test_root/out/.job_ids/runA/twas_cis.jids"
[[ $(tr ':' '\n' < "$test_root/out/.job_ids/runA/twas.jids" | wc -l) -eq 2 ]] || fail "EGRET TWAS IDs not recorded"
[[ $(tr ':' '\n' < "$test_root/out/.job_ids/runA/twas_cis.jids" | wc -l) -eq 2 ]] || fail "cis TWAS IDs not recorded"
PATH="$test_root/bin:$PATH" "$scripts_dir/5.0_run_transPCO_stage.sh" --tissue T --folds 0 --fdr 0.1 \
    --covariate_columns PC1 --gene_info genes --plink_path plink --output_dir "$test_root/out" \
    --genotype_prefix geno --conda_env fake --slurm_account acct --slurm_partition part \
    --slurm_log_dir "$test_root/logs" --matrix_job_id_file "$test_root/out/.job_ids/runA/transpco_matrix.jids" \
    --pco_job_id_file "$test_root/out/.job_ids/runA/transpco_pco.jids" \
    --job_id_file "$test_root/out/.job_ids/runA/transpco.jids"
[[ -s "$test_root/out/.job_ids/runA/transpco_matrix.jids" ]] || fail "transPCO MatrixeQTL IDs not recorded"
[[ -s "$test_root/out/.job_ids/runA/transpco_pco.jids" ]] || fail "transPCO PCO IDs not recorded"
[[ -s "$test_root/out/.job_ids/runA/transpco.jids" ]] || fail "transPCO final ID not recorded"

# Stable expression-based chunks remain disjoint after cleanup removes chunk-zero PLINK files.
chunk="$test_root/chunk"; mkdir -p "$chunk/bin" "$chunk/out/expression_files" "$chunk/out/fold_0_info/T" \
    "$chunk/out/plink_results/T/custom/fold_0"
printf 'gene_id\ts1\ts2\ng1\t1\t2\ng2\t2\t3\ng3\t3\t4\ng4\t4\t5\n' | gzip -c > "$chunk/out/expression_files/expression.gz"
printf 's1\ns2\n' > "$chunk/out/fold_0_info/T/train_individuals.txt"
for gene in g1 g2 g3 g4; do
    printf 'bed\n' > "$chunk/out/plink_results/T/custom/fold_0/$gene.bed"
    printf '1\trs1\t0\t100\tA\tG\n1\trs2\t0\t200\tC\tT\n' > "$chunk/out/plink_results/T/custom/fold_0/$gene.bim"
    printf '0\ts1\t0\t0\t0\t-9\n0\ts2\t0\t0\t0\t-9\n' > "$chunk/out/plink_results/T/custom/fold_0/$gene.fam"
done
cp "$test_root/bin/conda" "$chunk/bin/conda"
cat > "$chunk/bin/zcat" <<'ZCAT'
#!/bin/bash
printf '1\n' >> "$ZCAT_LOG"
exec /usr/bin/zcat "$@"
ZCAT
cat > "$chunk/bin/plink" <<'PLINK'
#!/bin/bash
set -euo pipefail
out=""; bfile=""
while (($#)); do case $1 in --out) out=$2; shift 2;; --bfile) bfile=$2; shift 2;; *) shift;; esac; done
mkdir -p "$(dirname "$out")"
cp "${bfile}.bim" "${out}.bim"; cp "${bfile}.fam" "${out}.fam"; printf 'bed\n' > "${out}.bed"
PLINK
cat > "$chunk/bin/Rscript" <<'RSCRIPT'
#!/bin/bash
set -euo pipefail
gene=""; weights=""
while (($#)); do case $1 in --gene) gene=$2; shift 2;; --weights_dir) weights=$2; shift 2;; *) shift;; esac; done
[[ -n "$gene" ]] || exit 0
printf '%s\n' "$gene" >> "$CHUNK_LOG"; mkdir -p "$weights"; printf 'weight\n' > "$weights/$gene.wgt.RDat"
RSCRIPT
chmod +x "$chunk/bin/"*
export CHUNK_LOG="$chunk/genes.log"
export ZCAT_LOG="$chunk/zcat.log"
common=(--tissue T --chunk_size 2 --output_dir "$chunk/out" --egret_subdir custom --folds 0 \
    --plink_path "$chunk/bin/plink" --gemma_path gemma --expression "$chunk/out/expression_files/expression.gz" \
    --gene_info gene_info --conda_env fake --cleanup_egret_plink TRUE)
PATH="$chunk/bin:$PATH" "$scripts_dir/6.3_train_EGRET_chunk.sh" "${common[@]}" --start_index 0
PATH="$chunk/bin:$PATH" "$scripts_dir/6.3_train_EGRET_chunk.sh" "${common[@]}" --start_index 2
[[ $(sort "$chunk/genes.log" | uniq | tr '\n' ' ') == 'g1 g2 g3 g4 ' ]] || fail "chunk assignment was incomplete"
[[ $(wc -l < "$chunk/genes.log") -eq 4 ]] || fail "chunk assignment contained duplicates"
[[ $(wc -l < "$chunk/zcat.log") -eq 2 ]] || fail "expression was decompressed more than once per chunk"

# Source accounting distinguishes missing, empty, unmatched, and genotype-matched inputs.
src="$test_root/source"; mkdir -p "$src/genotype_files" "$src/MatrixeQTL/T/fold_0/src" "$src/bin"
printf '1\trs1\t0\t100\tA\tG\n' > "$src/genotype_files/geno.bim"
printf '#chrom\tchromStart\tchromEnd\tgeneId\nchr1\t10\t20\tg1\nchr1\t20\t30\tg2\nchr1\t30\t40\tg3\nchr1\t40\t50\tg4\n' > "$src/gene_info.tsv"
printf 'gene_id\ts1\ng1\t1\ng2\t2\ng3\t3\ng4\t4\n' | gzip -c > "$src/expression.gz"
: > "$src/MatrixeQTL/T/fold_0/src/g2.txt"
printf 'not_in_bim\n' > "$src/MatrixeQTL/T/fold_0/src/g3.txt"
printf 'rs1\n' > "$src/MatrixeQTL/T/fold_0/src/g4.txt"
cat > "$src/bin/plink" <<'SOURCE_PLINK'
#!/bin/bash
set -euo pipefail
out=""; excluded=0
while (($#)); do
    case $1 in
        --out) out=$2; shift 2;;
        --exclude) excluded=1; shift 3;;
        *) shift;;
    esac
done
mkdir -p "$(dirname "$out")"; printf 'bed\n' > "$out.bed"
if ((excluded)); then
    printf '1\trs2\t0\t200\tA\tG\n' > "$out.bim"
else
    printf '1\trs1\t0\t100\tA\tG\n' > "$out.bim"
fi
printf '0\ts1\t0\t0\t0\t-9\n' > "$out.fam"
SOURCE_PLINK
chmod +x "$src/bin/plink"
Rscript "$scripts_dir/6.2_make_EGRET_beds_by_fold.R" --tissue T \
    --output_dir "$src" --matrixeqtl_bed_dir src --gbat_bed_dir none --transpco_bed_dir none \
    --egret_subdir custom --exclude_crossmap FALSE --crossmap_dir none --genotype_prefix geno \
    --plink_path "$src/bin/plink" --expression "$src/expression.gz" --bed_models MatrixeQTL \
    --fold 0 --gene_info "$src/gene_info.tsv" > "$src/accounting.out"
assert_grep 'MatrixeQTL files=3 listed_snps=2 matched_snps=1 contributing_genes=1' "$src/accounting.out"
if Rscript "$scripts_dir/6.2_make_EGRET_beds_by_fold.R" --tissue T \
    --output_dir "$src" --matrixeqtl_bed_dir src --gbat_bed_dir none --transpco_bed_dir none \
    --egret_subdir custom --exclude_crossmap FALSE --crossmap_dir none --genotype_prefix geno \
    --plink_path "$src/bin/plink" --expression "$src/expression.gz" --bed_models GBAT \
    --fold 0 --gene_info "$src/gene_info.tsv" > "$src/empty.out" 2>&1; then
    fail "zero-match requested source unexpectedly passed"
fi
assert_grep 'zero SNPs after final PLINK/cross-map filtering: GBAT' "$src/empty.out"

mkdir -p "$src/cross_mapped/none"
printf '1\t100\t101\n' > "$src/cross_mapped/none/g4.bed"
if Rscript "$scripts_dir/6.2_make_EGRET_beds_by_fold.R" --tissue T \
    --output_dir "$src" --matrixeqtl_bed_dir src --gbat_bed_dir none --transpco_bed_dir none \
    --egret_subdir custom --exclude_crossmap TRUE --crossmap_dir none --genotype_prefix geno \
    --plink_path "$src/bin/plink" --expression "$src/expression.gz" --bed_models MatrixeQTL \
    --fold 0 --gene_info "$src/gene_info.tsv" > "$src/postfilter.out" 2>&1; then
    fail "source removed by final PLINK/cross-map filtering unexpectedly passed"
fi
assert_grep 'MatrixeQTL files=3 listed_snps=2 matched_snps=0 contributing_genes=0' "$src/postfilter.out"

# Folds may have different module counts. Completeness is checked against the modules that
# actually exist in each fold, while expected_module_count describes their union.
pco="$test_root/pco"; mkdir -p "$pco/transPCO/T/fold_0/modules" \
    "$pco/transPCO/T/fold_1/modules" "$pco/transPCO/T/fold_0/PCO_association_results" \
    "$pco/transPCO/T/fold_1/PCO_association_results"
printf 'g1\ng3\n' > "$pco/transPCO/T/fold_0/modules/module_1.txt"
printf 'g2\ng4\n' > "$pco/transPCO/T/fold_0/modules/module_2.txt"
printf 'g1\ng3\n' > "$pco/transPCO/T/fold_1/modules/module_1.txt"
printf 'result\n' > "$pco/transPCO/T/fold_0/PCO_association_results/module_1_chr_1.txt.gz"
printf 'result\n' > "$pco/transPCO/T/fold_1/PCO_association_results/module_1_chr_1.txt.gz"
if PATH="$test_root/bin:$PATH" "$scripts_dir/5.7_run_results_analysis_job.sh" \
    --tissue T --fdr 0.1 --output_dir "$pco" --folds 1 --expected_module_count 2 \
    --genotype_prefix geno --conda_env fake > "$pco/missing.out" 2>&1; then
    fail "incomplete per-fold transPCO modules unexpectedly passed"
fi
assert_grep 'missing PCO results for fold_0:module_2' "$pco/missing.out"
printf 'result\n' > "$pco/transPCO/T/fold_0/PCO_association_results/module_2_chr_1.txt.gz"
PATH="$test_root/bin:$PATH" "$scripts_dir/5.7_run_results_analysis_job.sh" \
    --tissue T --fdr 0.1 --output_dir "$pco" --folds 1 --expected_module_count 2 \
    --genotype_prefix geno --conda_env fake

# External PLINK failures are fatal and cannot leave a reduced, apparently valid gene set.
plink_fail="$test_root/plink_fail"
printf '#!/bin/bash\nexit 7\n' > "$plink_fail"; chmod +x "$plink_fail"
if Rscript "$scripts_dir/1.1_setup_genotypes.R" --genotype_prefix input --out "$test_root/pruned" \
    --plink_path "$plink_fail" > "$test_root/plink_setup_fail.out" 2>&1; then
    fail "genotype setup ignored a PLINK failure"
fi
assert_grep 'PLINK LD pruning failed with exit status 7' "$test_root/plink_setup_fail.out"

mkdir -p "$test_root/cis_fail/expression_files"
printf 'gene_id\ts1\nG1\t1\n' | gzip -c > "$test_root/cis_fail/expression_files/T_expression.txt.gz"
printf '#chrom\tchromStart\tchromEnd\tgeneId\nchr1\t100\t200\tG1\n' > "$test_root/cis_fail/gene_info.txt"
if Rscript "$scripts_dir/2.1_setup_cis_windows.R" --gene_info "$test_root/cis_fail/gene_info.txt" \
    --tissue T --plink_path "$plink_fail" --output_dir "$test_root/cis_fail" \
    --genotype_prefix geno > "$test_root/plink_cis_fail.out" 2>&1; then
    fail "cis setup ignored a PLINK failure"
fi
assert_grep 'PLINK cis-window extraction failed for G1 with exit status 7' "$test_root/plink_cis_fail.out"

mkdir -p "$test_root/fold_fail"
printf 'gene_id\ts1\ts2\nG1\t1\t2\n' > "$test_root/fold_fail/expression.txt"
printf 's1\ns2\n' > "$test_root/fold_fail/individuals.txt"
if Rscript "$scripts_dir/1.3_setup_folds.R" --expression "$test_root/fold_fail/expression.txt" \
    --individuals "$test_root/fold_fail/individuals.txt" --tissue T --folds 2 \
    --plink_path "$plink_fail" --genotype_prefix input --output_dir "$test_root/fold_fail" \
    > "$test_root/plink_fold_fail.out" 2>&1; then
    fail "fold setup ignored a PLINK failure"
fi
assert_grep 'PLINK fold extraction failed for fold_0 with exit status 7' "$test_root/plink_fold_fail.out"

mkdir -p "$test_root/chr_fail/genotype_files" "$test_root/chr_fail/fold_0_info/T"
printf '1\trs1\t0\t100\tA\tG\n' > "$test_root/chr_fail/genotype_files/geno.bim"
printf 's1\n' > "$test_root/chr_fail/fold_0_info/T/train_individuals.txt"
if Rscript "$scripts_dir/5.2_initialize_chr_genotypes.R" --output_dir "$test_root/chr_fail" \
    --tissue T --folds 0 --genotype_prefix geno --plink_path "$plink_fail" \
    > "$test_root/plink_chr_fail.out" 2>&1; then
    fail "transPCO chromosome setup ignored a PLINK failure"
fi
assert_grep 'PLINK chromosome extraction failed for fold_0 chr1 with exit status 7' \
    "$test_root/plink_chr_fail.out"

# The source fixture already contains successful outputs, so this also proves that a
# failed retry cannot be hidden by stale per-gene PLINK files.
if Rscript "$scripts_dir/6.2_make_EGRET_beds_by_fold.R" --tissue T \
    --output_dir "$src" --matrixeqtl_bed_dir src --gbat_bed_dir none --transpco_bed_dir none \
    --egret_subdir custom --exclude_crossmap FALSE --crossmap_dir none --genotype_prefix geno \
    --plink_path "$plink_fail" --expression "$src/expression.gz" --bed_models Cis \
    --fold 0 --gene_info "$src/gene_info.tsv" > "$test_root/plink_egret_fail.out" 2>&1; then
    fail "EGRET BED setup reused stale output after a PLINK failure"
fi
assert_grep 'PLINK EGRET extraction failed for g1 fold_0 with exit status 7' \
    "$test_root/plink_egret_fail.out"

echo "Round 3 workflow smoke tests passed"
