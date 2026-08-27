#!/bin/bash
# --- test_grn_workflow.sh --------------------------------------------------
# Asserts the GRN dispatcher and the four chain scripts accept their flags and run
# end to end on synthetic inputs.
#
#   in : nothing -- self-contained
#   out: PASS/FAIL on stdout
#   run: bash scripts/tests/test_grn_workflow.sh
# --------------------------------------------------------------------------
set -euo pipefail

scripts_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
repo_dir=$(cd "${scripts_dir}/.." && pwd)
python_bin=${PYTHON_BIN:-python}
rscript_bin=${RSCRIPT_BIN:-Rscript}
test_root=$(mktemp -d)
trap 'rm -rf "$test_root"' EXIT

fail() { echo "FAIL: $*" >&2; exit 1; }

# The dispatcher must depend on the actual per-trait jobs, not merely on the two TWAS
# dispatchers that wrote these files.
mkdir -p "$test_root/bin" "$test_root/logs"
printf '101:102\n' > "$test_root/twas.jids"
printf '201:202\n' > "$test_root/cis_twas.jids"
cat > "$test_root/bin/sbatch" <<'STUB'
#!/bin/bash
printf '%s\n' "$*" >> "$SBATCH_LOG"
echo 303
STUB
chmod +x "$test_root/bin/sbatch"
SBATCH_LOG="$test_root/sbatch.log" PATH="$test_root/bin:$PATH" \
    "$scripts_dir/9.0_submit_GRN_analysis.sh" \
    --tissue T --output_dir "$test_root/out" --egret_subdir EGRET \
    --gene_info "$test_root/gene_info.txt" --eqtl_file "$test_root/eqtl.txt" \
    --gwas_sumstat_dir "$test_root/gwas" --gwas_traits trait.sumstats.gz \
    --conda_env fake --slurm_account acct --slurm_partition part \
    --slurm_log_dir "$test_root/logs" --twas_job_id_file "$test_root/twas.jids" \
    --cis_twas_job_id_file "$test_root/cis_twas.jids" \
    --pval_threshold 0.2 \
    --job_id_file "$test_root/grn.jids"
grep -q -- '--dependency=afterok:101:102:201:202' "$test_root/sbatch.log" ||
    fail "GRN worker does not depend on all TWAS child jobs"
[[ $(cat "$test_root/grn.jids") == 303 ]] || fail "GRN worker ID was not recorded"
grep -q -- '--pval_threshold 0.2' "$test_root/sbatch.log" ||
    fail "GRN threshold was not forwarded to the worker"

# rsID conversion must disambiguate variants sharing a position by unordered allele pair.
# Reversed ref/alt is valid; two rsIDs with the same pair are ambiguous and remain unmatched.
mkdir -p "$test_root/converter"
printf 'variant_id\tgene_id\tpval_nominal\nchr1_100_C_T_b38\tG1\t0.01\nchr1_100_G_A_b38\tG2\t0.02\nchr1_200_T_C_b38\tG3\t0.03\nchr1_300_A_C_b38\tG4\t0.04\n' \
    > "$test_root/converter/eqtl.txt"
printf '1 rs_position_first 0 100 A G\n1 rs_allele_match 0 100 C T\n1 rs_swapped 0 200 C T\n1 rs_ambiguous_a 0 300 A C\n1 rs_ambiguous_b 0 300 C A\n' \
    > "$test_root/converter/genotypes.bim"
"$python_bin" "$repo_dir/grn_analysis_files/convert_variant_id_to_rsid.py" \
    --tissue T --eqtl_file "$test_root/converter/eqtl.txt" \
    --bim_file "$test_root/converter/genotypes.bim" --out_dir "$test_root/converter"
"$python_bin" -c 'import pandas as pd, sys; d=pd.read_csv(sys.argv[1], sep="\t"); assert list(d.columns)==["variant_id","gene_id","pval_nominal","rsid"]; assert d["rsid"].iloc[:3].tolist()==["rs_allele_match","rs_position_first","rs_swapped"]; assert pd.isna(d["rsid"].iloc[3])' \
    "$test_root/converter/T_with_rsid.txt"

# An empty trans-component selection must remain a valid, header-only network rather than
# crashing the following stages.
mkdir -p "$test_root/out/results_sumstats/T/EGRET" \
    "$test_root/out/results_sumstats/T" "$test_root/out/expression_files"
printf 'gene\tr2_gw\tr2_cis_part\tr2_trans_part\tp_trans_part\nG1\t0.1\t0.2\t0.0\t1\n' \
    > "$test_root/out/results_sumstats/T/EGRET.txt"
printf 'gene_id\tS1\nG1\t1\n' | gzip > "$test_root/out/expression_files/T_expression.txt.gz"
printf 'geneId\tname\t#chrom\tchromStart\tchromEnd\nG1\tG1NAME\tchr1\t100\t200\n' \
    > "$test_root/gene_info.txt"
printf 'rsid\tgene_id\tpval_nominal\nrs1\tR1\t0.001\n' > "$test_root/eqtl.txt"

"$python_bin" "$repo_dir/grn_analysis_files/9.2_identify_trans_regulators.py" \
    --tissue T --output_dir "$test_root/out" --project EGRET \
    --gene_info "$test_root/gene_info.txt"
"$python_bin" "$repo_dir/grn_analysis_files/9.3_assign_trans_regulators_with_eqtl.py" \
    --tissue T --output_dir "$test_root/out" --eqtl_file_path "$test_root/eqtl.txt"
"$python_bin" "$repo_dir/grn_analysis_files/9.4_determine_grns.py" \
    --tissue T --output_dir "$test_root/out"
[[ $(wc -l < "$test_root/out/grn_analysis/trans_connections/T_trans_connections.txt") == 1 ]] ||
    fail "empty GRN output is not header-only"

# Exactly MIN_CLUSTER_SIZE targets should form a module, and the summary must ignore the
# dendrogram PNG stored beside module tables.
cat > "$test_root/out/grn_analysis/trans_connections_by_eqtl/T_trans_connections_by_eqtl.txt" <<'EDGES'
target_gene	snp	nearest_gene_id	egene
G1	rs1	R1	R1
G2	rs2	R1	R1
G3	rs3	R1	R1
EDGES
"$python_bin" "$repo_dir/grn_analysis_files/9.4_determine_grns.py" \
    --tissue T --output_dir "$test_root/out"
module_file=$(find "$test_root/out/grn_analysis/coregulation_grns/T" \
    -maxdepth 1 -name 'coregulation_cluster_*.txt' -print -quit)
[[ -n "$module_file" ]] || fail "three-target network did not produce a module"
[[ -s "$test_root/out/grn_analysis/coregulation_grns/T/grn_dendrogram.png" ]] ||
    fail "dendrogram was not produced"

printf 'gene\tp_trans_part\nG1\t0.001\nG2\t0.5\nG3\t0.5\n' \
    > "$test_root/out/results_sumstats/T/EGRET.txt"
printf 'geneId\tname\t#chrom\tchromStart\tchromEnd\nG1\tG1NAME\tchr1\t100\t200\nG2\tG2NAME\tchr1\t300\t400\nG3\tG3NAME\tchr1\t500\t600\n' \
    > "$test_root/gene_info.txt"
mkdir -p "$test_root/out/TWAS_association_results/T/EGRET" \
    "$test_root/out/TWAS_association_results/T/cis" "$test_root/gwas"
printf 'ID\tTWAS.Z\nG1\t2\nG2\t1\nG3\t3\n' \
    > "$test_root/out/TWAS_association_results/T/EGRET/trait.dat"
printf 'ID\tTWAS.Z\nG1\t1\nG2\t1\nG3\t2\n' \
    > "$test_root/out/TWAS_association_results/T/cis/trait.dat"
"$rscript_bin" "$repo_dir/grn_analysis_files/9.5_average_z_score_module.R" \
    --tissues T --output_dir "$test_root/out" --egret_subdir EGRET \
    --gene_info "$test_root/gene_info.txt" --gwas_sumstat_dir "$test_root/gwas" \
    --gwas_traits trait.sumstats.gz --pval_threshold 0.6 --out "$test_root/module_summary.txt"
[[ $(wc -l < "$test_root/module_summary.txt") -eq 2 ]] ||
    fail "module summary did not contain exactly one result row"
awk -F '\t' 'NR == 2 { exit !($6 == 3 && $7 == 1) }' "$test_root/module_summary.txt" ||
    fail "module summary did not use the configured trans-component threshold"

echo "GRN workflow tests passed"
