#!/bin/bash
# --- validate_egret.sh -----------------------------------------------------
# End-to-end checker. `prepare` builds a reduced input set; `check` asserts the run
# produced what it should; `usage` reports inode and disk consumption.
# Every assertion is gated so it reports SKIP rather than a false PASS when the stage
# that fills a tree never ran.
#
#   in : an output_dir produced by a validation run
#   out: a PASS/FAIL report on stdout; exit 0 only when FAIL=0
#   run: ./validate_egret.sh {prepare|check|usage}
# --------------------------------------------------------------------------
# End-to-end validation harness for the EGRET pipeline.
#
# The run itself is the normal orchestrator pointed at a small config; this script only
# builds its inputs and asserts the outcomes, adding no orchestration of its own.
#
#   ./validate_egret.sh prepare   # subset genotypes + expression into validation/inputs
#   sbatch run_EGRET_workflow.sh --config config_validation.sh
#   ./validate_egret.sh check     # assert the run produced what it should
#   ./validate_egret.sh usage     # report inode/byte footprint
set -euo pipefail
scripts_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); source "${scripts_dir}/workflow_lib.sh"

config_file="${scripts_dir}/config_validation.sh"
[[ -r "$config_file" ]] || cli_die "missing $config_file"
# shellcheck disable=SC1090
source "$config_file"
source_config="${scripts_dir}/config.sh"

n_genes=${VALIDATION_GENES:-1800}
n_snps=${VALIDATION_SNPS:-80000}

prepare() {
    activate_workflow_env "$conda_env"
    # grep exits 1 on a missing key, which under `set -euo pipefail` kills prepare with no
    # output. Read each key explicitly and name the one that is missing.
    local src_geno src_expr src_gene_info
    cfg_value() {
        local key=$1 val
        val=$(grep -m1 "^${key}=" "$source_config" | cut -d'"' -f2) || true
        [[ -n "$val" ]] || cli_die "could not read '${key}' from ${source_config}"
        printf '%s\n' "$val"
    }
    src_geno=$(cfg_value genotype_prefix)
    src_expr=$(cfg_value expression)
    src_gene_info=$(cfg_value gene_info)
    local indir="${scripts_dir}/validation/inputs"
    mkdir -p "$indir"

    echo "== building eligible gene universe =="
    zcat "$src_gene_info" | awk -F'\t' 'NR>1 && ($8=="protein_coding"||$8=="lincRNA"||$8=="antisense"){print $7}' \
        | sort -u > "${indir}/eligible.txt"
    zcat "$src_expr" | awk -F'\t' 'NR>1{print $4}' | sort -u > "${indir}/in_expr.txt"
    comm -12 "${indir}/eligible.txt" "${indir}/in_expr.txt" > "${indir}/both.txt"
    local total; total=$(wc -l < "${indir}/both.txt")
    (( total > 0 )) || cli_die "no genes are eligible in both gene_info and the expression file"

    # Even stride, not head -n, so every chromosome stays represented: transPCO's trans
    # signal only exists between genes and SNPs on *different* chromosomes.
    local stride=$(( total / n_genes )); (( stride >= 1 )) || stride=1
    awk -v s="$stride" 'NR % s == 1' "${indir}/both.txt" | head -n "$n_genes" > "${indir}/genes.txt"
    echo "   selected $(wc -l < "${indir}/genes.txt") of $total eligible genes (stride $stride)"

    echo "== subsetting expression =="
    zcat "$src_expr" | awk -F'\t' -v OFS='\t' -v gf="${indir}/genes.txt" '
        BEGIN { while ((getline g < gf) > 0) keep[g]=1 }
        NR==1 { print; next }
        $4 in keep { print }
    ' | gzip > "${indir}/validation_expression.bed.gz"
    echo "   rows: $(zcat "${indir}/validation_expression.bed.gz" | wc -l) (incl. header)"

    echo "== thinning genotypes to ~${n_snps} SNPs =="
    "$plink_path" --bfile "$src_geno" --thin-count "$n_snps" --seed 1234 \
        --make-bed --out "${indir}/validation_genotypes" >/dev/null 2>&1
    echo "   SNPs: $(wc -l < "${indir}/validation_genotypes.bim") across $(cut -f1 "${indir}/validation_genotypes.bim" | sort -un | wc -l) chromosomes"

    # Every chromosome 1..max must be non-empty: transPCO stages 5.2, 5.4, 5.6, and 5.8
    # all loop across that full range,
    # so a gap would make plink emit an empty per-chr fileset mid-run.
    local maxchr missing
    maxchr=$(cut -f1 "${indir}/validation_genotypes.bim" | sort -un | tail -1)
    missing=$(comm -23 <(seq 1 "$maxchr") <(cut -f1 "${indir}/validation_genotypes.bim" | sort -un) | tr '\n' ' ')
    [[ -z "${missing// }" ]] || cli_die "chromosomes with no SNPs after thinning: $missing"
    echo "   contiguous chromosomes 1..${maxchr}: OK"
    rm -f "${indir}/eligible.txt" "${indir}/in_expr.txt" "${indir}/both.txt"
    echo "prepare complete -> $indir"
}

check() {
    local fail=0 pass=0
    ok(){ echo "  PASS  $1"; pass=$((pass+1)); }
    no(){ echo "  FAIL  $1"; fail=$((fail+1)); }
    # count_files (workflow_lib.sh) instead of a raw find pipeline: under set -euo
    # pipefail, find on a missing directory exits 1 and silently kills this checker --
    # and a *passing* run is precisely when these directories have been removed by the
    # cleanup hooks. The checker must survive its own success.
    t(){ count_files "$1" -maxdepth 1 -type f -name "$2"; }

    echo "== setup =="
    [[ -s "${output_dir}/expression_files/${tissue}_expression_regressed.txt.gz" ]] \
        && ok "regressed expression written" || no "regressed expression missing"
    # The covariate regression must not produce an all-NA matrix.
    local expr_gz="${output_dir}/expression_files/${tissue}_expression_regressed.txt.gz"
    if [[ -s "$expr_gz" ]]; then
        # Read the WHOLE file and require gzip to exit clean. Piping through
        # `grep -c ... || true` would report PASS on a file it could not read.
        local decoded nacount
        if decoded=$(zcat "$expr_gz" 2>/dev/null); then
            nacount=$(printf '%s\n' "$decoded" | tail -n +2 | grep -c 'NA' || true)
            (( nacount == 0 )) && ok "regressed expression has no NA values (qr.resid fix)" \
                               || no "regressed expression contains NA in $nacount row(s)"
        else
            no "regressed expression is unreadable/truncated: $expr_gz"
        fi
    fi

    echo "== trans-eQTL sources =="
    for ((f=0; f<=folds; f++)); do
        (( $(t "${output_dir}/MatrixeQTL/${tissue}/fold_${f}/${matrixeqtl_bed_dir}" '*.txt') > 0 )) \
            && ok "MatrixeQTL fold_${f} beds" || no "MatrixeQTL fold_${f} beds empty"
        (( $(t "${output_dir}/GBAT/${tissue}/fold_${f}/${gbat_bed_dir}" '*.txt') > 0 )) \
            && ok "GBAT fold_${f} beds" || no "GBAT fold_${f} beds empty"
        # The producer and the consumer must agree on this directory name -- when they did
        # not, transPCO contributed zero SNPs to every model and nothing else complained.
        (( $(t "${output_dir}/transPCO/${tissue}/fold_${f}/${transpco_bed_dir}" '*.txt') > 0 )) \
            && ok "transPCO fold_${f} ${transpco_bed_dir}" \
            || no "transPCO fold_${f} ${transpco_bed_dir} EMPTY -- producer/consumer name mismatch?"
    done

    echo "== models and analysis =="
    (( $(t "${output_dir}/FUSION/${tissue}/cis" '*.wgt.RDat') > 0 )) && ok "cis weights" || no "cis weights missing"
    (( $(t "${output_dir}/EGRET_models/${tissue}/${egret_subdir}" '*.wgt.RDat') > 0 )) \
        && ok "EGRET weights" || no "EGRET weights missing"
    for s in "${egret_subdir}" cis; do
        [[ -s "${output_dir}/results_sumstats/${tissue}/${s}.txt" ]] && ok "sumstats ${s}.txt" || no "sumstats ${s}.txt missing"
    done
    for s in "${egret_subdir}" cis; do
        [[ -s "${output_dir}/pos_files/${tissue}/${s}.pos" ]] && ok "pos file ${s}.pos" || no "pos file ${s}.pos missing"
    done
    (( $(t "${output_dir}/TWAS_association_results/${tissue}/${egret_subdir}" '*.dat') > 0 )) \
        && ok "EGRET TWAS results" || no "EGRET TWAS results missing"
    (( $(t "${output_dir}/TWAS_association_results/${tissue}/cis" '*.dat') > 0 )) \
        && ok "cis TWAS results" || no "cis TWAS results missing"

    echo "== cleanup hooks =="
    # An empty tree only proves a hook fired if the stage that fills it actually ran;
    # without this gate every hook reports a false PASS on a run that stopped early.
    local models_ran=0
    (( $(t "${output_dir}/EGRET_models/${tissue}/${egret_subdir}" '*.wgt.RDat') > 0 )) && models_ran=1
    hook(){
        local label="$1" n="$2"
        if (( models_ran == 0 )); then echo "  SKIP  $label (model stage did not run; nothing to prove)"; return; fi
        (( n == 0 )) && ok "$label cleaned" || no "$label still holds $n files"
    }
    hook "working/"        "$(count_files "${output_dir}/working/${tissue}/${egret_subdir}" -type f)"
    hook "plink_results/"  "$(count_files "${output_dir}/plink_results/${tissue}/${egret_subdir}" -type f)"
    hook "z_matrices/"     "$(count_files "${output_dir}/z_matrices/${tissue}/${egret_subdir}" -type f)"
    hook "bed_files/"      "$(count_files "${output_dir}/bed_files/${tissue}/${egret_subdir}" -type f)"
    hook "MatrixeQTL raw"  "$(count_files "${output_dir}/MatrixeQTL/${tissue}" -maxdepth 1 -type f -name 'association_results_fold_*_*.txt')"
    hook "fold train_genotypes.txt" "$(count_files "${output_dir}" -maxdepth 3 -type f -name 'train_genotypes.txt')"
    hook "fold train_expression.txt" "$(count_files "${output_dir}" -maxdepth 3 -type f -name 'train_expression.txt')"
    hook "transPCO genotype_files" "$(count_files "${output_dir}/transPCO/${tissue}" -path '*/genotype_files/*' -type f)"
    # The per-gene hook in 2.2 is success-gated, so failed genes keep their scratch until
    # 7.0's terminal sweep. Gated on models_ran, since 7.0 is what fires that sweep.
    hook "working/cis"     "$(count_files "${output_dir}/working/${tissue}/cis" -type f)"
    hook "tmp/cis"         "$(count_files "${output_dir}/tmp/${tissue}/cis" -type f)"

    echo "== GRN (opt-in stage) =="
    # grn is absent from the default `stages` list, so an empty tree here means "not run",
    # not "broken". Report SKIP in that case rather than failing a run that never asked for it.
    if [[ -d "${output_dir}/grn_analysis" ]]; then
        (( $(t "${output_dir}/grn_analysis/trans_connections" '*.txt') > 0 )) \
            && ok "GRN trans_connections" || no "GRN trans_connections empty"
        (( $(t "${output_dir}/grn_analysis/trans_connections_by_eqtl" '*.txt') > 0 )) \
            && ok "GRN trans_connections_by_eqtl" || no "GRN trans_connections_by_eqtl empty"
        (( $(t "${output_dir}/grn_analysis/coregulation_grns/${tissue}" '*') > 0 )) \
            && ok "GRN coregulation_grns" || no "GRN coregulation_grns empty"
    else
        echo "  SKIP  GRN stage (grn_analysis/ absent; stage not enabled)"
    fi

    echo; echo "PASS=$pass FAIL=$fail"
    (( fail == 0 ))
}

usage_report() {
    echo "inodes: $(count_files "$output_dir" -type f)"
    du -sh --apparent-size "$output_dir" 2>/dev/null | awk '{print "apparent: "$1}'
    du -sh "$output_dir" 2>/dev/null | awk '{print "allocated: "$1}'
}

case "${1:-}" in
    prepare) prepare ;;
    check)   check ;;
    usage)   usage_report ;;
    *) echo "usage: $0 {prepare|check|usage}" >&2; exit 2 ;;
esac
