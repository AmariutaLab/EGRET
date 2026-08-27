"""Add an `rsid` column to a GTEx eQTL summary-statistics file.

9.3_assign_trans_regulators_with_eqtl.py joins the GRN trans connections to GTEx eQTLs on an
`rsid` column, but GTEx v8 distributes only `variant_id` (chr_pos_ref_alt_b38). This script
produces the `{tissue}_with_rsid.txt` file that config.sh's `eqtl_file` must point at, by
matching each variant to the study genotypes' .bim on chromosome, position, and unordered
allele pair.

    python convert_variant_id_to_rsid.py \
        --tissue Whole_Blood \
        --eqtl_file GTEx_Analysis_v8_eQTL/Whole_Blood.v8.signif_variant_gene_pairs.txt.gz \
        --bim_file <output_dir>/genotype_files/GTEX_v8_genotypes_pruned.bim \
        --out_dir GTEx_Analysis_v8_eQTL

Ref/alt order may be reversed between GTEx and PLINK, so allele order is ignored. Alleles
are still required to match: position-only matching can silently assign the wrong rsid at a
multi-allelic locus. If more than one rsid has the same locus and allele pair, that key is
treated as ambiguous and is left unmatched.
"""
# --- convert_variant_id_to_rsid.py -----------------------------------------
# Utility, not a chain step. Adds an rsid column to a GTEx eQTL file by matching
# variant_id against the study .bim on chromosome, position and unordered alleles.
# Produces the file config.sh's eqtl_file must point at.
#
#   in : a GTEx v8 eQTL sumstats file; the study genotype .bim
#   out: {tissue}_with_rsid.txt
#   run: by hand, once per tissue, before the grn stage
# --------------------------------------------------------------------------

import argparse
import os

import pandas as pd

BIM_COLUMNS = ["chr", "rsid", "cm", "pos", "ref", "alt"]


def normalize_chromosome(values):
    return values.astype("string").str.replace(r"^chr", "", case=False, regex=True)


def unordered_allele_pair(first, second):
    first = first.astype("string").str.upper()
    second = second.astype("string").str.upper()
    low = first.where(first <= second, second)
    high = second.where(first <= second, first)
    return low + "/" + high


def convert_variant_id_to_rsid(tissue, eqtl_file_path, bim_file, out_dir):
    df = pd.read_csv(eqtl_file_path, sep="\t")
    if "variant_id" not in df.columns:
        raise ValueError(
            f"{eqtl_file_path} has no 'variant_id' column; columns are: "
            f"{', '.join(df.columns)}"
        )

    # PLINK1 normally writes tabs, while PLINK2 and conversion tools may use spaces.
    bim = pd.read_csv(bim_file, sep=r"\s+", header=None)
    if bim.shape[1] != len(BIM_COLUMNS):
        raise ValueError(
            f"{bim_file} has {bim.shape[1]} columns; a PLINK .bim must have 6"
        )
    bim.columns = BIM_COLUMNS

    # Parse the GTEx-style variant ID: chr1_64764_C_T_b38
    parts = df["variant_id"].astype("string").str.split("_", expand=True)
    if parts.shape[1] != 5 or parts.isna().to_numpy().any():
        raise ValueError(
            "variant_id values must all have GTEx form chr_pos_ref_alt_build "
            "(for example chr1_64764_C_T_b38)"
        )
    parts.columns = ["chr", "pos", "ref", "alt", "build"]
    parts["chr"] = normalize_chromosome(parts["chr"])
    parts["pos"] = pd.to_numeric(parts["pos"], errors="raise").astype("Int64")

    bim["chr"] = normalize_chromosome(bim["chr"])
    bim["pos"] = pd.to_numeric(bim["pos"], errors="raise").astype("Int64")
    parts["alleles"] = unordered_allele_pair(parts["ref"], parts["alt"])
    bim["alleles"] = unordered_allele_pair(bim["ref"], bim["alt"])
    parts["match_key"] = (
        parts["chr"] + ":" + parts["pos"].astype("string") + ":" + parts["alleles"]
    )
    bim["match_key"] = (
        bim["chr"] + ":" + bim["pos"].astype("string") + ":" + bim["alleles"]
    )

    rsids_per_key = bim.groupby("match_key", dropna=False)["rsid"].nunique(dropna=True)
    ambiguous_keys = set(rsids_per_key[rsids_per_key > 1].index)
    bim_unique = (
        bim[~bim["match_key"].isin(ambiguous_keys)]
        .drop_duplicates(subset="match_key")
        .set_index("match_key")["rsid"]
    )

    merged = df.copy()
    merged["rsid"] = parts["match_key"].map(bim_unique)

    matched = int(merged["rsid"].notna().sum())
    missing = int(merged["rsid"].isna().sum())
    print(
        f"{matched} variants matched to an rsid, {missing} unmatched; "
        f"{len(ambiguous_keys)} ambiguous BIM allele keys excluded"
    )
    if matched == 0:
        raise ValueError(
            "no variant matched the .bim on chromosome, position, and alleles. Check that "
            "the files use the same genome build and allele representation."
        )

    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"{tissue}_with_rsid.txt")
    merged.to_csv(out_path, sep="\t", index=False)
    print(f"Wrote {len(merged)} rows to {out_path}")
    return out_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tissue", required=True, help="tissue name, used for the output filename")
    parser.add_argument("--eqtl_file", required=True,
                        help="GTEx *.signif_variant_gene_pairs.txt.gz for this tissue")
    parser.add_argument("--bim_file", required=True,
                        help="PLINK .bim providing the rsid for each chromosome/position")
    parser.add_argument("--out_dir", required=True,
                        help="directory to write {tissue}_with_rsid.txt into")
    args = parser.parse_args()
    convert_variant_id_to_rsid(args.tissue, args.eqtl_file, args.bim_file, args.out_dir)


if __name__ == "__main__":
    main()
