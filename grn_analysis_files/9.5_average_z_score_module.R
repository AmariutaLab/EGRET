# --- 9.5_average_z_score_module.R ------------------------------------------
# GRN step 4. Compares mean TWAS Z per co-regulation module between the EGRET and cis
# arms. Needs both arms' TWAS results, so 9.1 skips it when no GWAS directory is given.
#
#   in : grn_analysis/coregulation_grns/; TWAS_association_results/{tissue}/{EGRET,cis}/
#   out: a per-module Z-score summary table
#   run: 9.1_run_GRN_analysis.sh, only when --gwas_sumstat_dir is supplied
# --------------------------------------------------------------------------
library(data.table)
library(optparse)

option_list = list(
  make_option("--tissues", action = "store", default = NA, type = 'character',
              help = "comma-separated tissue names, or 'all' to sweep every tissue with a results_sumstats directory"),
  make_option("--output_dir", action = "store", default = NA, type = 'character',
              help = "Path to output directory"),
  make_option("--egret_subdir", action = "store", default = NA, type = 'character',
              help = "the name of the analysis"),
  make_option("--gene_info", action = "store", default = NA, type = 'character',
              help = "path to gene info file"),
  make_option("--gwas_sumstat_dir", action = "store", default = NA, type = 'character',
              help = "directory of GWAS *.sumstats.gz; the trait list is read from it"),
  make_option("--gwas_traits", action = "store", default = "", type = 'character',
              help = "optional comma-separated subset of trait file names, overriding --gwas_sumstat_dir's listing"),
  make_option("--pval_threshold", action = "store", default = 0.01, type = 'numeric',
              help = "trans-component significance threshold [default: %default]"),
  make_option("--out", action = "store", default = "", type = 'character',
              help = "output table [default: {output_dir}/grn_analysis/module_z_score_summary.txt]")
  )

opt = parse_args(OptionParser(option_list = option_list))

output_dir = opt$output_dir
project = opt$egret_subdir

for (required in c("tissues", "output_dir", "egret_subdir", "gene_info", "gwas_sumstat_dir")) {
  if (is.na(opt[[required]])) stop("missing required option: --", required)
}

# 'all' discovers tissues from what the analysis stage actually produced, rather than from a
# hardcoded list that silently drifts out of step with the run.
if (identical(opt$tissues, "all")) {
  tissues = list.files(paste0(output_dir, "/results_sumstats/"))
  if (length(tissues) == 0) stop("no tissues found under ", output_dir, "/results_sumstats/")
} else {
  tissues = trimws(strsplit(opt$tissues, ",")[[1]])
}

# The trait list is read from what is actually in gwas_sumstat_dir; --gwas_traits narrows
# it to a subset, mirroring the workflow's own gwas_traits variable.
if (nzchar(opt$gwas_traits)) {
  general_sumstats = trimws(strsplit(opt$gwas_traits, ",")[[1]])
} else {
  general_sumstats = list.files(opt$gwas_sumstat_dir, pattern = "\\.sumstats\\.gz$")
  if (length(general_sumstats) == 0) stop("no *.sumstats.gz found in ", opt$gwas_sumstat_dir)
}

out_file = opt$out
if (!nzchar(out_file)) out_file = paste0(output_dir, "/grn_analysis/module_z_score_summary.txt")
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

gene_info = fread(opt$gene_info, header = TRUE)
all_results = data.frame(
  tissue = character(),
  trait = character(),
  module = character(),
  avg_z_EGRET = numeric(),
  avg_z_cis = numeric(),
  num_trans_regulated = numeric(),
  prop_trans_regulated = numeric(),
  module_size = numeric(),
  t_test_p = numeric(),
  wilcoxon_p = numeric()
)

for (tissue in tissues) {
  # Guarded: with --tissues all, a tissue whose analysis stage has not run would otherwise
  # abort the whole sweep. The trait loop below skips missing files the same way.
  sumstats_path = paste0(output_dir, "/results_sumstats/", tissue, "/", project, ".txt")
  if (!file.exists(sumstats_path)) {
    message("Skipping ", tissue, " - no analysis summary at ", sumstats_path)
    next
  }
  model_performance = fread(sumstats_path, header = T)

  module_dir = paste0(output_dir, "/grn_analysis/coregulation_grns/", tissue, "/")
  if (!dir.exists(module_dir)) {
    message("Skipping ", tissue, " - no coregulation GRNs at ", module_dir)
    next
  }
  modules = list.files(
    module_dir,
    pattern = "^coregulation_cluster_[0-9]+\\.txt$",
    full.names = TRUE
  )
  if (length(modules) == 0) {
    message("Skipping ", tissue, " - no coregulation cluster tables at ", module_dir)
    next
  }

  for (trait in general_sumstats) {
    trait_file = paste0(substr(trait, 1, nchar(trait) - 12), ".dat")
	  print(trait_file)
    TWAS_path = paste0(output_dir, "/TWAS_association_results/", tissue, "/", project, "/", trait_file)
    cis_path = paste0(output_dir, "/TWAS_association_results/", tissue, "/cis/", trait_file)

    if (!file.exists(TWAS_path) || !file.exists(cis_path)) {
      message("Skipping ", trait, " - TWAS results not present")
      next
    }
    TWAS_results = fread(TWAS_path, header = TRUE)
    cis_TWAS_results = fread(cis_path, header = TRUE)

    TWAS_results[is.na(TWAS_results$TWAS.Z), "TWAS.Z"] = 0
    cis_TWAS_results[is.na(cis_TWAS_results$TWAS.Z), "TWAS.Z"] = 0

    for (module_path in modules) {
      module = basename(module_path)
      genes_in_module = fread(module_path, header = FALSE)
      gene_ids = gene_info$name[gene_info$geneId %in% genes_in_module$V1]

      TWAS_module_subset = merge(TWAS_results, genes_in_module, by.x = 'ID', by.y = 'V1', all.y = TRUE)
      TWAS_module_subset[is.na(TWAS_module_subset$TWAS.Z), "TWAS.Z"] = 0

      cis_module_subset = merge(cis_TWAS_results, genes_in_module, by.x = 'ID', by.y = 'V1', all.y = TRUE)
      cis_module_subset[is.na(cis_module_subset$TWAS.Z), "TWAS.Z"] = 0

      avg_z_EGRET = mean(abs(TWAS_module_subset$TWAS.Z), na.rm = TRUE)
      avg_z_cis = mean(abs(cis_module_subset$TWAS.Z), na.rm = TRUE)

      z_egret <- abs(TWAS_module_subset$TWAS.Z)
      z_cis   <- abs(cis_module_subset$TWAS.Z)

      # Ensure pairing is valid
      stopifnot(length(z_egret) == length(z_cis))

      d = z_egret - z_cis

      if (length(d) < 2 || sd(d) == 0) {
          t_test_results = list(
            statistic = NA,
            p.value = 1,
            estimate = mean(d),
            method = "paired t-test (skipped: zero variance)"
          )
        } else {
          t_test_results = t.test(
            z_egret,
            z_cis,
            paired = TRUE,
            alternative = "greater"
          )
        }

      if (length(unique(d)) > 1) {
        wilcoxon_test = wilcox.test(
          z_egret,
          z_cis,
          paired = TRUE,
          alternative = "greater",
          exact = FALSE
        )
        wilcoxon_p = wilcoxon_test$p.value
      } else {
        wilcoxon_p = 1
      }

      module_gene_performance = model_performance[gene %in% unlist(genes_in_module),]
      module_gene_performance$p_trans_part[is.na(module_gene_performance$p_trans_part)] = 1

      num_trans_regulated = sum(module_gene_performance$p_trans_part < opt$pval_threshold)
      prop_trans_regulated = num_trans_regulated/nrow(genes_in_module)


      all_results = rbind(all_results,
        data.frame(
          tissue = tissue,
		      trait = trait_file,
          module = module,
          avg_z_EGRET = avg_z_EGRET,
          avg_z_cis = avg_z_cis,
          num_trans_regulated = num_trans_regulated,
          prop_trans_regulated = prop_trans_regulated,
          module_size = length(gene_ids),
          t_test_p = t_test_results$p.value,
          wilcoxon_p = wilcoxon_p
        )
      )

    }
  }
}

fwrite(all_results, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote ", nrow(all_results), " rows to ", out_file)
