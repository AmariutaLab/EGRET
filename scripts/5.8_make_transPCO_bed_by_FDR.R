# --- 5.8_make_transPCO_bed_by_FDR.R ----------------------------------------
# FDR-filters the PCO p-values across all modules and writes one BED file of
# significant trans SNPs per module gene.
#
#   in : transPCO/{tissue}/fold_{i}/PCO_association_results/
#   out: transPCO/{tissue}/fold_{i}/{transpco_bed_dir}/{gene}.txt; significant_results_*.txt
#   run: 5.7_run_results_analysis_job.sh
# --------------------------------------------------------------------------
library(data.table)
library(optparse)

# Define command-line options
option_list = list(
  make_option("--tissue", type = "character", default = NA,
              help = "Tissue name to process", metavar = "character"),
  make_option("--FDR", type = "numeric", default = NA,
              help = "false discovery rate", metavar = "numeric"),
  make_option("--output_dir", type = "character", default = NA,
              help = "Base output directory", metavar = "character"),
  make_option("--bed_dir", type = "character", default = NA,
              help = "Subdirectory name for bed output files", metavar = "character"),
  make_option("--folds", type = "integer", default = NA,
              help = "Number of cross-validation folds", metavar = "integer"),
  make_option("--genotype_prefix", type = "character", default = NA,
              help = "Prefix for genotype files", metavar = "character"),
  make_option("--PCO_association_dir",type = 'character',default = NA,
	      help = "place where association results are stored"),
  make_option("--module_dir", type = "character", default = NA,
              help = "dir containing modules", metavar = "character")
)

# Parse command-line arguments
opt = parse_args(OptionParser(option_list = option_list))


required = c("tissue", "FDR", "output_dir", "bed_dir", "folds", "genotype_prefix", "PCO_association_dir", "module_dir")
missing = required[!vapply(required, function(o) length(opt[[o]]) == 1L && !is.na(opt[[o]]), logical(1))]
if (length(missing)) stop("missing required argument(s): --", paste(missing, collapse = ", --"), call. = FALSE)

# Extract arguments
tissue = opt$tissue
FDR = opt$FDR
output_dir = opt$output_dir
PCO_association_dir = opt$PCO_association_dir

# Initialize results data.table
all_snp_info = fread(paste0(output_dir, "/genotype_files/", opt$genotype_prefix, ".bim"), header = F)
# Same .bim-derived chromosome count as 9.1/9.3/5.4_MatrixeQTL_for_transPCO.R, replacing a hardcoded 1:22.
number_chr = max(all_snp_info$V1)

for (fold in 0:opt$folds) {
  empty_sig_results = data.table(module = character(), rsIDs = character(), pval = numeric(), FDR = numeric())
  sig_result_parts = list()
  # Determine number of modules
  module_dir = paste0(output_dir, "/transPCO/", tissue, "/fold_", fold, "/",opt$module_dir,"/")
  modules = list.files(module_dir)

  PCO_association_path = paste0(output_dir, "/transPCO/", tissue, "/fold_", fold, "/",opt$PCO_association_dir,"/")

  # Cleared, not just created: the per-gene writes below append so a gene in several
  # modules accumulates SNPs from each, which would double on a re-run.
  bed_dir = paste0(output_dir, "/transPCO/", tissue, "/fold_", fold, "/",opt$bed_dir,"/")
  if (dir.exists(bed_dir)) {
    unlink(bed_dir, recursive = TRUE)
  }
  dir.create(bed_dir, recursive = TRUE)
  largest_sig_pval = c()
  for (module in modules) {

    module = strsplit(module,"\\.")[[1]][1]
    module_sig_parts = list()
    print(module) 
    for (chr in 1:number_chr) {
      result_file = paste0(PCO_association_path, module, "_chr_", chr, ".txt.gz")
      # Skip if result file doesn't exist
      print(result_file)
      if (!file.exists(result_file)) {
        next
      }
		
      # Read results and filter by p-value
      results = fread(result_file, header = TRUE)
      results$FDR = p.adjust(results$pval,n = nrow(all_snp_info), method = 'fdr')
      sig_rows = results[results$FDR < opt$FDR,]
      largest_sig_pval = c(largest_sig_pval,sig_rows$pval[nrow(sig_rows)])
      if (nrow(sig_rows) > 0) {
        module_sig_parts[[length(module_sig_parts) + 1L]] = data.table(module = module, sig_rows)
      }
    }
    module_sig_results = if (length(module_sig_parts)) rbindlist(module_sig_parts) else copy(empty_sig_results)
    if (nrow(module_sig_results) > 0) {
      sig_result_parts[[length(sig_result_parts) + 1L]] = module_sig_results
    }
    # Do NOT rebuild sig_results here: accumulate into sig_result_parts and bind once
    # after the loop. Binding per module to print restores the quadratic cost.
    print(sprintf("module %s: %d significant SNP-gene rows", module, nrow(module_sig_results)))
    # Read genes in the current module
    genes_in_module = fread(paste0(module_dir, module,".txt"), header = FALSE)
    genes_in_module = unlist(genes_in_module$V1)

    # Identify significant SNPs for each gene in this module
    for (gene in genes_in_module) {
      sig_snps_per_gene = as.matrix(module_sig_results[, c('rsIDs','pval')])
      if (length(sig_snps_per_gene) > 0) {
        fwrite(sig_snps_per_gene, 
               paste0(bed_dir, gene, ".txt"), 
               col.names = FALSE, row.names = FALSE, sep = '\t', quote = FALSE,append = TRUE)
      }
    }
  }
  sig_results = if (length(sig_result_parts)) rbindlist(sig_result_parts) else copy(empty_sig_results)
  print(mean(largest_sig_pval,na.rm = T))
  fwrite(sig_results,paste0(output_dir, "/transPCO/",tissue,"/fold_",fold,"/significant_results_",opt$module_dir,"_FDR_",FDR,".txt"), col.names = T, row.names = F, sep = '\t', quote = F)
}


