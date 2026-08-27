library(data.table)
library(optparse)

# Define command-line options
option_list = list(
  make_option("--tissue", type = "character", default = "Whole_Blood",
              help = "Tissue name to process", metavar = "character"),
  make_option("--FDR", type = "numeric", default = 0.1,
              help = "false discovery rate", metavar = "numeric"),
  make_option("--output_dir", type = "character", default = ".",
              help = "Base output directory", metavar = "character"),
  make_option("--bed_dir", type = "character", default = "bed_files",
              help = "Subdirectory name for bed output files", metavar = "character"),
  make_option("--folds", type = "integer", default = 5,
              help = "Number of cross-validation folds", metavar = "integer"),
  make_option("--genotype_prefix", type = "character", default = "GTEX_v8_genotypes_pruned",
              help = "Prefix for genotype files", metavar = "character"),
  make_option("--PCO_association_dir",type = 'character',default = 'PCO_association_results',
	      help = "place where association results are stored"),
  make_option("--module_dir", type = "character", default = "modules",
              help = "dir containing modules", metavar = "character")
)

# Parse command-line arguments
opt = parse_args(OptionParser(option_list = option_list))

# Extract arguments
tissue = opt$tissue
FDR = opt$FDR
output_dir = opt$output_dir
PCO_association_dir = opt$PCO_association_dir

# Initialize results data.table
all_snp_info = fread(paste0(output_dir, "/genotype_files/", opt$genotype_prefix, ".bim"), header = F)

for (fold in 0:opt$folds) {
  sig_results = data.table(module = numeric(), rsIDs = character(), pval = numeric(),FDR = numeric())
  # Determine number of modules
  module_dir = paste0(output_dir, "/transPCO/", tissue, "/fold_", fold, "/",opt$module_dir,"/")
  modules = list.files(module_dir)

  PCO_association_path = paste0(output_dir, "/transPCO/", tissue, "/fold_", fold, "/",opt$PCO_association_dir,"/")

  # Create bed directory
  bed_dir = paste0(output_dir, "/transPCO/", tissue, "/fold_", fold, "/",opt$bed_dir,"/")
  dir.create(bed_dir, recursive = TRUE)
  largest_sig_pval = c()
  for (module in modules) {

    module = strsplit(module,"\\.")[[1]][1]
    print(module) 
    for (chr in 1:22) {
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
        sig_results = rbind(sig_results, data.table(module = module, sig_rows))
      }
    }
    print(sig_results)
    # Read genes in the current module
    genes_in_module = fread(paste0(module_dir, module,".txt"), header = FALSE)
    genes_in_module = unlist(genes_in_module$V1)

    # Identify significant SNPs for each gene in this module
    for (gene in genes_in_module) {
      module_name = module
      sig_snps_per_gene = as.matrix(sig_results[sig_results$module == module_name,c('rsIDs','pval')])
      if (length(sig_snps_per_gene) > 0) {
        fwrite(sig_snps_per_gene, 
               paste0(bed_dir, gene, ".txt"), 
               col.names = FALSE, row.names = FALSE, sep = '\t', quote = FALSE,append = TRUE)
      }
    }
  }
  print(mean(largest_sig_pval,na.rm = T))
  fwrite(sig_results,paste0(output_dir, "/transPCO/",tissue,"/fold_",fold,"/significant_results_",opt$module_dir,"_FDR_",FDR,".txt"), col.names = T, row.names = F, sep = '\t', quote = F)
}


