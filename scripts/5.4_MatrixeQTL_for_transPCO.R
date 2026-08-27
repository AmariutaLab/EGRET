library('MatrixEQTL')
library('data.table')
library('optparse')
library('plink2R')

# Parse arguments
option_list = list(
  make_option("--tissue", action="store", default=NA, type='character',
              help="tissue"),
  make_option("--module", action="store", default=NA, type='character',
              help="module file name"),
  make_option("--module_dir", action="store", default="modules", type='character',
              help="subdirectory name containing module gene lists (within transPCO/tissue/fold_N/)"),
  make_option("--output_dir", action="store", default=NA, type='character',
              help="base output directory for the pipeline"),
  make_option("--folds", action="store", default=5, type='integer',
              help="number of cross-validation folds (iterates fold_0 through fold_N)"),
  make_option("--plink_path", action="store", default=NA, type='character',
              help="path to plink executable (used as fallback if genotype txt files are missing)"),
  make_option("--genotype_prefix", action="store", default="GTEX_v8_genotypes_pruned", type='character',
              help="basename of LD-pruned plink genotype files in output_dir/genotype_files/")
)
opt = parse_args(OptionParser(option_list=option_list))

# Define input variables
tissue = opt$tissue
module = opt$module
module_dir = opt$module_dir
base_dir = opt$output_dir

for (fold in 0:opt$folds) {
  module_path = file.path(base_dir, "transPCO", tissue, paste0("fold_", fold), module_dir, module)

  if (!file.exists(module_path)) next

  genes_in_module = fread(module_path, header = FALSE)

  # determine number of chromosomes
  bim_file = fread(paste0(opt$output_dir,"/genotype_files/",opt$genotype_prefix,".bim"), header = F, sep = '\t')
  number_chr = max(bim_file$V1)

  for (chr in 1:number_chr) {
    module_name = strsplit(module, "\\..")[[1]][1]
    result_file = file.path(base_dir, "transPCO", tissue, paste0("fold_", fold), "association_results",
                           paste0(module_name, "_chr_", chr, ".txt.gz"))
    if (file.exists(result_file)) next

    # Prepare gene expression file
    gene_expression_path = file.path(base_dir, paste0("fold_", fold, "_info"), tissue, "train_expression.txt")
    gene_expression = fread(gene_expression_path, header = TRUE)

    gene_expression_module = gene_expression[match(genes_in_module$V1, gene_expression$gene_id), ]
    expr_output_dir = file.path(base_dir, "transPCO", tissue, paste0("fold_", fold), "expression_files")
    dir.create(expr_output_dir, recursive = TRUE, showWarnings = FALSE)
    expr_output_file = file.path(expr_output_dir, paste0(module_name, ".txt"))
    fwrite(gene_expression_module, expr_output_file, sep = '\t', quote = FALSE, col.names = TRUE)

    # Load gene expression
    gene = SlicedData$new()
    gene$fileDelimiter = "\t"
    gene$fileOmitCharacters = "NA"
    gene$fileSkipRows = 1
    gene$fileSkipColumns = 1
    gene$fileSliceSize = 2000
    gene$LoadFile(expr_output_file)

    # Prepare genotype file
    geno_output_dir = file.path(base_dir, "transPCO", tissue, paste0("fold_", fold), "genotype_files")
    geno_output_file = file.path(geno_output_dir, paste0("train_genotypes_chr_", chr, ".txt"))

    if (!file.exists(geno_output_file)) {
      dir.create(geno_output_dir, recursive = TRUE, showWarnings = FALSE)
      train_individuals_path = file.path(base_dir, paste0("fold_", fold, "_info"), tissue, "train_individuals.txt")
      plink_cmd = paste(opt$plink_path,
                        "--bfile", file.path(base_dir, "genotype_files", opt$genotype_prefix),
                        "--make-bed --keep", train_individuals_path,
                        "--chr", chr,
                        "--out", file.path(geno_output_dir, paste0("train_genotypes_chr_", chr)))
      system(plink_cmd)

      genotypes = read_plink(file.path(geno_output_dir, paste0("train_genotypes_chr_", chr)), impute='avg')
      genotypes$bed = scale(genotypes$bed)
      genotypes_transposed = as.data.table(t(genotypes$bed), keep.rownames = "SNP")
      fwrite(genotypes_transposed, geno_output_file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = '\t')
    }

    # Load genotypes
    snps = SlicedData$new()
    snps$fileDelimiter = "\t"
    snps$fileOmitCharacters = "NA"
    snps$fileSkipRows = 1
    snps$fileSkipColumns = 1
    snps$fileSliceSize = 20000
    snps$LoadFile(geno_output_file)

    # Run Matrix eQTL
    assoc_output_dir = file.path(base_dir, "transPCO", tissue, paste0("fold_", fold), "association_results")
    dir.create(assoc_output_dir, recursive = TRUE, showWarnings = FALSE)
    assoc_output_file = file.path(assoc_output_dir, paste0(module_name, "_chr_", chr, ".txt"))
    
    Matrix_eQTL_main(
      gene = gene,
      snps = snps,
      output_file_name = assoc_output_file,
      pvOutputThreshold = 1,
      useModel = modelLINEAR,
      errorCovariance = numeric(),
      verbose = TRUE,
      pvalue.hist = TRUE,
      min.pv.by.genesnp = FALSE,
      noFDRsaveMemory = FALSE
    )

    # Transform results
    results = fread(assoc_output_file, header = TRUE)
    matrix_form = dcast(results, SNP ~ gene, value.var = "t-stat", fun.aggregate = mean, drop = FALSE)
    matrix_output_file = file.path(assoc_output_dir, paste0("matrix_",module_name, "_chr_", chr, ".txt.gz"))
    fwrite(matrix_form, matrix_output_file, sep = "\t", row.names = FALSE, col.names = TRUE)

    # Cleanup
    file.remove(assoc_output_file)
    rm(list = c("matrix_form", "results", "snps", "gene", "genotypes", "genotypes_transposed"))
  }
}

