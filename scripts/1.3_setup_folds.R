# --- 1.3_setup_folds.R -----------------------------------------------------
# Assigns individuals to cross-validation folds and materialises per-fold training
# genotype and expression matrices.
# `--folds N` writes N+1 directories. fold_0 is NOT a CV fold: its test set is empty
# and it trains on everyone, supplying the final deployed fit. folds 1..N are the
# real holdouts, each training on (N-1)/N of the cohort.
#
#   in : expression_files/{tissue}_expression_regressed.txt.gz; individuals; genotype_files/
#   out: fold_{i}_info/{tissue}/{train,test}_individuals.txt, train_expression.txt, train_genotypes.txt
#   run: 1.0_setup_genotype_and_expression.sh
# --------------------------------------------------------------------------
library('data.table')
library('optparse')
library('plink2R')

option_list = list(
	make_option("--expression", action="store",default=NA, type='character',
			help="Path to gene expression matrix"),
	make_option("--individuals", action="store",default=NA, type='character',
            help="Path to file containing list of individuals"),
    make_option("--tissue", action = "store", default = NA, type = 'character',
			help="tissue that will set up for"),
	make_option("--folds", action="store",default=NA, type='numeric',
            help="Number of folds"),
	make_option("--plink_path", action="store",default=NA, type='character',
            help="Path to plink"),
	 make_option("--genotype_prefix", action="store",default=NA, type='character',
            help="Path to genotype files"),
    make_option("--output_dir", action="store",default=NA, type='character',
        	help="directory where output files are placed")
  )

opt = parse_args(OptionParser(option_list=option_list))


required = c("expression", "individuals", "tissue", "folds", "plink_path", "genotype_prefix", "output_dir")
missing = required[!vapply(required, function(o) length(opt[[o]]) == 1L && !is.na(opt[[o]]), logical(1))]
if (length(missing)) stop("missing required argument(s): --", paste(missing, collapse = ", --"), call. = FALSE)

tissue = opt$tissue
cross_val_folds = opt$folds

gene_expression = fread(opt$expression,header = T)

### SETUP CROSS-VALIDATION INFO ###
all_individuals = fread(opt$individuals,header = F, colClasses = "character")
N = nrow(all_individuals)
folds = cut(seq(1,N),breaks = cross_val_folds,labels=FALSE)

cross_val_results = matrix(NA,nrow=N,ncol=1)
for (i in 0:cross_val_folds) {
	train_individuals = all_individuals[which(folds != i,arr.ind=TRUE),]
	test_individuals = all_individuals[which(folds == i,arr.ind=TRUE),]

	output_dir = paste0(opt$output_dir,"/fold_",i,"_info/",tissue,"/")
	dir.create(output_dir,recursive = T)

	fwrite(train_individuals, paste0(output_dir,"/train_individuals.txt"),row.names = F, col.names = F, quote = F, sep = '\t')
	fwrite(test_individuals, paste0(output_dir,"/test_individuals.txt"),row.names = F, col.names = F, quote = F, sep = '\t')

	train_gene_expression = gene_expression[, c(1,match(as.character(unlist(train_individuals)), as.character(colnames(gene_expression)))), with = FALSE]
	fwrite(train_gene_expression,paste0(output_dir,"/train_expression.txt"),col.names = T,row.names = F,quote = F, sep = '\t')

	plink_out = paste0(output_dir, "/train_genotypes")
	status = system2(opt$plink_path, args = c("--bfile", shQuote(opt$genotype_prefix),
		"--make-bed", "--keep", shQuote(paste0(output_dir, "/train_individuals.txt")),
		"--out", shQuote(plink_out)))
	if (status != 0) {
		stop(sprintf("PLINK fold extraction failed for fold_%s with exit status %s", i, status))
	}
	expected = paste0(plink_out, c(".bed", ".bim", ".fam"))
	if (any(!file.exists(expected)) || any(file.info(expected)$size == 0)) {
		stop(sprintf("PLINK fold extraction did not create a complete fileset for fold_%s", i))
	}

	genotypes = read_plink(plink_out,impute='avg')
	genotypes$bed = scale(genotypes$bed)

	rownames = colnames(genotypes$bed)
	genotypes_transposed = t(genotypes$bed)
	rownames(genotypes_transposed) = rownames
	genotypes_transposed_dt = as.data.table(genotypes_transposed, keep.rownames = "SNP")

	fwrite(genotypes_transposed_dt,paste0(output_dir,"/train_genotypes.txt"),row.names = F,col.names = T,quote = F, sep = '\t')

	# The binary set exists only for the read_plink() above; the dense-text
	# train_genotypes.txt is what MatrixeQTL consumes. 3.0 cleans that one up later.
	unlink(paste0(output_dir,"/train_genotypes",c(".bed",".bim",".fam",".log",".nosex")))
}



