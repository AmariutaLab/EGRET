# --- 5.1_run_WGCNA_clustering.R --------------------------------------------
# Builds WGCNA co-expression modules per fold. Modules are the unit the PCO test
# aggregates over.
# Note this regresses only `covariate_columns_for_coexpression`, a subset, so modules
# are defined on a differently adjusted matrix than the association stages use.
#
#   in : expression_files/{tissue}_expression.txt.gz (raw); covariates; individuals
#   out: transPCO/{tissue}/fold_{i}/modules/module_{k}_gene_list.txt
#   run: 5.0_run_transPCO_stage.sh
# --------------------------------------------------------------------------
library(WGCNA)
library(data.table)
library(optparse)

option_list = list(
        make_option("--folds", action="store", default=NA, type='numeric',
                help="number of cross-validation folds"),
        make_option("--expression", action="store",default=NA, type='character',
                help="expression file"),
        make_option("--output_dir", action="store",default=NA, type='character',
                help="Path to store output"),
        make_option("--individuals", action="store",default=NA, type='character',
                help="Path to file containing list of individuals to keep"),
        make_option("--covariates", action="store",default=NA, type='character',
                help="Path to file of expresion covariates"),
        make_option("--tissue", action="store",default=NA, type='character',
                help="tissue for analysis"),
		make_option("--covariate_columns_for_coexpression", action = 'store', default = NA, type = 'character',
		    	help = 'column names of covariate file to keep for regressing on expression for coexpression analysis'),
		make_option("--gene_info", action="store",default=NA, type='character',
				help="file containg gene info such as gene id, gene name, start, chromosome, and gene type")
  )

opt = parse_args(OptionParser(option_list=option_list))


required = c("folds", "expression", "output_dir", "individuals", "covariates", "tissue", "covariate_columns_for_coexpression", "gene_info")
missing = required[!vapply(required, function(o) length(opt[[o]]) == 1L && !is.na(opt[[o]]), logical(1))]
if (length(missing)) stop("missing required argument(s): --", paste(missing, collapse = ", --"), call. = FALSE)

tissue = opt$tissue
print(tissue)

gene_expression = fread(opt$expression,header = T)
covar = fread(opt$covariates, header = T)

individuals = fread(opt$individuals,header = F, colClasses = "character")
individual_ids = individuals$V1

gene_ids = gene_expression$gene_id
gene_expression = gene_expression[, ..individual_ids]

gene_info = fread(opt$gene_info,header = T)
genes_on_x_y_chr = gene_info[gene_info$`#chrom` == 'chrX' | gene_info$`#chrom` == 'chrY',]

covariate_cols = unlist(strsplit(opt$covariate_columns_for_coexpression, "\\s+"))
print("Regressing covariates now")
expression_matrix = as.matrix(gene_expression)
design_matrix = cbind(1, as.matrix(covar[, ..covariate_cols]))
design_qr = qr(design_matrix)
# qr.resid, not qr.coef -- see 1.2_setup_expression.R. covariate_columns_for_coexpression
# includes `sex`, which is rank-deficient in the sex-specific tissues.
residuals_by_gene = qr.resid(design_qr, t(expression_matrix))
residual_gene_expression = t(scale(residuals_by_gene))
colnames(residual_gene_expression) = unlist(individual_ids)
rownames(residual_gene_expression) = gene_ids

for (fold in 0:opt$folds) {
	individual_ids = fread(paste0(opt$output_dir,"/fold_",fold,"_info/",tissue,"/train_individuals.txt"),header = F, colClasses = "character")
	residual_gene_expression_fold_individuals = residual_gene_expression[,match(individual_ids$V1,colnames(residual_gene_expression))]

	residual_gene_expression_no_xy = residual_gene_expression_fold_individuals[which(!rownames(residual_gene_expression) %in% genes_on_x_y_chr$geneId),]

	power = 5
	print(paste0("The power is ", power))

	temp_cor = cor
	cor = WGCNA::cor

	bwnet = blockwiseModules(t(residual_gene_expression_no_xy),maxBlockSize = 30000,TOMType = 'unsigned',power = power,mergeCutHeight = 0.15,minModuleSize=10,numericLabels = T,randomSeed = 1234,verbose = 3)


	print(dim(table(bwnet$colors)))

	output_dir = paste0(opt$output_dir,"/transPCO/",tissue,"/fold_",fold,"/modules/")
	dir.create(paste0(output_dir),recursive = T)

	# Modules are labelled 0 (grey/unassigned) then 1..K; write the K real ones.
	# seq_len, not 1:(n-1): with only grey, "1:0" counts DOWN and writes an empty module_1.
	n_modules = dim(table(bwnet$colors)) - 1
	for (i in seq_len(n_modules)) {
  		fwrite(as.data.frame(names(bwnet$colors[bwnet$colors == i])),paste0(output_dir,"/module_",i,"_gene_list.txt"),sep = '\t',quote = F,col.names = F)
	}
}
