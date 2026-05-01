library("data.table")
library("optparse")

option_list = list(
        make_option("--expression", action="store",default=NA, type='character',
                help="expression file"),
        make_option("--individuals", action="store",default=NA, type='character',
                help="Path to file containing list of individuals to keep"),
	make_option("--covariates", action="store",default=NA, type='character',
                help="Path to file of expresion covariates"),
	make_option("--tissue", action="store",default=NA, type='character',
                help="tissue for analysis"),
        make_option("--gene_info", action="store",default=NA, type='character',
                help="file containg gene info such as gene id, gene name, start, chromosome, and gene type"),
        make_option("--output_dir", action="store",default=NA, type='character',
                help="directory where output files are placed")
  )

opt = parse_args(OptionParser(option_list=option_list))

tissue = opt$tissue

expression_file = fread(opt$expression, header = TRUE)

individuals = fread(opt$individuals, header = FALSE)
individuals = unlist(individuals$V1)

expression_file = expression_file[, c("gene_id", individuals), with = FALSE]

#filter out genes which do not have gene info and genes which are not protein coding, lincRNA, or antisense RNA
all_gene_info = fread(opt$gene_info,header = T)
expression_file$geneType = all_gene_info$geneType[match(expression_file$gene_id,all_gene_info$geneId)]
expression_file = expression_file[expression_file$geneType == 'protein_coding' | expression_file$geneType == 'lincRNA' | expression_file$geneType == 'antisense']
expression_file = expression_file[,-c("geneType")]

# create directory to write expression
dir.create(paste0(opt$output_dir,"/expression_files/"),recursive = T)

#write expression file to expression folder
fwrite(expression_file, paste0(opt$output_dir,"/expression_files/",tissue,"_expression.txt.gz"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')

# read in covariate file
covariates = fread(opt$covariates, header = TRUE)

# match individuals
covariates = covariates[match(individuals, covariates$ID)]

dir.create(paste0(opt$output_dir,"/covariate_files/"),recursive = T)

fwrite(covariates,paste0(opt$output_dir,"/covariate_files/",tissue,"_covariates.txt.gz"),col.names = T, row.names = F, quote = F, sep = '\t')

residual_gene_expression = matrix(nrow = nrow(expression_file), ncol = ncol(expression_file))
colnames(residual_gene_expression) = colnames(expression_file)

for (row in 1:nrow(expression_file)) {
        reg = summary(lm( t(as.matrix(expression_file[row,-c('gene_id')])) ~ as.matrix(covariates[,c(3:ncol(covariates)),with = FALSE]) ))
        residual_gene_expression[row,2:ncol(residual_gene_expression)] = scale(reg$residuals) # scaled residuals
}

residual_gene_expression_df = as.data.frame(residual_gene_expression)
residual_gene_expression_df$gene_id = expression_file$gene_id

# Reorder columns to place gene_id first
residual_gene_expression_df = residual_gene_expression_df[, c("gene_id", setdiff(names(residual_gene_expression_df), "gene_id"))]

# Write to file with fwrite
fwrite(residual_gene_expression_df, paste0(opt$output_dir,"/expression_files/",tissue,"_expression_regressed.txt.gz"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')





