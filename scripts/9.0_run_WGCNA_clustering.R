library(WGCNA)
library(data.table)
library(optparse)

option_list = list(
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
		make_option("--num_PCs", action = 'store', default = 10, type = 'numeric',
		    	help = 'number of gene expression PCs to regress out'),
		make_option("--gene_info", action="store",default=NA, type='character',
				help="file containg gene info such as gene id, gene name, start, chromosome, and gene type")
  )

opt = parse_args(OptionParser(option_list=option_list))

tissue = opt$tissue
print(tissue)

gene_expression = fread(opt$expression,header = T)
covar = fread(opt$covar, header = T)

individuals = fread(opt$individuals,header = F)
individual_ids = individuals$V1

gene_ids = gene_expression$gene_id
gene_expression = gene_expression[, ..individual_ids]

gene_info = fread(opt$gene_info,header = T)
genes_on_x_y_chr = gene_info[gene_info$`#chrom` == 'chrX' | gene_info$`#chrom` == 'chrY',]

residual_gene_expression = matrix(nrow = nrow(gene_expression), ncol = ncol(gene_expression))
colnames(residual_gene_expression) = unlist(individual_ids)
rownames(residual_gene_expression) = t(gene_ids)

PC_column = 7 + opt$num_PCs
print("Regressing covariates now")
for (row in 1:nrow(gene_expression)) {
	reg = summary(lm( as.matrix(t(gene_expression[row,])) ~ as.matrix(covar[,c(3:..PC_column,68:70)]) ))
  	residual_gene_expression[row,] = scale(reg$residuals) # scaled residuals
}

for (fold in 0:5) {
	individual_ids = fread(paste0(opt$output_dir,"/fold_",fold,"_info/",tissue,"/train_individuals.txt"),header = F)
	residual_gene_expression_fold_individuals = residual_gene_expression[,match(individual_ids$V1,colnames(residual_gene_expression))]

	residual_gene_expression_no_xy = residual_gene_expression_fold_individuals[which(!rownames(residual_gene_expression) %in% genes_on_x_y_chr$geneId),]

	#power = c(c(1:10),seq(from = 12, to = 50, by = 2))
	#sft = pickSoftThreshold(t(residual_gene_expression_no_xy),powerVector = power,networkType = 'signed',verbose = 5)
	#power = sft$powerEstimate
	power = 5
	print(paste0("The power is ", power))

	temp_cor = cor
	cor = WGCNA::cor

	bwnet = blockwiseModules(t(residual_gene_expression_no_xy),maxBlockSize = 30000,TOMType = 'unsigned',power = power,mergeCutHeight = 0.15,minModuleSize=10,numericLabels = T,randomSeed = 1234,verbose = 3)


	print(dim(table(bwnet$colors)))

	output_dir = paste0(opt$output_dir,"/transPCO/",tissue,"/fold_",fold,"/modules/")
	dir.create(paste0(output_dir),recursive = T)

	for (i in 1:(dim(table(bwnet$colors))-1)) {
  		fwrite(as.data.frame(names(bwnet$colors[bwnet$colors == i])),paste0(output_dir,"/module_",i,"_gene_list.txt"),sep = '\t',quote = F,col.names = F)
	}
}
