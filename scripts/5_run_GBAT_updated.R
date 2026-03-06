library('data.table')
library('optparse')
library('plink2R')

option_list = list(
  make_option("--tissue", action="store",default=NA, type='character',
              help="tissue for analysis"),
  make_option("--cis_model_dir", action="store",default=NA, type='character',
              help="directory containing previously trained FUSION models"),
  make_option("--gene_info", action="store",default=NA, type='character',
        	  help="file containg gene info such as gene id, gene name, start, chromosome, and gene type"),
  make_option("--plink_path", action="store",default=NA, type='character',
        	  help="path to plink executable"),
  make_option("--genotype_file_path", action="store",default=NA, type='character',
        	  help="path to bed/bim/fam to use for cis model imputation"),		  
  make_option("--output_dir", action="store",default=NA, type='character',
              help="directory where results are stored")
  )


opt = parse_args(OptionParser(option_list=option_list))

print(opt$output_dir)
if (is.na(opt$tissue)) {
	print("no tissue specified")
	q()
} else {
	tissue = opt$tissue
}

gene_info = fread(opt$gene_info, header = T, sep = '\t')

trans_genes = list.files(paste0(opt$cis_model_dir),pattern = ".RDat")
cis_predicted_expression = list()

individuals_path = paste0(opt$output_dir,"/fold_0_info/",tissue,"/train_individuals.txt")
individual_ids = fread(individuals_path, header = F)

for (trans_gene in trans_genes) {
	gene_name = paste0(strsplit(trans_gene,"\\.")[[1]][1:2],collapse = ".")
	print(gene_name)

    load(paste0(opt$cis_model_dir,trans_gene))
	model_pval = cv.performance[2,colnames(cv.performance) == "lasso"]

	if (is.na(model_pval) | model_pval > 0.01) {
		next
	}

	print(paste("the r2 is",cv.performance[1,1]))

 	weights = wgt.matrix[, c('lasso')]
    non_zero_weights = weights[weights != 0]

    if (length(non_zero_weights) == 0) {
        cat("Skipping", gene, "- no non-zero weights\n")
        next
    }

	bed_file_dir = paste0(opt$output_dir,"/bed_files/",tissue,"/cis/")
	plink_output_dir = paste0(opt$output_dir,"/plink_results/",tissue,"/cis/")
	dir.create(bed_file_dir,recursive = T)
	dir.create(plink_output_dir,recursive = T)

    if(is.null(names(weights))) {next}
    fwrite(as.matrix(names(non_zero_weights)), paste0(bed_file_dir,"/",gene_name,".bed"), row.names = F, col.names = F, quote = F)

	arg = paste0(opt$plink_path," --bfile ",opt$genotype_file_path," --keep ",individuals_path," --extract ",bed_file_dir,gene_name, ".bed --make-bed --out ",plink_output_dir,gene_name)
	system(arg,ignore.stdout = T)
	print("PLINK is DONE")

    if (!file.exists(paste0(plink_output_dir, gene_name,".bed"))) next
    plink_data = read_plink(paste0(plink_output_dir, gene_name),impute='avg')
	
    genos = scale(plink_data$bed)

	matching_indices = match(names(non_zero_weights), colnames(genos), nomatch=0)
    valid_indices = which(matching_indices > 0)

	if (length(valid_indices) != 0) {
        genotypes_subset = genos[, matching_indices[valid_indices], drop=FALSE]

        # Compute gene expression
        expression = as.numeric(genotypes_subset %*% non_zero_weights[valid_indices])
        cis_predicted_expression[[gene_name]] = expression  # Store as a named list

    } else {
        cat("Skipping total expression of ", gene_name, "- no matching SNPs found in genotypes\n")
    }

}
cis_predicted_expression_matrix = do.call(cbind, cis_predicted_expression)
genes = colnames(cis_predicted_expression_matrix)
cis_predicted_expression_matrix = t(cis_predicted_expression_matrix)
colnames(cis_predicted_expression_matrix) = unlist(individual_ids)

cis_predicted_expression_matrix = as.data.frame(cis_predicted_expression_matrix)

cis_predicted_expression_matrix = cbind(genes,cis_predicted_expression_matrix)


dir.create(paste0(opt$output_dir,"/GBAT/",tissue,"/"), recursive = T)
fwrite(cis_predicted_expression_matrix,paste0(opt$output_dir,"/GBAT/",tissue,"/cis_predicted_expression.txt"), col.names = T, row.names = F, sep = '\t', quote = F)

