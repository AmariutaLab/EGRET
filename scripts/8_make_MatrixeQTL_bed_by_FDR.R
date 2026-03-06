library('data.table')
library('optparse')

option_list = list(
  make_option("--tissue", action="store", default=NA, type='character',
              help="tissue"),
  make_option("--folds", action="store",default=NA, type='numeric',
              help="number of folds for analysis"),
  make_option("--FDR", action="store", default=0.1, type='numeric',
              help="FDR rate to include snps into the model"),
  make_option("--gene_info", action="store",default=NA, type='character',
        	  help="file containg gene info such as gene id, gene name, start, chromosome, and gene type"),
  make_option("--genotype_output_prefix", action="store",default=NA, type='character',
        	  help="prefix to access bim file within genotypes folder "),
  make_option("--output_dir", action="store",default=NA, type='character',
              help="directory where results are stored")
  )

opt = parse_args(OptionParser(option_list=option_list))
tissue = opt$tissue

all_gene_info = fread(opt$gene_info, header = T)
all_snp_info = fread(paste0(opt$output_dir,"/genotype_files/", opt$genotype_output_prefix ,".bim"), header = F)


gene_expression = fread(paste0(opt$output_dir,"/expression_files/",tissue,"_expression_regressed.txt.gz"),header = T)

for (fold in 0:5) {
	results = fread(paste0(opt$output_dir,"/MatrixeQTL/",tissue,"/association_results_fold_",fold,"_0.0001_threshold.txt"),header = T)
	genes = unique(results$gene)
	output_dir = paste0(opt$output_dir,"/MatrixeQTL/",tissue,"/fold_",fold,"/results_FDR_",opt$FDR,"/")
	largest_sig_pval = c()
	for (current_gene in genes) {
		subset = results[gene == current_gene,]
		subset$new_FDR = p.adjust(subset$'p-value',n = nrow(all_snp_info), method = 'fdr')
		subset = subset[subset$new_FDR < opt$FDR,]
		largest_sig_pval = c(largest_sig_pval, subset$'p-value'[nrow(subset)])
		gene_info = all_gene_info[match(current_gene,all_gene_info$geneId),]
		chr = strsplit(gene_info$'#chrom','chr')[[1]][2]
		u_bound = max(0,gene_info$chromStart - 500000)
		l_bound = gene_info$chromStart + 500000
		snp_info = all_snp_info[match(subset$SNP,all_snp_info$V2),]
		cis_snps = which(snp_info$V1 == chr & snp_info$V4 > u_bound & snp_info$V4 < l_bound)
		subset = subset[!cis_snps,]
		
		dir.create(output_dir,recursive = T)
		snp_matrix = as.matrix(subset[, c("SNP", "p-value")])
		if (nrow(snp_matrix) == 0) { next }
		fwrite(snp_matrix,paste0(output_dir,current_gene,".txt") ,quote = F, col.names = F, row.names = F,sep = '\t')
		
	}
}
