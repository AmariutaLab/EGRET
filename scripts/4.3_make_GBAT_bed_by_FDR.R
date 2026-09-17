# --- 4.3_make_GBAT_bed_by_FDR.R --------------------------------------------
# FDR-filters the GBAT associations and writes one BED file of significant trans
# SNPs per gene, taking each regulator's SNPs from its cis model.
#
#   in : GBAT/{tissue}/fold_{i}/ association results; FUSION/{tissue}/cis/*.wgt.RDat
#   out: GBAT/{tissue}/fold_{i}/results_FDR_{fdr}/{gene}.txt
#   run: 4.0_run_GBAT_stage.sh
# --------------------------------------------------------------------------
library('data.table')
library('optparse')

option_list = list(
  make_option("--folds", action="store", default=NA, type='numeric',
              help="number of cross-validation folds"),
  make_option("--tissue", action="store",default=NA, type='character',
              help="tissue for analysis"),
  make_option("--gene_info", action="store",default=NA, type='character',
        	  help="file containg gene info such as gene id, gene name, start, chromosome, and gene type"),
  make_option("--cis_model_dir", action="store",default=NA, type='character',
              help="directory containing previously trained FUSION models"),
  make_option("--FDR", action="store", default=NA, type='numeric',
              help="FDR to include gene pairs"),
  make_option("--output_dir", action="store",default=NA, type='character',
              help="directory where results are stored"),
  make_option("--cis_window_size", action="store", default=1000000, type='numeric',
              help="total width of the cis region in bp; halved to a +/- half-width around the TSS. Gene pairs closer than the half-width are cis and are skipped [default: 1000000]")
  )

opt = parse_args(OptionParser(option_list=option_list))


required = c("folds", "tissue", "gene_info", "cis_model_dir", "FDR", "output_dir")
missing = required[!vapply(required, function(o) length(opt[[o]]) == 1L && !is.na(opt[[o]]), logical(1))]
if (length(missing)) stop("missing required argument(s): --", paste(missing, collapse = ", --"), call. = FALSE)

tissue = opt$tissue

# Normalise the trailing slash: uses below are paste0(cis_model_dir, gene, ...), so
# ".../cis" would build ".../cisENSG....RDat". Do not rely on the caller.
cis_model_dir = sub("/*$", "/", opt$cis_model_dir)
all_gene_info = fread(opt$gene_info, header = T)
genes = fread(paste0(opt$output_dir,"/expression_files/",tissue,"_expression_regressed.txt.gz"), header = T)
for (fold in 0:opt$folds) {
	output_dir = paste0(opt$output_dir,"/GBAT/",tissue,"/fold_",fold,"/results_FDR_",opt$FDR,"/")
	dir.create(output_dir)
 	association_results = fread(paste0(opt$output_dir,"/GBAT/",tissue,"/fold_",fold,"/association_results.txt"),header = T)
    if (nrow(association_results) == 0 ) { next }
	unique_genes = unique(association_results$gene)
	largest_sig_pval = c()
	for (gene_name in unique_genes) {
		print(gene_name)
		subset = association_results[gene == gene_name,]

		subset$new_FDR = p.adjust(subset$'p-value',n = nrow(genes), method = 'fdr')
		subset = subset[subset$new_FDR < opt$FDR,]
		largest_sig_pval = c(largest_sig_pval, subset$'p-value'[nrow(subset)])
		if (nrow(subset) == 0) {
			next
		}

		gene_info = all_gene_info[ all_gene_info$geneId %in% gene_name]
		gene_chr = gene_info$'#chrom'
		gene_start = gene_info$chromStart
		sig_snp_parts = list()
		for (row in 1:nrow(subset)) {
			if (subset$SNP[row] == gene_name) { next }
			trans_gene_info = all_gene_info[all_gene_info$geneId %in% subset$SNP[row]]
			trans_gene_chr = trans_gene_info$'#chrom'
			trans_gene_start = trans_gene_info$chromStart

			if (any((gene_chr == trans_gene_chr) & abs(gene_start - trans_gene_start) < floor(opt$cis_window_size/2))) { next }
			
			load(paste0(cis_model_dir,subset$SNP[row],".wgt.RDat"))
			snps = snps$V2[wgt.matrix[,colnames(wgt.matrix) == "lasso"] != 0]
			pvals = rep(subset$'p-value'[row],length(snps))
			sig_snp_parts[[length(sig_snp_parts) + 1L]] = data.table(snps,pvals)
		}
		sig_snps = rbindlist(sig_snp_parts)
		if (nrow(sig_snps) == 0) {next}
		sig_snps = as.matrix(sig_snps)
		fwrite(sig_snps,paste0(output_dir,gene_name,".txt"),row.names = F, col.names = F, quote = F, sep = '\t')
    }
	print(mean(largest_sig_pval, na.rm = T))
}

