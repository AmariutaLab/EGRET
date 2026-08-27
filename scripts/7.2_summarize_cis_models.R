library('data.table')
library('optparse')

option_list = list(
  make_option("--weight_dir", action="store",default=NA, type='character',
              help='Path to weight files'),
  make_option("--output_file", action="store",default=NA, type='character',
              help='Path to output'),
  make_option("--tissue", action="store", default=NA, type='character',
              help="tissue")

  )
opt = parse_args(OptionParser(option_list=option_list))
tissue = opt$tissue
weight_dir = opt$weight_dir
output_file = opt$output_file

weight_files = list.files(weight_dir)

model_results = c()

for (gene in weight_files) {

	gene_name = substr(gene,1,nchar(gene)-9)
	weight_file_path = file.path(weight_dir,gene)
	load(weight_file_path)

	h2 = hsq[1]
	h2_pval = hsq.pv[1]

	r2 = max(cv.performance[1,])
	r2_pval = min(cv.performance[2,])

	best_model = colnames(cv.performance)[which.max(cv.performance[1,])]

	model_results = rbind(model_results,c(gene_name,h2,h2_pval,r2,r2_pval,best_model))
}
colnames(model_results) <- c("gene","h2","h2_pval","r2","r2_pval","best_model")

dir.create(paste0("results_sumstats/",opt$tissue,"/"),recursive = T)
output_file = paste0("results_sumstats/",opt$tissue,"/",output_file)

write.table(model_results,file = output_file,col.names = TRUE,row.names = FALSE, sep = "\t",quote= FALSE)






