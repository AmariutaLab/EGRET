# --- 7.2_summarize_cis_models.R --------------------------------------------
# The cis-arm counterpart of 7.1: aggregates the cis FUSION weights into the baseline
# sumstats that the cis-only TWAS and the EGRET-vs-cis comparison use.
#
#   in : FUSION/{tissue}/cis/*.wgt.RDat
#   out: results_sumstats/{tissue}/cis.txt
#   run: 7.0_run_model_analysis.sh
# --------------------------------------------------------------------------
library('data.table')
library('optparse')

option_list = list(
  make_option("--weight_dir", action="store",default=NA, type='character',
              help='Path to weight files'),
  make_option("--output_file", action="store",default=NA, type='character',
              help='Path to output'),
  make_option("--tissue", action="store", default=NA, type='character',
              help="tissue"),
  make_option("--output_dir", action="store", default=NA, type='character',
              help="directory where outputs from EGRET are stored")

  )
opt = parse_args(OptionParser(option_list=option_list))

required = c("weight_dir", "output_file", "tissue", "output_dir")
missing = required[!vapply(required, function(o) length(opt[[o]]) == 1L && !is.na(opt[[o]]), logical(1))]
if (length(missing)) stop("missing required argument(s): --", paste(missing, collapse = ", --"), call. = FALSE)

tissue = opt$tissue
weight_dir = opt$weight_dir
output_file = opt$output_file
output_dir = opt$output_dir

weight_files = list.files(weight_dir)

model_results = c()

for (gene in weight_files) {
	print(gene)
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

dir.create(paste0(output_dir,"/results_sumstats/",opt$tissue,"/"),recursive = T)
output_file = paste0(output_dir,"/results_sumstats/",opt$tissue,"/",output_file)

write.table(model_results,file = output_file,col.names = TRUE,row.names = FALSE, sep = "\t",quote= FALSE)






