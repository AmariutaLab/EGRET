# --- 7.1_summarize_EGRET_models.R ------------------------------------------
# Aggregates every EGRET model into one table, reporting R2 and p separately for the
# gw, cis and trans partitions. p_trans_part is what the GRN stage filters on.
#
#   in : EGRET_models/{tissue}/{egret_subdir}/*.wgt.RDat
#   out: results_sumstats/{tissue}/{egret_subdir}.txt (r2_gw, r2_cis_part, r2_trans_part, p_trans_part)
#   run: 7.0_run_model_analysis.sh
# --------------------------------------------------------------------------
library('data.table')
library('optparse')


option_list = list(
	make_option("--tissue", action="store", default=NA, type='character',
              help="tissue for analysis"),
	make_option("--egret_subdir", action = "store", default=NA, type = 'character',
		    help = "the name of the analysis"),
	make_option("--output_dir", action = "store", default=NA, type = 'character',
		    help = "Path to output directory")
	)

opt = parse_args(OptionParser(option_list=option_list))

required = c("tissue", "egret_subdir", "output_dir")
missing = required[!vapply(required, function(o) length(opt[[o]]) == 1L && !is.na(opt[[o]]), logical(1))]
if (length(missing)) stop("missing required argument(s): --", paste(missing, collapse = ", --"), call. = FALSE)

tissue = opt$tissue


project = opt$egret_subdir
output_dir = opt$output_dir
genes = list.files(paste0(output_dir,"/EGRET_models/",tissue,"/",project,"/"))

result_rows = list()


for (gene in unlist(genes)) {
	gene = substr(gene, 1, nchar(gene) - 9)
	print(gene)

	if (file.exists(paste0(output_dir,"/EGRET_models/",tissue,"/",project,"/",gene,".wgt.RDat"))) {
                load(paste0(output_dir,"/EGRET_models/",tissue,"/",project,"/",gene,".wgt.RDat"))
		gw_best = which.max(ifelse(is.na(cv.performance_gw[1,]), 0, cv.performance_gw[1,]))
		r2_gw = cv.performance_gw[1,gw_best]
		p_gw = cv.performance_gw[2,gw_best]
		cis_part_best = which.max(ifelse(is.na(cv.performance_cis[1,]), 0, cv.performance_cis[1,]))
		r2_cis_part = cv.performance_cis[1,cis_part_best]
		p_cis_part = cv.performance_cis[2,cis_part_best]
		trans_part_best = which.max(ifelse(is.na(cv.performance_trans[1,]), 0, cv.performance_trans[1,]))
		r2_trans_part = cv.performance_trans[1,trans_part_best]
		p_trans_part = cv.performance_trans[2,trans_part_best]
        	best_model = colnames(cv.performance)[which.max(ifelse(is.na(cv.performance[1,]), 0, cv.performance[1,]))]
	} else { r2_trans = "NA"}
	result_rows[[length(result_rows) + 1L]] = data.table(gene,r2_gw,p_gw, r2_cis_part,p_cis_part, r2_trans_part, p_trans_part, best_model)
}
results = if (length(result_rows)) {
	rbindlist(result_rows)
} else {
	data.table(gene = character(), r2_gw = numeric(), p_gw = numeric(), r2_cis_part = numeric(),
		p_cis_part = numeric(), r2_trans_part = numeric(), p_trans_part = numeric(), best_model = character())
}
dir.create(paste0(output_dir,"/results_sumstats/",tissue,"/"), recursive = T)
fwrite(results,paste0(output_dir,"/results_sumstats/",tissue,"/",project,".txt"), col.names = T, row.names = F, quote = F, sep = '\t')



