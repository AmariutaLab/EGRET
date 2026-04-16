library('data.table')
library('optparse')


option_list = list(
	make_option("--tissue", action="store", default=NA, type='character',
              help="tissue for analysis"),
	make_option("--project", action = "store", default=NA, type = 'character',
		    help = "the name of the analysis"),
	make_option("--output_dir", action = "store", default=NA, type = 'character',
		    help = "Path to output directory")
	)

opt = parse_args(OptionParser(option_list=option_list))
tissue = opt$tissue


project = opt$project
output_dir = opt$output_dir
genes = list.files(paste0(output_dir,"/xtune_fusion_models/",tissue,"/",project,"/"))

results = data.frame(gene = as.character(),  r2_gw = as.numeric(),p_gw = as.numeric(), r2_cis_part = as.numeric(),p_cis_part = as.numeric(), r2_trans_part = as.numeric(), p_trans_part = as.numeric(),best_model = as.character())


for (gene in unlist(genes)) {
	gene = substr(gene, 1, nchar(gene) - 9)
	print(gene)

	if (file.exists(paste0(output_dir,"/xtune_fusion_models/",tissue,"/",project,"/",gene,".wgt.RDat"))) {
                load(paste0(output_dir,"/xtune_fusion_models/",tissue,"/",project,"/",gene,".wgt.RDat"))
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
	results = rbind(results,data.frame(gene,r2_gw,p_gw, r2_cis_part,p_cis_part, r2_trans_part, p_trans_part, best_model))
}
dir.create(paste0(output_dir,"/results_sumstats/",tissue,"/"), recursive = T)
fwrite(results,paste0(output_dir,"/results_sumstats/",tissue,"/",project,".txt"), col.names = T, row.names = F, quote = F, sep = '\t')




