# --- 8.2_create_EGRET_pos.R ------------------------------------------------
# Writes the FUSION .pos index of EGRET models that passed the significance filter.
# Filters on the best partition's p-value, matching how 6.4 selects the deployed model.
#
#   in : results_sumstats/{tissue}/{egret_subdir}.txt; gene_info
#   out: pos_files/{tissue}/{egret_subdir}.pos
#   run: 8.0_run_TWAS_EGRET.sh
# --------------------------------------------------------------------------
library(data.table)
library(optparse)

option_list = list(
  make_option("--tissue", action="store",default=NA, type='character',
              help="tissue for analysis"),
  make_option("--output_dir", action="store", default=NA, type='character',
              help="Path to output directory"),
  make_option("--gene_info", action="store", default=NA, type='character',
              help="Path to gene info file"),
  make_option("--egret_subdir", action="store", default=NA, type='character',
              help="EGRET model namespace")
  )

opt = parse_args(OptionParser(option_list=option_list))


required = c("tissue", "output_dir", "gene_info", "egret_subdir")
missing = required[!vapply(required, function(o) length(opt[[o]]) == 1L && !is.na(opt[[o]]), logical(1))]
if (length(missing)) stop("missing required argument(s): --", paste(missing, collapse = ", --"), call. = FALSE)

tissue = opt$tissue
output_dir = opt$output_dir

gene_info_file_path = opt$gene_info
all_gene_info = fread(gene_info_file_path, header = T)
sumstats = fread(paste0(output_dir,"/results_sumstats/",tissue,"/",opt$egret_subdir,".txt"),header = T)
sumstats$best_r2_pval = pmin(sumstats$p_gw, sumstats$p_cis_part, sumstats$p_trans_part, na.rm = TRUE)
print(sumstats)
wgt_dir = paste0(output_dir,"/EGRET_models/",tissue,"/",opt$egret_subdir,"/")

pos_rows = list()

for (row in 1:nrow(sumstats)) {
    if (is.na(sumstats$'best_r2_pval'[row]) || sumstats$'best_r2_pval'[row] > 0.01) {
        next
    }

    # match() returns NA for an unknown gene and all_gene_info[NA,] is an all-NA row,
    # which would write CHR=NA into the .pos file. Fail loudly instead.
    gene_row = match(sumstats$gene[row],all_gene_info$geneId)
    if (is.na(gene_row)) {
        stop(sprintf("gene %s is absent from %s", sumstats$gene[row], gene_info_file_path))
    }
    gene_info = all_gene_info[gene_row,]
    gene_chr = strsplit(gene_info$'#chrom',"chr")[[1]][2]
    gene_start = gene_info$chromStart
    gene_end = gene_info$chromEnd

    pos_rows[[length(pos_rows) + 1L]] = data.table(WGT = paste0(sumstats$gene[row],".wgt.RDat"),ID = sumstats$gene[row],CHR = gene_chr,P0 = gene_start, P1 = gene_end)

}
pos_matrix = if (length(pos_rows)) rbindlist(pos_rows) else data.table(WGT = character(), ID = character(), CHR = numeric(), P0 = numeric(), P1 = numeric())

pos_dir = paste0(output_dir,"/pos_files/",tissue,"/")
dir.create(pos_dir,recursive = T)

fwrite(pos_matrix,paste0(pos_dir,opt$egret_subdir,".pos"),sep = '\t', quote = F, col.names = T, row.names = F)
