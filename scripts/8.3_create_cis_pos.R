# --- 8.3_create_cis_pos.R --------------------------------------------------
# The cis-arm counterpart of 8.2.
#
#   in : results_sumstats/{tissue}/cis.txt; gene_info
#   out: pos_files/{tissue}/cis.pos
#   run: 8.1_run_TWAS_cis.sh
# --------------------------------------------------------------------------
library(data.table)
library(optparse)

option_list = list(
  make_option("--tissue", action="store",default=NA, type='character',
              help="tissue for analysis"),
  make_option("--output_dir", action="store", default=NA, type='character',
              help="Path to output directory"),
  make_option("--gene_info", action="store", default=NA, type='character',
              help="Path to gene info file")
)
opt = parse_args(OptionParser(option_list=option_list))


required = c("tissue", "output_dir", "gene_info")
missing = required[!vapply(required, function(o) length(opt[[o]]) == 1L && !is.na(opt[[o]]), logical(1))]
if (length(missing)) stop("missing required argument(s): --", paste(missing, collapse = ", --"), call. = FALSE)

tissue = opt$tissue
output_dir = opt$output_dir

gene_info_file_path = opt$gene_info
all_gene_info = fread(gene_info_file_path, header = T)

sumstats = fread(paste0(output_dir,"/results_sumstats/",tissue,"/cis.txt"),header = T)

wgt_dir = paste0(output_dir,"/FUSION/",tissue,"/cis/")  

pos_rows = list()

for (row in 1:nrow(sumstats)) {
    if (is.na(sumstats$'r2_pval'[row]) || sumstats$'r2_pval'[row] > 0.01) {
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

fwrite(pos_matrix,paste0(pos_dir,"cis.pos"),sep = '\t', quote = F, col.names = T, row.names = F)
