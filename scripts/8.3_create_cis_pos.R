library(data.table)
library(optparse)

option_list = list(
  make_option("--tissue", action="store",default=NA, type='character',
              help="tissue for analysis"),
  make_option("--output_dir", action="store", default=NA, type='character',
              help="Path to output directory"),
  make_option("--gene_info_file_path", action="store", default="../data/GTEx_V8.txt.gz", type='character',
              help="Path to gene info file")
)
opt = parse_args(OptionParser(option_list=option_list))

if (is.na(opt$tissue)) {
        print("no tissue specified")
        q()
} else {
        tissue = opt$tissue
}
output_dir = opt$output_dir

gene_info_file_path = opt$gene_info_file_path
all_gene_info = fread(gene_info_file_path, header = T)

sumstats = fread(paste0(output_dir,"/results_sumstats/",tissue,"/cis.txt"),header = T)

wgt_dir = paste0(output_dir,"/FUSION/",tissue,"/cis/")  

pos_matrix = data.frame(WGT = as.character(),ID = as.character(),CHR = as.numeric(), P0 = as.numeric(),P1 = as.numeric())

for (row in 1:nrow(sumstats)) {
    print(row)
    if (is.na(sumstats$'r2_pval'[row]) || sumstats$'r2_pval'[row] > 0.01) {
        next
    }

    gene_info = all_gene_info[grep(sumstats$gene[row],all_gene_info$geneId),]
    gene_chr = strsplit(gene_info$'#chrom',"chr")[[1]][2]
    gene_start = gene_info$chromStart
    gene_end = gene_info$chromEnd

    pos_matrix = rbind(pos_matrix,data.frame(WGT = paste0(sumstats$gene[row],".wgt.RDat"),ID = sumstats$gene[row],CHR = gene_chr,P0 = gene_start, P1 = gene_end))

}

pos_dir = paste0(output_dir,"/pos_files/",tissue,"/")
dir.create(pos_dir,recursive = T)

fwrite(pos_matrix,paste0(pos_dir,"cis.pos"),sep = '\t', quote = F, col.names = T, row.names = F)
