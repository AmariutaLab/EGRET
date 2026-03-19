library("data.table")
library("optparse")

option_list = list(
        make_option("--gene_information", action="store", default=NA, type='character',
                help="Path to file containing gene information"),
        make_option("--cis_window_size", action="store", default=1000000, type='numeric',
                help="size of cis window to create"),
        make_option("--tissue", action="store", default=NA, type='character',
                help="tissue that will set up for"),
        make_option("--plink_path", action="store", default=NA, type='character',
                help="Path to plink"),
        make_option("--output_dir", action="store", default=NA, type='character',
                help="Base output directory for the pipeline"),
        make_option("--genotype_prefix", action="store", default="GTEX_v8_genotypes_pruned", type='character', help="Basename of LD-pruned plink genotype files in output_dir/genotype_files/")
  )

opt = parse_args(OptionParser(option_list=option_list))

all_gene_info = fread(opt$gene_information, header = T)
expression = fread(file.path(opt$output_dir, "expression_files", paste0(opt$tissue, "_expression.txt.gz")), header = T)

genes = unlist(expression$gene_id)

bfile = file.path(opt$output_dir, "genotype_files", opt$genotype_prefix)
bed_dir = file.path(opt$output_dir, "bed_files", opt$tissue, "cis")
plink_dir = file.path(opt$output_dir, "plink_results", opt$tissue, "cis")

dir.create(bed_dir, recursive = T)
dir.create(plink_dir, recursive = T)

for (gene in genes) {
	gene_info = all_gene_info[grep(gene, all_gene_info$geneId)]
	chr = strsplit(gene_info$'#chrom', split = 'chr')[[1]][2]
	tss = gene_info$chromStart

	print(gene)
	if (chr == 'X' | chr == 'Y') {
  		next
	}

	bounds_u = max(1, tss - floor(opt$cis_window_size/2))
	bounds_l = tss + floor(opt$cis_window_size/2)
	bed = c(chr, bounds_u, bounds_l, gene, 0, "+")
	bed = matrix(bed, nrow = 1, ncol = 6)

	bed_file = file.path(bed_dir, paste0(gene, ".bed"))
	write.table(bed, bed_file, row.names = F, col.names = F, sep = "\t", quote = F)

	plink_out = file.path(plink_dir, gene)
	arg = paste(opt$plink_path, "--bfile", bfile, "--extract bed1", bed_file, "--make-bed --out", plink_out)
	system(arg)
	file.remove(paste0(plink_out, ".log"))
}

