# --- 2.1_setup_cis_windows.R -----------------------------------------------
# Extracts a plink fileset per gene covering the cis window (default +/-500kb of
# the TSS), which is what FUSION trains the cis model on.
#
#   in : gene_info; genotype_files/{genotype_prefix}
#   out: plink_results/{tissue}/cis/{gene}.bed/.bim/.fam
#   run: 2.0_run_cis_models.sh
# --------------------------------------------------------------------------
library("data.table")
library("optparse")

option_list = list(
        make_option("--gene_info", action="store", default=NA, type='character',
                help="Path to file containing gene information"),
        make_option("--cis_window_size", action="store", default=1000000, type='numeric',
                help="size of cis window to create"),
        make_option("--tissue", action="store", default=NA, type='character',
                help="tissue that will set up for"),
        make_option("--plink_path", action="store", default=NA, type='character',
                help="Path to plink"),
        make_option("--output_dir", action="store", default=NA, type='character',
                help="Base output directory for the pipeline"),
        make_option("--genotype_prefix", action="store", default=NA, type='character', help="Basename of LD-pruned plink genotype files in output_dir/genotype_files/")
  )

opt = parse_args(OptionParser(option_list=option_list))


required = c("gene_info", "tissue", "plink_path", "output_dir", "genotype_prefix")
missing = required[!vapply(required, function(o) length(opt[[o]]) == 1L && !is.na(opt[[o]]), logical(1))]
if (length(missing)) stop("missing required argument(s): --", paste(missing, collapse = ", --"), call. = FALSE)

all_gene_info = fread(opt$gene_info, header = T)
expression = fread(file.path(opt$output_dir, "expression_files", paste0(opt$tissue, "_expression.txt.gz")), header = T)

genes = unlist(expression$gene_id)

bfile = file.path(opt$output_dir, "genotype_files", opt$genotype_prefix)
bed_dir = file.path(opt$output_dir, "bed_files", opt$tissue, "cis")
plink_dir = file.path(opt$output_dir, "plink_results", opt$tissue, "cis")

dir.create(bed_dir, recursive = T)
dir.create(plink_dir, recursive = T)
manifest = file.path(plink_dir, "cis_genes.txt")
if (file.exists(manifest)) file.remove(manifest)
completed_genes = character(length(genes))
num_completed = 0L

for (gene in genes) {
	gene_info = all_gene_info[match(gene, all_gene_info$geneId)]
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
	status = system2(opt$plink_path, args = c("--bfile", shQuote(bfile), "--extract", "bed1",
		shQuote(bed_file), "--make-bed", "--out", shQuote(plink_out)))
	if (status != 0) stop("PLINK cis-window extraction failed for ", gene,
		" with exit status ", status, call. = FALSE)
	outputs = paste0(plink_out, c(".bed", ".bim", ".fam"))
	missing_outputs = outputs[!file.exists(outputs) | file.info(outputs)$size == 0]
	if (length(missing_outputs)) stop("PLINK cis-window extraction for ", gene,
		" did not produce: ", paste(missing_outputs, collapse = ", "), call. = FALSE)
	num_completed = num_completed + 1L
	completed_genes[num_completed] = paste0(gene, ".bed")
	file.remove(paste0(plink_out, ".log"))
}

if (num_completed == 0) stop("no autosomal cis-window genotype files were produced", call. = FALSE)
writeLines(sort(unique(completed_genes[seq_len(num_completed)])), manifest)
