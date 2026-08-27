library(data.table)
library(optparse)

option_list = list(
  make_option("--gene_info", action="store", default=NA, type='character',
              help="Path to gene info file (e.g. GTEx_V8.txt.gz)"),
  make_option("--crossmap_file", action="store", default=NA, type='character',
              help="Path to cross-mappability strength file (cross_mappability_strength.txt.gz)"),
  make_option("--output_dir", action="store", default=NA, type='character',
              help="Base output directory; BED files written to output_dir/cross_mapped/background_mismatches/"),
  make_option("--cis_window", action="store", default=1000000, type='integer',
              help="Window (bp) around focal gene TSS used to exclude same-region cross-map partners [default: 1000000]"),
  make_option("--crossmap_window", action="store", default=100000, type='integer',
              help="Window (bp) around each cross-mapped gene TSS for BED region bounds [default: 100000]")
)

opt = parse_args(OptionParser(option_list=option_list))

# gene info
all_gene_info = fread(opt$gene_info, header = TRUE)

# cross mappability downloaded from saha et al.
cross_mappable_genes = fread(opt$crossmap_file, header = FALSE)

focal_genes = unique(cross_mappable_genes$V2)

out_dir = file.path(opt$output_dir, "cross_mapped", "background_mismatches_new")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
num_cross_mapped_genes = length(focal_genes)

# Strip "chr" prefix once for all rows (vectorised)
all_gene_info[, chr := sub("^chr", "", `#chrom`)]

# what share of crossmappable genes are mapped to this gene (computed once per focal gene)
cross_mappable_genes[, bg := sum(V3) / num_cross_mapped_genes, by = V2]
kept = cross_mappable_genes[V3 > bg]

# Key both tables so per-gene subsets are O(log n) binary-search lookups
setkey(kept, V2)
setkey(all_gene_info, geneId)

for (gene in focal_genes) {
  print(gene)

  #check if output file already exists
  if (file.exists(file.path(out_dir, paste0(gene, ".bed")))) {next}

  gene_rows = kept[.(gene), nomatch = NULL]
  if (nrow(gene_rows) == 0) next

  print(gene_rows$bg[1])

  focal_gene_info = all_gene_info[.(gene), nomatch = NULL]
  if (nrow(focal_gene_info) == 0) next

  crossmap_gene_info = all_gene_info[.(gene_rows$V1),
                                     .(`#chrom`, chromStart, geneId, chr),
                                     nomatch = NULL]

  crossmap_gene_info = crossmap_gene_info[
    `#chrom` != focal_gene_info$`#chrom` |
    chromStart < focal_gene_info$chromStart - opt$cis_window |
    chromStart > focal_gene_info$chromStart + opt$cis_window, ]

  if (nrow(crossmap_gene_info) == 0) next

  bed_start = pmax(0, crossmap_gene_info$chromStart - opt$crossmap_window)
  bed_end   = crossmap_gene_info$chromStart + opt$crossmap_window
  crossmap_bed = cbind(crossmap_gene_info$chr, bed_start, bed_end, crossmap_gene_info$geneId, "0", "+")
  fwrite(crossmap_bed, file.path(out_dir, paste0(gene, ".bed")), quote = FALSE, col.names = FALSE, row.names = FALSE, sep = '\t')
}
