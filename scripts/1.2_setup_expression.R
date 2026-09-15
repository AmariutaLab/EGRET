# --- 1.2_setup_expression.R ------------------------------------------------
# Subsets expression to the tissue and individual list, then writes it twice: raw,
# and residualized on every covariate column and scaled.
# This is the ONLY place covariates are regressed out. Training stages read the
# regressed file and pass no --covariates; re-adding it would regress twice.
#
#   in : expression, covariates, individuals, gene_info
#   out: expression_files/{tissue}_expression.txt.gz and _expression_regressed.txt.gz; covariate_files/
#   run: 1.0_setup_genotype_and_expression.sh
# --------------------------------------------------------------------------
library("data.table")
library("optparse")

option_list = list(
        make_option("--expression", action="store",default=NA, type='character',
                help="expression file"),
        make_option("--individuals", action="store",default=NA, type='character',
                help="Path to file containing list of individuals to keep"),
	make_option("--covariates", action="store",default=NA, type='character',
                help="Path to file of expresion covariates"),
	make_option("--tissue", action="store",default=NA, type='character',
                help="tissue for analysis"),
        make_option("--gene_info", action="store",default=NA, type='character',
                help="file containg gene info such as gene id, gene name, start, chromosome, and gene type"),
        make_option("--output_dir", action="store",default=NA, type='character',
                help="directory where output files are placed")
  )

opt = parse_args(OptionParser(option_list=option_list))


required = c("expression", "individuals", "covariates", "tissue", "gene_info", "output_dir")
missing = required[!vapply(required, function(o) length(opt[[o]]) == 1L && !is.na(opt[[o]]), logical(1))]
if (length(missing)) stop("missing required argument(s): --", paste(missing, collapse = ", --"), call. = FALSE)

tissue = opt$tissue

expression_file = fread(opt$expression, header = TRUE)

individuals = fread(opt$individuals, header = FALSE, colClasses = "character")
individuals = as.character(unlist(individuals$V1))

expression_file = expression_file[, c("gene_id", individuals), with = FALSE]

#filter out genes which do not have gene info and genes which are not protein coding, lincRNA, or antisense RNA
all_gene_info = fread(opt$gene_info,header = T)
gene_row = match(expression_file$gene_id, all_gene_info$geneId)
expression_file$geneType = all_gene_info$geneType[gene_row]

# Autosomes only. 2.1 already skips X and Y when building cis windows and calls its expected
# output "autosomal", but nothing kept non-autosomal genes out of the matrix in the first place.
# chrM mattered most: those genes are protein_coding, so they passed the type filter, and a
# nuclear genotype panel carries no chrM variants -- so they reached the models with an empty cis
# partition and a trans-only fit. They then dominated the trans results (r2 ~30-80x the
# genome-wide median) and formed co-regulation modules whose cis comparator was zero by
# construction, which 9.5 scores as EGRET beating cis. Dropping them here keeps them out of every
# downstream stage -- trans screens, WGCNA modules, EGRET models, TWAS and the GRN chain -- rather
# than half-present. Accepts either "#chrom" or "chrom", with or without the "chr" prefix.
chrom_col = grep("^#?chrom$", names(all_gene_info), value = TRUE)
if (length(chrom_col) != 1) {
	stop("gene info file needs exactly one chrom column, found: ",
		paste(names(all_gene_info), collapse = ", "), call. = FALSE)
}
expression_file$chrom = sub("^chr", "", as.character(all_gene_info[[chrom_col[1]]][gene_row]))

# An autosome is a numerically-named chromosome. Testing the label rather than whitelisting 1:22
# keeps this correct for any assembly -- mouse has 19 autosomes, dog 38 -- while still dropping
# sex chromosomes (X/Y, or Z/W), the mitochondrion (M/MT), unplaced scaffolds, alt haplotypes,
# and any gene missing from the annotation (chrom is NA, which fails the test).
is_autosome = !is.na(expression_file$chrom) & grepl("^[0-9]+$", expression_file$chrom)
dropped_chroms = sort(unique(expression_file$chrom[!is_autosome & !is.na(expression_file$chrom)]))

n_before = nrow(expression_file)
expression_file = expression_file[expression_file$geneType %in% c('protein_coding', 'lincRNA', 'antisense') &
	is_autosome]
message(sprintf("gene filter: %d of %d genes retained (protein_coding/lincRNA/antisense on autosomes)",
	nrow(expression_file), n_before))
if (length(dropped_chroms)) {
	message("  non-autosomal chromosomes dropped: ", paste(dropped_chroms, collapse = ", "))
}
if (nrow(expression_file) == 0) stop("no genes remain after the gene type / autosome filter", call. = FALSE)
expression_file = expression_file[,-c("geneType","chrom")]

# create directory to write expression
dir.create(paste0(opt$output_dir,"/expression_files/"),recursive = T)

#write expression file to expression folder
fwrite(expression_file, paste0(opt$output_dir,"/expression_files/",tissue,"_expression.txt.gz"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')

# read in covariate file
covariates = fread(opt$covariates, header = TRUE)

# match individuals
covariates = covariates[match(individuals, covariates$ID)]

dir.create(paste0(opt$output_dir,"/covariate_files/"),recursive = T)

fwrite(covariates,paste0(opt$output_dir,"/covariate_files/",tissue,"_covariates.txt.gz"),col.names = T, row.names = F, quote = F, sep = '\t')

expression_matrix = as.matrix(expression_file[, setdiff(names(expression_file), "gene_id"), with = FALSE])
design_matrix = cbind(1, as.matrix(covariates[, 3:ncol(covariates), with = FALSE]))
design_qr = qr(design_matrix)
# qr.resid, NOT `y - X %*% qr.coef(...)`: qr.coef returns NA on a rank-deficient design,
# which silently turns the whole regressed matrix into NA. Reachable in practice -- `sex`
# is constant in the sex-specific GTEx tissues, so it is aliased with the intercept.
residuals_by_gene = qr.resid(design_qr, t(expression_matrix))
residual_gene_expression = t(scale(residuals_by_gene))

residual_gene_expression_df = data.frame(
        gene_id = expression_file$gene_id,
        residual_gene_expression,
        check.names = FALSE
)

# Write to file with fwrite
fwrite(residual_gene_expression_df, paste0(opt$output_dir,"/expression_files/",tissue,"_expression_regressed.txt.gz"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')





