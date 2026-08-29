# --- 6.2_make_EGRET_beds_by_fold.R -----------------------------------------
# Unions each gene's cis window with the trans SNPs found by MatrixeQTL, GBAT and
# transPCO, subtracts cross-mappable regions, and writes the variant set EGRET trains on.
# Also writes the z-score matrix that xtune uses as its prior.
# Hard error when exclude_crossmap=TRUE but no cross-mappable beds exist.
#
#   in : {matrixeqtl,gbat,transpco}_bed_dir trees; cross_mapped/{crossmap_dir}/
#   out: bed_files/{tissue}/{egret_subdir}/fold_{i}/{gene}; z_matrices/{tissue}/{egret_subdir}/
#   run: 6.1_run_make_EGRET_beds.sh
# --------------------------------------------------------------------------
library('data.table')
library('optparse')

option_list = list(
  make_option("--tissue", action="store", default=NA, type='character',
              help="tissue"),
  make_option("--output_dir", action="store", default=NA, type='character',
              help="base directory prefixed to all file paths"),
  make_option("--matrixeqtl_bed_dir", action="store", default=NA, type='character',
              help="the dir containing MatrixeQTL bed files"),
  make_option("--gbat_bed_dir", action="store", default=NA, type='character',
              help="the dir containing GBAT bed files"),
  make_option("--transpco_bed_dir", action="store", default=NA, type='character',
              help="the dir which contains transPCO bed files"),
  make_option("--exclude_crossmap", action="store", default=TRUE, type='logical',
              help="if True, will exclude regions which are cross mappable to that gene"),
  make_option("--crossmap_dir", action="store", default=NA, type='character',
              help="subdirectory under bed_files/cross_mapped/ for cross-mappable regions"),
  make_option("--egret_subdir", action = 'store', default=NA, type = 'character',
	          help = 'directory to store bed files, plink results, and z matrices'),
  make_option("--bed_models", action = 'store', default=NA, type = 'character',
	          help = 'Cis, MatrixeQTL, GBAT, transPCO'),
  make_option("--fold", action = 'store', default=NA, type = 'numeric',
              help = 'should be 0:num_folds'),
  make_option("--genotype_prefix", action="store", default=NA,type='character',
              help="prefix for genotype files (.bim/.bed/.fam)"),
  make_option("--plink_path", action="store", default=NA, type='character',
              help="path to plink2 executable"),
  make_option("--expression", action="store", default=NA, type='character',
              help="optional expression file to define gene list; if not provided, uses cis bed dir"),
  make_option("--gene_info", action="store", default=NA, type='character',
              help="Path to file containing gene information"),
  make_option("--cis_window_size", action="store", default=1000000, type='numeric',
              help="size of cis window to create"),
  make_option("--cleanup_bed_lists", action="store", default=FALSE, type='logical',
              help="delete each per-gene SNP list once plink has consumed it")
  )

opt = parse_args(OptionParser(option_list=option_list))

required = c("tissue", "output_dir", "matrixeqtl_bed_dir", "gbat_bed_dir", "transpco_bed_dir", "crossmap_dir", "egret_subdir", "bed_models", "fold", "genotype_prefix", "plink_path", "expression", "gene_info")
missing = required[!vapply(required, function(o) length(opt[[o]]) == 1L && !is.na(opt[[o]]), logical(1))]
if (length(missing)) stop("missing required argument(s): --", paste(missing, collapse = ", --"), call. = FALSE)

tissue = opt$tissue
fold = opt$fold
base_dir = opt$output_dir

expression = fread(opt$expression, header = T)
genes = expression[[1]]
all_gene_info = fread(opt$gene_info, header = T)

allsnps = fread(paste0(base_dir, "/genotype_files/", opt$genotype_prefix, ".bim"), header = F)

output_dir = opt$egret_subdir

models = unique( c(unlist(strsplit(opt$bed_models,','))) )

crossmap_dir = opt$crossmap_dir

# exclude_crossmap=TRUE with no cross_mapped/ tree would skip the filter for every gene and
# still exit 0, producing a complete-looking set of models that are not the published
# method. The beds are not bundled, so a fresh output_dir hits this. Checked at directory
# level; a single gene with no bed of its own is a normal skip, handled per-gene below.
crossmap_bed_dir = paste0(base_dir, "/cross_mapped/", crossmap_dir)
if (opt$exclude_crossmap) {
	n_crossmap_beds = if (dir.exists(crossmap_bed_dir)) {
		length(list.files(crossmap_bed_dir, pattern = "\\.bed$"))
	} else {
		0L
	}
	if (n_crossmap_beds == 0) {
		stop("exclude_crossmap is TRUE but no cross-mappability BED files were found in ",
		     crossmap_bed_dir, ".\n",
		     "  Generate them first:\n",
		     "    Rscript 0.1_make_crossmappable_gene_beds.R --gene_info <gene_info> \\\n",
		     "        --crossmap_file <cross_mappability.tsv.gz> --output_dir ", base_dir,
		     " \\\n        --crossmap_dir ", crossmap_dir, "\n",
		     "  ...or set exclude_crossmap=FALSE to run without cross-mappability filtering.")
	}
	cat("cross-mappability exclusion enabled:", n_crossmap_beds, "BED files in",
	    crossmap_bed_dir, "\n")
}

trans_sources = c("MatrixeQTL", "GBAT", "transPCO")
source_stats = data.table(source=trans_sources, files=0L, listed_snps=0L,
                          matched_snps=0L, contributing_genes=0L)

read_source = function(path, source_name) {
	row = match(source_name, source_stats$source)
	if (!file.exists(path)) return(character())
	source_stats[row, files := files + 1L]
	if (file.info(path)$size == 0) return(character())
	tab = tryCatch(fread(path, header=FALSE), error=function(e) data.table())
	if (ncol(tab) < 1 || nrow(tab) < 1) return(character())
	ids = unique(as.character(tab[[1]]))
	ids = ids[!is.na(ids) & nzchar(ids)]
	source_stats[row, listed_snps := listed_snps + length(ids)]
	ids
}

record_final_matches = function(ids, source_name, final_ids) {
	row = match(source_name, source_stats$source)
	matched = intersect(unique(ids), unique(as.character(final_ids)))
	source_stats[row, matched_snps := matched_snps + length(matched)]
	if (length(matched) > 0) {
		source_stats[row, contributing_genes := contributing_genes + 1L]
	}
}

plink_left_nothing = function(log_path) {
	if (!file.exists(log_path)) return(FALSE)
	any(grepl("No variants remaining after main filters",
		readLines(log_path, warn = FALSE), fixed = TRUE))
}
num_empty = 0L

matched_bed = function(ids) {
	idx = match(ids, allsnps$V2, nomatch=0L)
	idx = idx[idx > 0L]
	if (length(idx) == 0) return(matrix(ncol=6, nrow=0))
	hits = allsnps[idx]
	cbind(hits$V1, hits$V4, hits$V4 + 1, hits$V2, "0", "+")
}

for (gene in genes) {
	bed_file = matrix(ncol = 6,nrow = 0)

	gene_name = gene

	gene_info = all_gene_info[all_gene_info$geneId == gene_name, ]
	chr = strsplit(gene_info$'#chrom', split = 'chr')[[1]][2]
	tss = gene_info$chromStart
	if (chr == 'X' | chr == 'Y') {
		next
	}

	# 2.1 writes a six-column window here. Reject anything else rather than rbind-ing a
	# one-column SNP list into the bed at the wrong width.
	cis_bed_path = paste0(base_dir, "/bed_files/",tissue,"/cis/",gene,".bed")
	cis_bed = NULL
	if (file.exists(cis_bed_path)) {
		candidate = fread(cis_bed_path,header = F)
		if (ncol(candidate) == 6) {
			cis_bed = candidate
		} else {
			warning(sprintf("%s has %d columns, expected 6; recomputing the cis window", cis_bed_path, ncol(candidate)))
		}
	}
	if (is.null(cis_bed)) {
		bounds_u = max(1, tss - floor(opt$cis_window_size/2))
		bounds_l = tss + floor(opt$cis_window_size/2)
		bed = c(chr, bounds_u, bounds_l, gene_name, 0, "+")
		cis_bed = matrix(bed, nrow = 1, ncol = 6)
	}

	#GBAT
	GBAT_bed_path = paste0(base_dir, "/GBAT/",tissue,"/fold_",fold,"/",opt$gbat_bed_dir,"/",gene_name,".txt")
	GBAT_ids = read_source(GBAT_bed_path, "GBAT")
	# Matrix eQTL
	MatrixeQTL_bed_path = paste0(base_dir, "/MatrixeQTL/",tissue,"/fold_",fold,"/",opt$matrixeqtl_bed_dir,"/",gene_name,".txt")
	MatrixeQTL_ids = read_source(MatrixeQTL_bed_path, "MatrixeQTL")

	# trans-PCO
	transPCO_bed_path = paste0(base_dir, "/transPCO/",tissue,"/fold_",fold,"/",opt$transpco_bed_dir,"/",gene_name,".txt")
	transPCO_ids = read_source(transPCO_bed_path, "transPCO")
	GBAT_bed = matched_bed(GBAT_ids)
	MatrixeQTL_bed = matched_bed(MatrixeQTL_ids)
	transPCO_bed = matched_bed(transPCO_ids)

	if ("Cis" %in% models) { bed_file = rbind(bed_file,cis_bed) }
	if ("MatrixeQTL" %in% models ) { bed_file = rbind(bed_file,MatrixeQTL_bed) }
	if ("GBAT" %in% models) { bed_file = rbind(bed_file,GBAT_bed) }
	if ("transPCO" %in% models) { bed_file = rbind(bed_file,transPCO_bed) }

	bed_file_dir = paste0(base_dir, "/bed_files/",tissue,"/",output_dir,"/fold_",fold,"/")
	dir.create(bed_file_dir,recursive = T,showWarnings = FALSE)
	fwrite(bed_file,paste0(bed_file_dir,gene),row.names = F , col.names = F,quote = F, sep = '\t')

	plink_output_dir = paste0(base_dir, "/plink_results/",tissue,"/",output_dir,"/fold_",fold,"/")
	dir.create(plink_output_dir,recursive = T,showWarnings = FALSE)

	crossmap_bed_path = paste0(base_dir, "/cross_mapped/",crossmap_dir,"/",gene,".bed")
	plink_args = c("--bfile", shQuote(paste0(base_dir, "/genotype_files/", opt$genotype_prefix)),
		"--extract", "bed1", shQuote(paste0(bed_file_dir, gene)), "--make-bed")
	if(file.exists(crossmap_bed_path) & opt$exclude_crossmap) {
		print("exluding cross mappable genes")
		plink_args = c(plink_args, "--exclude", "bed1", shQuote(crossmap_bed_path))
	} else {
		print("no cross mappable genes to exlude")
	}
	plink_out = paste0(plink_output_dir, gene_name)
	status = system2(opt$plink_path, args = c(plink_args, "--out", shQuote(plink_out)))
	if (status != 0) {
		if (plink_left_nothing(paste0(plink_out, ".log"))) {
			message(sprintf("Skipping %s fold_%s: no SNPs remained after extraction", gene_name, fold))
			num_empty = num_empty + 1L
			if (opt$cleanup_bed_lists) { unlink(paste0(bed_file_dir,gene)) }
			next
		}
		stop(sprintf("PLINK EGRET extraction failed for %s fold_%s with exit status %s",
		             gene_name, fold, status))
	}
	expected = paste0(plink_out, c(".bed", ".bim", ".fam"))
	if (any(!file.exists(expected)) || any(file.info(expected)$size == 0)) {
		stop(sprintf("PLINK EGRET extraction did not create a complete fileset for %s fold_%s",
		             gene_name, fold))
	}

	# Cleanup: this per-gene SNP list exists only for the `plink --extract` above.
	if (opt$cleanup_bed_lists) { unlink(paste0(bed_file_dir,gene)) }

	plink_results = fread(paste0(plink_out, ".bim"), header = FALSE)

	# Count after extraction and cross-map exclusion, not against the input BIM: a source
	# whose SNPs are all removed here has an all-zero z-matrix column and contributes nothing.
	record_final_matches(MatrixeQTL_ids, "MatrixeQTL", plink_results$V2)
	record_final_matches(GBAT_ids, "GBAT", plink_results$V2)
	record_final_matches(transPCO_ids, "transPCO", plink_results$V2)

	MatrixeQTL_matches = match(MatrixeQTL_ids,plink_results$V2)
	MatrixeQTL_matches = MatrixeQTL_matches[!is.na(MatrixeQTL_matches)]
	plink_results$MatrixeQTL = 0
	plink_results$MatrixeQTL[MatrixeQTL_matches] = 1

	GBAT_matches = match(GBAT_ids,plink_results$V2)
	GBAT_matches = GBAT_matches[!is.na(GBAT_matches)]
	plink_results$GBAT = 0
	plink_results$GBAT[GBAT_matches] = 1

	transPCO_matches = match(transPCO_ids,plink_results$V2)
	transPCO_matches = transPCO_matches[!is.na(transPCO_matches)]
	plink_results$transPCO = 0
	plink_results$transPCO[transPCO_matches] = 1

	z_matrix = data.table()
	if ("MatrixeQTL" %in% models ) { z_matrix = cbind(z_matrix,plink_results$MatrixeQTL) }
	if ("GBAT" %in% models) { z_matrix = cbind(z_matrix,plink_results$GBAT) }
	if ("transPCO" %in% models) { z_matrix = cbind(z_matrix,plink_results$transPCO) }

	z_matrix_dir = paste0(base_dir, "/z_matrices/",tissue,"/",output_dir,"/fold_",fold,"/")
	dir.create(z_matrix_dir, recursive = T,showWarnings = FALSE)
	fwrite(z_matrix,paste0(z_matrix_dir,gene_name,"_z_matrix.txt"),sep = '\t',col.names = F,row.names = F,quote=F)
}

message(sprintf("fold %s: %d gene(s) skipped for an empty SNP set", fold, num_empty))

missing_sources = character(0)
for (src in trans_sources) {
	if (!(src %in% models)) next
	s = source_stats[source == src]
	cat(sprintf("fold %s: %s files=%d listed_snps=%d matched_snps=%d contributing_genes=%d\n",
	            fold, src, s$files, s$listed_snps, s$matched_snps, s$contributing_genes))
	if (s$matched_snps == 0L) missing_sources = c(missing_sources, src)
}
if (length(missing_sources) > 0) {
	stop(sprintf(paste0("fold %s: requested source(s) contributed zero SNPs after final PLINK/cross-map filtering: %s. ",
	                    "Check upstream outputs, source directories, and variant identifiers."),
	             fold, paste(missing_sources, collapse = ", ")))
}
