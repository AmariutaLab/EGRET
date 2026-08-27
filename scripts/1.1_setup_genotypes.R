# --- 1.1_setup_genotypes.R -------------------------------------------------
# LD-prunes the genotype panel with plink so downstream association tests run on
# an approximately independent SNP set.
#
#   in : genotype_prefix .bed/.bim/.fam
#   out: genotype_files/{genotype_output_prefix}.bed/.bim/.fam
#   run: 1.0_setup_genotype_and_expression.sh
# --------------------------------------------------------------------------
library("data.table")
library("optparse")

option_list = list(
	make_option("--genotype_prefix", action="store",default=NA, type='character',
                help="Path to bed, bim, fam files"),
	make_option("--LD_r2", action="store",default=0.9, type='numeric',
                help="LD r2 threshold used for pruning snps"),
	make_option("--LD_window", action="store",default=100, type='numeric',
                help="LD window size used for pruning"),
	make_option("--LD_chunk", action="store",default=1, type='numeric',
                help="LD chunk size used for pruning"),	   
        make_option("--out", action="store",default=NA, type='character',
                help="Path to store output"),
	make_option("--plink_path", action="store",default=NA, type='character',
                help="Path to plink")
  )

opt = parse_args(OptionParser(option_list=option_list))



required = c("genotype_prefix", "out", "plink_path")
missing = required[!vapply(required, function(o) length(opt[[o]]) == 1L && !is.na(opt[[o]]), logical(1))]
if (length(missing)) stop("missing required argument(s): --", paste(missing, collapse = ", --"), call. = FALSE)

run_plink = function(args, label) {
	status = system2(opt$plink_path, args = args)
	if (status != 0) stop(label, " failed with exit status ", status, call. = FALSE)
}
require_outputs = function(paths, label) {
	missing = paths[!file.exists(paths) | file.info(paths)$size == 0]
	if (length(missing)) stop(label, " did not produce: ", paste(missing, collapse = ", "), call. = FALSE)
}

run_plink(c("--bfile", shQuote(opt$genotype_prefix), "--indep-pairwise",
	opt$LD_window, opt$LD_chunk, opt$LD_r2, "--out", shQuote(opt$out)), "PLINK LD pruning")
require_outputs(paste0(opt$out, ".prune.in"), "PLINK LD pruning")

run_plink(c("--bfile", shQuote(opt$genotype_prefix), "--extract", shQuote(paste0(opt$out, ".prune.in")),
	"--maf", "0.05", "--make-bed", "--out", shQuote(opt$out)), "PLINK genotype pruning")
require_outputs(paste0(opt$out, c(".bed", ".bim", ".fam")), "PLINK genotype pruning")
