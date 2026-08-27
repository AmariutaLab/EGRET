# --- 5.2_initialize_chr_genotypes.R ----------------------------------------
# Splits the fold's genotype matrix by chromosome so the per-module association jobs
# can stream one chromosome at a time instead of the whole panel.
#
#   in : fold_{i}_info/{tissue}/train_genotypes.txt
#   out: transPCO/{tissue}/fold_{i}/genotype_files/train_genotypes_chr_{c}.txt
#   run: 5.0_run_transPCO_stage.sh
# --------------------------------------------------------------------------
library(data.table)
library(plink2R)
library(optparse)

option_list = list(
        make_option("--output_dir", action="store",default=NA, type='character',
                help="Path to store output"),
        make_option("--tissue", action="store",default=NA, type='character',
                help="tissue for analysis"),
        make_option("--folds", action="store",default=NA, type='integer',
                help="Number of cross-validation folds (creates fold_0 through fold_N)"),
        make_option("--genotype_prefix", action="store",default=NA, type='character',
                help="Basename of LD-pruned plink genotype files in output_dir/genotype_files/"),
        make_option("--plink_path", action="store",default=NA, type='character',
        	help="path to plink executable")
  )

opt = parse_args(OptionParser(option_list=option_list))


required = c("output_dir", "tissue", "folds", "genotype_prefix", "plink_path")
missing = required[!vapply(required, function(o) length(opt[[o]]) == 1L && !is.na(opt[[o]]), logical(1))]
if (length(missing)) stop("missing required argument(s): --", paste(missing, collapse = ", --"), call. = FALSE)

tissue = opt$tissue

bim_file = fread(paste0(opt$output_dir,"/genotype_files/",opt$genotype_prefix,".bim"), header = F, sep = '\t')
number_chr = max(bim_file$V1)

for (fold in 0:opt$folds) {
	for (chr in 1:number_chr) {

                geno_dir = paste0(opt$output_dir,"/transPCO/",tissue,"/fold_",fold,"/genotype_files")
                dir.create(geno_dir, recursive=TRUE, showWarnings=FALSE)

                train_individuals_path = paste0(opt$output_dir,"/fold_",fold,"_info/",tissue,"/train_individuals.txt")
		plink_out = paste0(geno_dir, "/train_genotypes_chr_", chr)
		status = system2(opt$plink_path, args = c("--bfile",
			shQuote(paste0(opt$output_dir, "/genotype_files/", opt$genotype_prefix)),
			"--make-bed", "--keep", shQuote(train_individuals_path), "--chr", chr,
			"--out", shQuote(plink_out)))
		if (status != 0) {
			stop(sprintf("PLINK chromosome extraction failed for fold_%s chr%s with exit status %s",
			             fold, chr, status))
		}
		expected = paste0(plink_out, c(".bed", ".bim", ".fam"))
		if (any(!file.exists(expected)) || any(file.info(expected)$size == 0)) {
			stop(sprintf("PLINK chromosome extraction did not create a complete fileset for fold_%s chr%s",
			             fold, chr))
		}

		genotypes = read_plink(plink_out,impute='avg')
                genotypes$bed = scale(genotypes$bed)

                rownames = colnames(genotypes$bed)
                genotypes_transposed = t(genotypes$bed)
                rownames(genotypes_transposed) = rownames
                genotypes_transposed_dt = as.data.table(genotypes_transposed, keep.rownames = "SNP")

                fwrite(genotypes_transposed_dt,paste0(geno_dir,"/train_genotypes_chr_",chr,".txt"),row.names = F,col.names = T,quote = F, sep = '\t')

	}


}
