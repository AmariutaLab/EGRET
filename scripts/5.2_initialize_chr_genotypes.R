library(data.table)
library(plink2R)
library(optparse)

option_list = list(
        make_option("--output_dir", action="store",default=NA, type='character',
                help="Path to store output"),
        make_option("--tissue", action="store",default=NA, type='character',
                help="tissue for analysis"),
        make_option("--folds", action="store",default=5, type='integer',
                help="Number of cross-validation folds (creates fold_0 through fold_N)"),
        make_option("--genotype_prefix", action="store",default="GTEX_v8_genotypes_pruned", type='character',
                help="Basename of LD-pruned plink genotype files in output_dir/genotype_files/"),
        make_option("--plink_path", action="store",default=NA, type='character',
        	help="path to plink executable")
  )

opt = parse_args(OptionParser(option_list=option_list))

tissue = opt$tissue

bim_file = fread(paste0(opt$output_dir,"/genotype_files/",opt$genotype_prefix,".bim"), header = F, sep = '\t')
number_chr = max(bim_file$V1)

for (fold in 0:opt$folds) {
	for (chr in 1:number_chr) {

                geno_dir = paste0(opt$output_dir,"/transPCO/",tissue,"/fold_",fold,"/genotype_files")
                dir.create(geno_dir, recursive=TRUE, showWarnings=FALSE)

                train_individuals_path = paste0(opt$output_dir,"/fold_",fold,"_info/",tissue,"/train_individuals.txt")
                arg = paste0(opt$plink_path," --bfile ",opt$output_dir,"/genotype_files/",opt$genotype_prefix," --make-bed --keep ",train_individuals_path," --chr ",chr," --out ",geno_dir,"/train_genotypes_chr_",chr)
                system(arg)

                genotypes = read_plink(paste0(geno_dir,"/train_genotypes_chr_",chr),impute='avg')
                genotypes$bed = scale(genotypes$bed)

                rownames = colnames(genotypes$bed)
                genotypes_transposed = t(genotypes$bed)
                rownames(genotypes_transposed) = rownames
                genotypes_transposed_dt = as.data.table(genotypes_transposed, keep.rownames = "SNP")

                fwrite(genotypes_transposed_dt,paste0(geno_dir,"/train_genotypes_chr_",chr,".txt"),row.names = F,col.names = T,quote = F, sep = '\t')

	}


}
