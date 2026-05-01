library("data.table")
library("optparse")

option_list = list(
	make_option("--bfile", action="store",default=NA, type='character',
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


arg = paste0(opt$plink_path, " --bfile ", opt$bfile, " --indep-pairwise ",opt$LD_window, " ",opt$LD_chunk, " ", opt$LD_r2, " --out ",opt$out)
system(arg)

arg = paste0(opt$plink_path, " --bfile ", opt$bfile, " --extract ", opt$out, ".prune.in --maf 0.05 --make-bed --out ", opt$out)
system(arg)

