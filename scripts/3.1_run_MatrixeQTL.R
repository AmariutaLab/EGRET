# --- 3.1_run_MatrixeQTL.R --------------------------------------------------
# Univariate trans-eQTL scan: tests every SNP against every gene, per fold, keeping
# associations below a permissive p-value threshold for later FDR control.
#
#   in : fold_{i}_info/{tissue}/train_expression.txt and train_genotypes.txt
#   out: MatrixeQTL/{tissue}/association_results_fold_{i}_0.0001_threshold.txt
#   run: 3.0_run_MatrixeQTL_stage.sh
# --------------------------------------------------------------------------
library('MatrixEQTL')
library('optparse')

option_list = list(
  make_option("--tissue", action="store",default=NA, type='character',
              help="tissue for analysis"),
  make_option("--folds", action="store",default=NA, type='numeric',
              help="number of folds for analysis"),
  make_option("--output_dir", action="store",default=NA, type='character',
              help="directory where results are stored")
)

opt = parse_args(OptionParser(option_list=option_list))


required = c("tissue", "folds", "output_dir")
missing = required[!vapply(required, function(o) length(opt[[o]]) == 1L && !is.na(opt[[o]]), logical(1))]
if (length(missing)) stop("missing required argument(s): --", paste(missing, collapse = ", --"), call. = FALSE)

tissue = opt$tissue
output_dir = opt$output_dir

for (fold in 0:opt$folds) {

	expression_file_name = paste0(output_dir,"/fold_",fold,"_info/",tissue,"/train_expression.txt")
        gene = SlicedData$new()
        gene$fileDelimiter = "\t"    
        gene$fileOmitCharacters = "NA" 
        gene$fileSkipRows = 1         
        gene$fileSkipColumns = 1     
        gene$fileSliceSize = 2000     
        gene$LoadFile( expression_file_name )

	SNP_file_name = paste0(output_dir,"/fold_",fold,"_info/",tissue,"/train_genotypes.txt")
	snps = SlicedData$new()
	snps$fileDelimiter = "\t"   
	snps$fileOmitCharacters = "NA" 
	snps$fileSkipRows = 1        
	snps$fileSkipColumns = 1    
	snps$fileSliceSize = 20000    
	snps$LoadFile( SNP_file_name )

	errorCovariance = numeric()

	dir.create(paste0(output_dir,"/MatrixeQTL/",tissue,"/"),recursive = T)
	me = Matrix_eQTL_main(
                        gene = gene,
                        snps = snps,
                        output_file_name = paste0(output_dir,"/MatrixeQTL/",tissue,"/association_results_fold_",fold,"_0.0001_threshold.txt"),
                        pvOutputThreshold = 0.0001, 
                        useModel = modelLINEAR,
                        errorCovariance = errorCovariance,
                        verbose = TRUE,
                        pvalue.hist = TRUE,
                        min.pv.by.genesnp = FALSE,
                        noFDRsaveMemory = FALSE)
	
}


