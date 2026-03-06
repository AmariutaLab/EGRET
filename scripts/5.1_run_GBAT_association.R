library('MatrixEQTL')
library('data.table')
library('optparse')

option_list = list(
  make_option("--tissue", action="store",default=NA, type='character',
              help="tissue for analysis"),
  make_option("--output_dir", action="store",default=NA, type='character',
              help="directory where results are stored")
  )

opt = parse_args(OptionParser(option_list=option_list))

if (is.na(opt$tissue)) {
        print("no tissue specified")
        q()
} else {
        tissue = opt$tissue
}

cis_predicted_expression = fread(paste0(opt$output_dir,"/GBAT/",tissue,"/cis_predicted_expression.txt"), header = T)

for (fold in 0:5) {

	target_gene_expression_file_name <- paste0(opt$output_dir,"/fold_", fold, "_info/",tissue,"/train_expression.txt")

    	target_gene_fold <- SlicedData$new()


    	target_gene_fold$fileDelimiter <- "\t"   
    	target_gene_fold$fileOmitCharacters <- "NA" 
    	target_gene_fold$fileSkipRows <- 1    
    	target_gene_fold$fileSkipColumns <- 1  
    	target_gene_fold$fileSliceSize <- 2000     

    	# Load the file
    	target_gene_fold$LoadFile(target_gene_expression_file_name)

        fold_individuals = fread(paste0(opt$output_dir,"/fold_",fold,"_info/",tissue,"/train_individuals.txt"),header = F)
	matched_cols <- match(paste0("0:", unlist(fold_individuals$V1)), colnames(cis_predicted_expression))

	# Remove NA values from matched indices
	matched_cols <- c(1,matched_cols[!is.na(matched_cols)])
	
	prediction_subset <- cis_predicted_expression[, matched_cols, with = FALSE]
	
        if (!file.exists(paste0("GBAT/",tissue,"/fold_",fold,"/"))) {
                dir.create(paste0("GBAT/",tissue,"/fold_",fold,"/"), recursive = T)
        }
        fwrite(prediction_subset,paste0(opt$output_dir,"/GBAT/",tissue,"/fold_",fold,"/cis_predicted_expression.txt"),row.names = F, col.names = T,quote = F, sep = '\t')

        predicted_gene_expression_file_name = paste0(opt$output_dir,"/GBAT/",tissue,"/fold_",fold,"/cis_predicted_expression.txt")
        predicted_gene_expression = SlicedData$new()
        predicted_gene_expression$fileDelimiter = "\t"     
        predicted_gene_expression$fileOmitCharacters = "NA" 
        predicted_gene_expression$fileSkipRows = 1         
        predicted_gene_expression$fileSkipColumns = 1      
        predicted_gene_expression$fileSliceSize = 2000     
        predicted_gene_expression$LoadFile( predicted_gene_expression_file_name )


        errorCovariance = numeric()


        me = Matrix_eQTL_main(
                gene = target_gene_fold,
                snps = predicted_gene_expression,
                output_file_name = paste0(opt$output_dir,"/GBAT/",tissue,"/fold_",fold,"/association_results.txt"),
                pvOutputThreshold = 0.0001,
                useModel = modelLINEAR,
                errorCovariance = errorCovariance,
                verbose = FALSE,
                pvalue.hist = FALSE,
                min.pv.by.genesnp = FALSE,
                noFDRsaveMemory = FALSE)

}
