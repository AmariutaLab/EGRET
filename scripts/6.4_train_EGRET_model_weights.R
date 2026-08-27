# --- 6.4_train_EGRET_model_weights.R ---------------------------------------
# Trains one gene's model. Fits three variant partitions -- gw (cis+trans), cis only,
# trans only -- across lasso/enet/blup/xtune, and deploys whichever partition wins on
# cross-validated R2. That selection is what stops null trans features from hurting:
# when trans adds nothing, the cis partition wins and no accuracy is lost.
# Reads the pre-regressed expression, so --covariates is deliberately not passed.
#
#   in : working/{tissue}/{egret_subdir}/fold_{i}/{gene}.*; z_matrices/
#   out: EGRET_models/{tissue}/{egret_subdir}/{gene}.wgt.RDat (wgt.matrix, cv.performance, best_model)
#   run: 6.3_train_EGRET_chunk.sh
# --------------------------------------------------------------------------
library('data.table')
library('optparse')
library('glmnet')
library('xtune')
library('plink2R')

options(scipen = 999)

option_list = list(
  make_option("--gene", action="store", default=NA, type='character',
	help="gene name"),
  make_option("--working_dir", action="store", default=NA, type='character',
        help="location where bed/bim/fam files are stored"),
  make_option("--fusion_models", action="store", default=NA, type='character',
        help="Comma-separated list of prediction models"),
  make_option("--weights_dir",action ='store',default=NA,type='character',
	help="path to weights dir"),
  make_option("--tissue",action = 'store',default=NA,type='character',
        help="the tissue for analysis"),
  make_option("--plink_path",action = 'store',default=NA,type='character',
        help="path to plink"),
  make_option("--gemma_path",action = 'store',default=NA,type='character',
        help="path to gemma"),
  make_option("--covariates",action = 'store',default=NA,type='character',
        help="covariates"),
  make_option("--z_matrix_dir", action = 'store', default = NA, type = 'character',
	help = 'directory containing z_matrix for each fold'),
  make_option("--clean", action="store", default=FALSE,
        help="whether to clean working dir after run"),
  make_option("--gene_info", action="store", default=NA, type='character',
        help="path to gene info file (e.g. GTEx_V8.txt.gz)"),
  make_option("--folds", action="store", default=NA, type='numeric',
        help="number of cross-validation folds"),
  make_option("--output_dir", action="store", default=NA, type='character',
        help="Path to output dir"))

opt = parse_args(OptionParser(option_list=option_list))

required = c("gene", "working_dir", "fusion_models", "weights_dir", "tissue", "plink_path", "gemma_path", "z_matrix_dir", "gene_info", "folds", "output_dir")
missing = required[!vapply(required, function(o) length(opt[[o]]) == 1L && !is.na(opt[[o]]), logical(1))]
if (length(missing)) stop("missing required argument(s): --", paste(missing, collapse = ", --"), call. = FALSE)

tissue = opt$tissue

models = unique( c(unlist(strsplit(opt$fusion_models,','))) )
M = length(models)

run_plink = function(args, out_prefix, context) {
	status = system2(opt$plink_path, args = args)
	if (status != 0) {
		stop(sprintf("PLINK %s failed with exit status %s", context, status))
	}
	expected = paste0(out_prefix, c(".bed", ".bim", ".fam"))
	if (any(!file.exists(expected)) || any(file.info(expected)$size == 0)) {
		stop(sprintf("PLINK %s did not create a complete fileset", context))
	}
}

# ---- GLMNET LASSO
weights.enet = function( genos , pheno , alpha=1 ) {   #alpha has been changed to 1
        eff.wgt = matrix( 0 , ncol=1 , nrow=ncol(genos) )
        # remove monomorphics
        sds = apply( genos  , 2 , sd )
        keep = sds != 0 & !is.na(sds)
        enet = cv.glmnet( x=genos[,keep] , y=pheno , alpha=alpha , nfold=5 , intercept=T , standardize=F )
        eff.wgt[ keep ] = coef( enet , s = "lambda.min")[2:(sum(keep)+1)]
        return( eff.wgt )
}

#BLUP/BSLMM
weights.bslmm = function( input , bv_type , snp , fold , tag="" ) {
        outdir = paste0(opt$working_dir, "/fold_", fold)
        # tag keeps the three partitions from overwriting each other's gemma output
        prefix = paste0(opt$gene, ".BSLMM", tag)
        arg = paste( opt$gemma_path , " -miss 1 -maf 0 -r2 1 -rpace 1000 -wpace 1000 -bfile " , input , " -bslmm " , bv_type , " -outdir " , outdir , " -o " , prefix , sep='' )
	system(arg)
        out = paste0(outdir, "/", prefix)
	print(file.exists(paste(out,".param.txt",sep='')))
        eff = read.table( paste(out,".param.txt",sep=''),head=T,as.is=T)
        eff.wgt = rep(NA,length(snp))
        m = match( snp , eff$rs )
        m.keep = !is.na(m)
        m = m[m.keep]
        eff.wgt[m.keep] = (eff$alpha + eff$beta * eff$gamma)[m]
        return( eff.wgt )
}


# ---- XTUNE LASSO ----
weights.xtune = function(genos, pheno, alpha=1, fold, idx=NULL) {
    eff.wgt = matrix(0, ncol=1, nrow=ncol(genos))
    # remove monomorphics
    sds = apply(genos, 2, sd)
    keep = sds != 0 & !is.na(sds)
    print(paste0("The current fold is ",fold))

    #system(z_matrix_args)
    z_matrix_file = paste0(opt$z_matrix_dir, "/fold_", fold, "/", opt$gene, "_z_matrix.txt")
    if (file.exists(z_matrix_file)) {
	Z = read.table(paste0(opt$z_matrix_dir, "/fold_", fold, "/", opt$gene, "_z_matrix.txt"), header = FALSE)
	# z-matrix rows are one-per-bim-row in bim order, so the partition is a row subset
	if (!is.null(idx)) {
		Z = Z[idx, , drop=FALSE]
		keepZ = apply(Z, 2, function(x) length(unique(x)) > 1)
		# an all-constant Z is guaranteed to error out into the no-Z fallback below
		if (!any(keepZ)) Z = matrix() else Z = Z[, keepZ, drop=FALSE]
	}
    }
    else {
	    Z = matrix()
    }
    alphas = NULL
    xtune_success = FALSE
    print("Trying with z matrix")
    fit.xtune = tryCatch({
        result = xtune(genos, pheno, Z, c = alpha, epsilon = 15, sigma.square = 1, family = 'linear')
        alphas = result$alpha.est
        eff.wgt = result$beta.est
        eff.wgt = eff.wgt[-1, ]  # this returns a weight for the intercept which needs to be removed
        xtune_success = TRUE
    }, error = function(e) {
        message("Was not able to run xtune with z-matrix")
        NULL
    })

    if (!xtune_success) {
        fit.xtune = tryCatch({
            result = xtune(genos, pheno, c = alpha, epsilon = 15, sigma.square = 1, family = 'linear')
            eff.wgt = result$beta.est
            eff.wgt = eff.wgt[-1, ]
        }, error = function(e) {
            message("Was not able to run xtune without z-matrix")
            eff.wgt = rep(NA,ncol(genos))
            print(eff.wgt)
        })
    }

    return(list(weights = eff.wgt, alphas = alphas))
}

# Marginal Z-scores (used for top1)
weights.marginal = function( genos , pheno , beta=F ) {
        if ( beta ) eff.wgt = t( genos ) %*% (pheno) / ( nrow(pheno) - 1)
        else eff.wgt = t( genos ) %*% (pheno) / sqrt( nrow(pheno) - 1 )
        return( eff.wgt )
}


#clean up function
cleanup = function() {
        if (opt$clean ) {
                arg = paste0("rm -rf ", opt$working_dir, "/*")
                system(arg)
        }
}

num_folds = opt$folds
fold_individuals = list()
for (f in 0:num_folds) {
    fold_individuals[[paste0("fold_", f)]] = fread(
        paste0(opt$output_dir,"/fold_", f, "_info/", tissue, "/train_individuals.txt"), header = F, colClasses = "character"
    )
}

if(!is.na(opt$covariates)) {
	covar = fread(opt$covariates,header = T)
} else {
	covar = NA
}

# each partition is fit separately, so performance and held-out calls are tracked per partition
partitions = c("gw","cis","trans")
cv.performance = list()
cv.calls = list()
for (p in partitions) {
	cv.performance[[p]] = matrix(NA,nrow=2,ncol=M,dimnames=list(c("rsq","pval"),models))
	cv.calls[[p]] = matrix(NA,nrow=nrow(fold_individuals[["fold_0"]]),ncol=M)
}
alphas = c()

#new stuff to do cis/trans analysis
all_gene_info = fread(opt$gene_info, header = T)
gene_info = all_gene_info[all_gene_info$geneId == opt$gene,]
if (nrow(gene_info) != 1) {
    stop(sprintf("Expected exactly 1 row for gene %s in gene_info_file; found %d", opt$gene, nrow(gene_info)))
}
u_bound = max(0,gene_info$chromStart - 2500000)
l_bound = gene_info$chromStart + 2500000
chr = strsplit(gene_info$'#chrom','chr')[[1]][2]

#stays TRUE only if no fold has a single trans SNP, in which case gw and cis are the same model
gw_is_cis = TRUE

for (fold in 1:num_folds) {
	print("starting cross val")
	# read in genotypes
	genos = read_plink(paste0(opt$working_dir,"/fold_",fold,"/",opt$gene),impute="avg")

	# important : genotypes are standardized and scaled here:
	genos$bed = scale(genos$bed)
	pheno = genos$fam[,c(1,2,6)]

	#regress out covariates
	if (!is.na(opt$covariates)) {
		m = match( unlist(pheno$V2) ,unlist(covar$ID) )
                m.keep = !is.na(m)
                pheno = pheno[m.keep,]
                m = m[m.keep]
                covar_fold = covar[m,]
                reg = summary(lm(as.matrix(pheno[,3]) ~ as.matrix(covar_fold[,3:ncol(covar_fold)]) ))
                pheno[,3] = scale(reg$resid)
                pheno[,3] = scale(pheno[,3])   #is this necessary?? Should we scale just the train portion separately

        }

	#new stuff for cis/trans analysis
	bim = genos$bim
	cis = (bim$V1 == chr & bim$V4 > u_bound & bim$V4 < l_bound)
        cis_snps = bim[cis,]
	trans_snps = bim[!cis,]

	#checking for missing snps
	nasnps = apply( is.na(genos$bed) , 2 , sum )
	if ( sum(nasnps) != 0 ) {
	cat( "WARNING :",sum(nasnps != 0),"SNPs could not be scaled and were zeroed out, make sure all SNPs are polymorphic\n" , file=stderr())
	genos$bed[,nasnps != 0] = 0
	}

	#subset geno and pheno files to train individuals
	keep_geno = match(paste0("0:",unlist(fold_individuals[[paste0("fold_",fold)]])),rownames(genos$bed))
	train_geno = genos$bed[keep_geno,]

	keep_pheno = match(unlist(fold_individuals[[paste0("fold_",fold)]]),pheno[,2])
        train_pheno = pheno[keep_pheno,]


	#call models
	N = nrow(pheno)
	#cv.calls = matrix(NA,nrow=N,ncol=M)
	pheno_file = paste0(opt$working_dir,"/fold_",fold,"/",opt$gene)
	cv.file = paste(opt$working_dir,"/fold_",fold,"/",opt$gene,".cv",sep='')
	write.table( train_pheno , quote=F , row.names=F , col.names=F , file=paste(cv.file,".keep",sep=''))
	run_plink(c("--allow-no-sex", "--bfile", shQuote(pheno_file), "--keep",
	            shQuote(paste0(cv.file, ".keep")), "--out", shQuote(cv.file), "--make-bed"),
	          cv.file, sprintf("cross-validation extraction for %s fold_%s", opt$gene, fold))

	#each partition gets its own fit rather than a masked copy of the genome-wide one
	part_idx = list( gw = seq_len(ncol(genos$bed)) , cis = which(cis) , trans = which(!cis) )

	#when one side is empty the other side IS the genome-wide design matrix, so fit it
	#once and copy rather than fitting the identical problem twice
	do_parts = partitions
	twin = NULL
	if ( length(part_idx$trans) == 0 ) {
		do_parts = "cis"
		twin = "cis"
	} else if ( length(part_idx$cis) == 0 ) {
		do_parts = "trans"
		twin = "trans"
	}
	gw_is_cis = gw_is_cis && length(part_idx$trans) == 0

	for (part in do_parts) {
		idx = part_idx[[part]]
		if ( length(idx) < 1 ) next   # nothing to predict from
		geno_part = train_geno[ , idx , drop=FALSE ]

		#gemma reads from disk, so blup needs a bfile restricted to this partition
		part_bfile = cv.file
		if ( "blup" %in% models && part != "gw" && length(idx) > 1 ) {
			snplist = paste0(cv.file,".",part,".snps")
			write.table( bim[idx,2] , quote=F , row.names=F , col.names=F , file=snplist )
			part_bfile = paste0(cv.file,".",part)
			run_plink(c("--allow-no-sex", "--bfile", shQuote(cv.file), "--extract",
			            shQuote(snplist), "--make-bed", "--out", shQuote(part_bfile)),
			          part_bfile, sprintf("%s partition extraction for %s fold_%s",
			                             part, opt$gene, fold))
		}

		for (mod in 1:M) {
			print(paste0(part," ",models[mod]))
			#a fit that fails on one small partition should not cost us the whole gene
			pred.wgt = tryCatch({
				if ( length(idx) == 1 ) {
					#cv.glmnet rejects a single-column matrix, but CV on one SNP is well
					#defined: every model collapses to the same simple regression. genos
					#are scaled, so the marginal coefficient is the OLS slope.
					weights.marginal( geno_part , as.matrix(train_pheno[,3]) , beta=T )
				}
				else if (models[mod] == 'xtune') {
					results = weights.xtune( geno_part , as.matrix(train_pheno[,3]) , alpha=1,fold,idx)
					alphas = results$alphas
					results$weights
				}
				else if (models[mod] == 'lasso') {
					weights.enet( geno_part , as.matrix(train_pheno[,3]) , alpha=1 )
				}
				else if (models[mod] == 'enet') {
					weights.enet( geno_part , as.matrix(train_pheno[,3]) , alpha=0.5 )
				}
				else if ( models[mod] == "blup" ) {
					weights.bslmm( part_bfile , bv_type=2 , snp=bim[idx,2] , fold=fold , tag=part )
				}
				else if ( models[mod] == "top1" ) {
					w = weights.marginal( geno_part , as.matrix(train_pheno[,3]) , beta=T )
					w[ - which.max( w^2 ) ] = 0
					w
				}
			}, error = function(e) {
				message(paste0("Fit failed for ",part," / ",models[mod],": ",conditionMessage(e)))
				NULL
			})
			if ( is.null(pred.wgt) ) next   # leave this partition/model as NA

			pred.wgt[ is.na(pred.wgt) ] = 0
			cv.calls[[part]][ -keep_geno , mod ] = genos$bed[ -keep_geno , idx , drop=FALSE ] %*% as.matrix(pred.wgt)
		}
	}

	#exact copy, not an approximation: part_idx[[twin]] == part_idx$gw here
	if ( !is.null(twin) ) cv.calls[["gw"]][ -keep_geno , ] = cv.calls[[twin]][ -keep_geno , ]

	#geno_part/pred.wgt do not exist if no partition was fittable
	rm(list=intersect(c("genos","train_geno","train_pheno","keep_geno","keep_pheno","bim","cis","cis_snps","trans_snps","pred.wgt","part_idx","geno_part"),ls()))
	gc()

}


#cells start as NA, so partitions that were never fit simply stay NA
for ( part in partitions ) {
	for ( mod in 1:M ) {
		ok = which(!is.na(cv.calls[[part]][,mod]))
		s = if ( length(ok) ) sd(cv.calls[[part]][ok,mod]) else NA
		if ( !is.na(s) && s != 0 ) {
			reg = summary(lm( pheno[ok,3] ~ cv.calls[[part]][ok,mod] ))
			cv.performance[[part]][ 1, mod ] = reg$adj.r.sq
			cv.performance[[part]][ 2, mod ] = reg$coef[2,4]
		}
	}
}

print('gw')
print(cv.performance[["gw"]])
print('cis')
print(cv.performance[["cis"]])
print('trans')
print(cv.performance[["trans"]])

#pick the winning partition from the separate CV R2s, before the final fit
best_r2 = sapply( partitions , function(p) {
	v = suppressWarnings(max( cv.performance[[p]][1,] , na.rm=TRUE ))
	if ( is.finite(v) ) v else -Inf
} )
best_model = partitions[ which.max(best_r2) ]
print(paste0("Best model is ",best_model))

# --- FULL ANALYSES

genos = read_plink(paste0(opt$working_dir,"/fold_0/",opt$gene),impute="avg")
genos$bed = scale(genos$bed)
print(dim(genos$bed))
print(dim(pheno[,3]))

pheno_file = paste0(opt$working_dir,"/fold_0/",opt$gene)
geno.file = paste(opt$working_dir,"/fold_0/",opt$gene,".cv",sep='')
write.table( pheno , quote=F , row.names=F , col.names=F , file=paste(geno.file,".keep",sep=''))
run_plink(c("--allow-no-sex", "--bfile", shQuote(pheno_file), "--keep",
            shQuote(paste0(geno.file, ".keep")), "--out", shQuote(geno.file), "--make-bed"),
          geno.file, sprintf("final extraction for %s", opt$gene))

bim = genos$bim
cis = (bim$V1 == chr & bim$V4 > u_bound & bim$V4 < l_bound)

#fold_0's SNP set can differ from the CV folds', so recompute the partitions here
part_idx = list( gw = seq_len(ncol(genos$bed)) , cis = which(cis) , trans = which(!cis) )

#with no trans SNPs anywhere, gw and cis are the same SNP set and fit, so report the
#more specific label. Both conditions required, or a gw model gets reported as cis.
if ( best_model == "gw" && gw_is_cis && length(part_idx$trans) == 0 ) {
	best_model = "cis"
	print("No trans SNPs in any fold, reporting best model as cis")
}

idx = part_idx[[best_model]]
if ( length(idx) < 1 ) {
	print(paste0("Partition ",best_model," is empty at fold_0, falling back to gw"))
	best_model = "gw"
	idx = part_idx[["gw"]]
}

snps = genos$bim[idx,]
wgt.matrix = matrix(0,nrow=length(idx),ncol=M)
colnames(wgt.matrix) = models
rownames(wgt.matrix) = genos$bim[idx,2]

alphas = matrix()

part_bfile = geno.file
if ( "blup" %in% models && best_model != "gw" && length(idx) > 1 ) {
        snplist = paste0(geno.file,".",best_model,".snps")
        write.table( bim[idx,2] , quote=F , row.names=F , col.names=F , file=snplist )
        part_bfile = paste0(geno.file,".",best_model)
        run_plink(c("--allow-no-sex", "--bfile", shQuote(geno.file), "--extract",
                    shQuote(snplist), "--make-bed", "--out", shQuote(part_bfile)),
                  part_bfile, sprintf("%s final partition extraction for %s",
                                     best_model, opt$gene))
}

for (mod in 1:M) {
        if ( length(idx) == 1 ) {
                #same single-SNP fallback as the CV loop above: glmnet cannot take one
                #column, and every model collapses to the same simple regression here
                wgt.matrix[,mod] = weights.marginal( genos$bed[,idx,drop=FALSE] , as.matrix(pheno[,3]) , beta=T )
        }
        else if (models[mod] == 'xtune') {
                results = weights.xtune( genos$bed[,idx,drop=FALSE] , as.matrix(pheno[,3]) , alpha=1,0,idx)
                wgt.matrix[,mod] = results$weights
                alphas = results$alphas
        }
        else if (models[mod] == 'lasso') {
                wgt.matrix[,mod] = weights.enet( genos$bed[,idx,drop=FALSE] , as.matrix(pheno[,3]) , alpha=1 )
        }
        else if (models[mod] == 'enet') {
                wgt.matrix[,mod] = weights.enet( genos$bed[,idx,drop=FALSE] , as.matrix(pheno[,3]) , alpha=0.5 )
        }
        else if ( models[mod] == "blup" ) {
                wgt.matrix[,mod] = weights.bslmm( part_bfile , bv_type=2 , snp=snps[,2] , fold=0 , tag=best_model )
        }
        else if ( models[mod] == "top1" ) {
                wgt.matrix[,mod] = weights.marginal( genos$bed[,idx,drop=FALSE] , as.matrix(pheno[,3]) , beta=F )
        }

}

#7.1_summarize_EGRET_models.R reads all three by name; cv.performance is the winning partition's, column-aligned with wgt.matrix
cv.performance_gw    = cv.performance[["gw"]]
cv.performance_cis   = cv.performance[["cis"]]
cv.performance_trans = cv.performance[["trans"]]
cv.performance       = cv.performance[[best_model]]


print(paste0("The weight dir is ",opt$weights_dir ,"/",opt$gene, ".wgt.RDat" ))
save( wgt.matrix , snps , cv.performance,cv.performance_gw,cv.performance_cis, cv.performance_trans,best_model,alphas, file = paste0( opt$weights_dir ,"/",opt$gene, ".wgt.RDat" ) )

cleanup()
