library('data.table')
library('optparse')
library('glmnet')
library('xtune')
library('plink2R')

option_list = list(
  make_option("--gene", action="store", default=NA, type='character',
	help="gene name"),
  make_option("--working_dir", action="store", default=NA, type='character',
        help="location where bed/bim/fam files are stored"),
  make_option("--models", action="store", default=NA, type='character',
        help="Comma-separated list of prediction models"),
  make_option("--output_dir", action="store", default=NA, type='character',
        help="Path to output files"),
  make_option("--weights_dir",action ='store',default=NA,type='character',
  	help="path to weights dir"),
  make_option("--tissue",action = 'store',default=NA,type='character',
        help="the tissue for analysis"),
  make_option("--PATH_plink",action = 'store',default=NA,type='character',
        help="path to plink"),
  make_option("--PATH_gemma",action = 'store',default=NA,type='character',
        help="path to gemma"),
  make_option("--covar",action = 'store',default=NA,type='character',
        help="covariates"),
  make_option("--z_matrix_dir", action = 'store', default = NA, type = 'character',
	help = 'directory containing z_matrix for each fold'),
  make_option("--clean", action="store", default=FALSE,
        help="whether to clean working dir after run"),
  make_option("--gene_info_file", action="store", default="../data/GTEx_V8.txt.gz", type='character',
        help="path to gene info file (e.g. GTEx_V8.txt.gz)"),
  make_option("--folds", action="store", default=5, type='numeric',
        help="number of cross-validation folds"),
  make_option("--base_dir", action="store", default=NA, type='character',
        help="Path to output dir"))

opt = parse_args(OptionParser(option_list=option_list))
tissue = opt$tissue

models = unique( c(unlist(strsplit(opt$models,','))) )
M = length(models)

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
weights.bslmm = function( input , bv_type , snp , fold ) {
        outdir = paste0(opt$working_dir, "/fold_", fold)
        prefix = paste0(opt$gene, ".BSLMM")
        arg = paste( opt$PATH_gemma , " -miss 1 -maf 0 -r2 1 -rpace 1000 -wpace 1000 -bfile " , input , " -bslmm " , bv_type , " -outdir " , outdir , " -o " , prefix , sep='' )
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
weights.xtune = function(genos, pheno, alpha=1, fold) {
    eff.wgt = matrix(0, ncol=1, nrow=ncol(genos))
    # remove monomorphics
    sds = apply(genos, 2, sd)
    keep = sds != 0 & !is.na(sds)
    print(paste0("The current fold is ",fold))
    
    #z_matrix_args = paste0("Rscript make_z_matrix.R --gene ",opt$gene," --tissue ",opt$tissue," --working_dir ",opt$working_dir," --fold ", fold)
    #system(z_matrix_args)
    z_matrix_file = paste0(opt$z_matrix_dir, "/fold_", fold, "/", opt$gene, "_z_matrix.txt")
    if (file.exists(z_matrix_file)) {
    	Z = read.table(paste0(opt$z_matrix_dir, "/fold_", fold, "/", opt$gene, "_z_matrix.txt"), header = FALSE)
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
                arg = paste("rm -f " , opt$working_dir , "*", sep='')
                system(arg)
        }
}

num_folds = opt$folds
fold_individuals = list()
for (f in 0:num_folds) {
    fold_individuals[[paste0("fold_", f)]] = fread(
        paste0(opt$base_dir,"/fold_", f, "_info/", tissue, "/train_individuals.txt"), header = F
    )
}

if(!is.na(opt$covar)) {
	covar = fread(opt$covar,header = T)
} else {
	covar = NA
}

cv.performance = matrix(NA,nrow=2,ncol=M)
rownames(cv.performance) = c("rsq","pval")
colnames(cv.performance) = models

cv.performance_cis = matrix(NA,nrow=2,ncol=M)
rownames(cv.performance_cis) = c("rsq","pval")
colnames(cv.performance_cis) = models

cv.performance_trans = matrix(NA,nrow=2,ncol=M)
rownames(cv.performance_trans) = c("rsq","pval")
colnames(cv.performance_trans) = models


cv.calls = matrix(NA,nrow=nrow(fold_individuals[["fold_0"]]),ncol=M)

#new
cv.calls_cis = matrix(NA,nrow=nrow(fold_individuals[["fold_0"]]),ncol=M)
cv.calls_trans = matrix(NA,nrow=nrow(fold_individuals[["fold_0"]]),ncol=M)
alphas = c()

print('here')
#new stuff to do cis/trans analysis
all_gene_info = fread(opt$gene_info_file, header = T)
print('here2')
gene_info = all_gene_info[grep(opt$gene, all_gene_info$geneId),]
u_bound = max(0,gene_info$chromStart - 2500000)
l_bound = gene_info$chromStart + 2500000
chr = strsplit(gene_info$'#chrom','chr')[[1]][2]

for (fold in 1:num_folds) {
	print("starting cross val")
	# read in genotypes
	genos = read_plink(paste0(opt$working_dir,"/fold_",fold,"/",opt$gene),impute="avg") 

	# important : genotypes are standardized and scaled here:
	genos$bed = scale(genos$bed)
	pheno = genos$fam[,c(1,2,6)] 
	
	#regress out covariates
	if (!is.na(opt$covar)) {
		m = match( unlist(pheno$V2) ,unlist(covar$ID) )
                m.keep = !is.na(m)
                pheno = pheno[m.keep,]
                m = m[m.keep]
                covar = covar[m,]
                reg = summary(lm(as.matrix(pheno[,3]) ~ as.matrix(covar[,3:ncol(covar)]) ))
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
        arg = paste( opt$PATH_plink ," --allow-no-sex --bfile ",pheno_file," --keep ",cv.file,".keep --out ",cv.file," --make-bed",sep='')
        system(arg)

	for (mod in 1:M) {
		print(models[mod])
		if (models[mod] == 'xtune') {
			results = weights.xtune( train_geno , as.matrix(train_pheno[,3]) , alpha=1,fold)
			pred.wgt = results$weights
			alphas = results$alphas
		} 
		else if (models[mod] == 'lasso') {
			pred.wgt = weights.enet( train_geno , as.matrix(train_pheno[,3]) , alpha=1 )
		}
		else if (models[mod] == 'enet') {
			pred.wgt = weights.enet( train_geno , as.matrix(train_pheno[,3]) , alpha=0.5 )
		}
		else if ( models[mod] == "blup" ) {
                        pred.wgt = weights.bslmm( cv.file , bv_type=2 , snp=genos$bim[,2] , fold=fold )
                }
		else if ( models[mod] == "top1" ) {
                        pred.wgt = weights.marginal( train_geno , as.matrix(train_pheno[,3]) , beta=T )
                        pred.wgt[ - which.max( pred.wgt^2 ) ] = 0
                }

		pred.wgt[ is.na(pred.wgt) ] = 0
                cv.calls[ -keep_geno , mod ] = genos$bed[ -keep_geno , ] %*% pred.wgt
		
		pred.wgt = as.matrix(pred.wgt)
		#have to only keep the genotypes and weights which correspond to cis snps
		cv.calls_cis[ -keep_geno , mod ] = genos$bed[-keep_geno ,cis] %*% as.matrix(pred.wgt[cis,])
		cv.calls_trans[ -keep_geno , mod] = genos$bed[-keep_geno ,!cis] %*% as.matrix(pred.wgt[!cis,])

	}

}

wgt.matrix = matrix(0,nrow=nrow(genos$bim),ncol=M)
colnames(wgt.matrix) = models
rownames(wgt.matrix) = genos$bim[,2]


for ( mod in 1:M ) {
	if ( !is.na(sd(cv.calls[which(!is.na(cv.calls[,mod])),mod])) && sd(cv.calls[which(!is.na(cv.calls[,mod])),mod]) != 0 ) {
 		reg = summary(lm( pheno[which(!is.na(cv.calls[,mod])),3] ~ cv.calls[which(!is.na(cv.calls[,mod])),mod] ))
                cv.performance[ 1, mod ] = reg$adj.r.sq
                cv.performance[ 2, mod ] = reg$coef[2,4]
        } else {
                cv.performance[ 1, mod ] = NA
                cv.performance[ 2, mod ] = NA
        }
	wgt.matrix[,mod] = pred.wgt
}

#no longer write the results
#write.table(cv.performance,paste0(opt$output_dir,"/",opt$gene,".txt"),quote=F,sep='\t')
#new
for ( mod in 1:M ) {
        if ( !is.na(sd(cv.calls_cis[which(!is.na(cv.calls_cis[,mod])),mod])) && sd(cv.calls_cis[which(!is.na(cv.calls_cis[,mod])),mod]) != 0 ) {
                reg = summary(lm( pheno[which(!is.na(cv.calls_cis[,mod])),3] ~ cv.calls_cis[which(!is.na(cv.calls_cis[,mod])),mod] ))
                cv.performance_cis[ 1, mod ] = reg$adj.r.sq
                cv.performance_cis[ 2, mod ] = reg$coef[2,4]
        } else {
                cv.performance_cis[ 1, mod ] = NA
                cv.performance_cis[ 2, mod ] = NA
        }
}


for ( mod in 1:M ) {
        if ( !is.na(sd(cv.calls_trans[which(!is.na(cv.calls_trans[,mod])),mod])) && sd(cv.calls_trans[which(!is.na(cv.calls_trans[,mod])),mod]) != 0 ) {
                reg = summary(lm( pheno[which(!is.na(cv.calls_trans[,mod])),3] ~ cv.calls_trans[which(!is.na(cv.calls_trans[,mod])),mod] ))
                cv.performance_trans[ 1, mod ] = reg$adj.r.sq
                cv.performance_trans[ 2, mod ] = reg$coef[2,4]
        } else {
                cv.performance_trans[ 1, mod ] = NA
                cv.performance_trans[ 2, mod ] = NA
        }
}

# --- FULL ANALYSES

genos = read_plink(paste0(opt$working_dir,"/fold_0/",opt$gene),impute="avg")
genos$bed = scale(genos$bed)
print(dim(genos$bed))
print(dim(pheno[,3]))

pheno_file = paste0(opt$working_dir,"/fold_0/",opt$gene)
geno.file = paste(opt$working_dir,"/fold_0/",opt$gene,".cv",sep='')
write.table( pheno , quote=F , row.names=F , col.names=F , file=paste(geno.file,".keep",sep=''))
arg = paste( opt$PATH_plink ," --allow-no-sex --bfile ",pheno_file," --keep ",geno.file,".keep --out ",geno.file," --make-bed",sep='')
system(arg)

bim = genos$bim
cis = (bim$V1 == chr & bim$V4 > u_bound & bim$V4 < l_bound)
cis_snps = bim[cis,]
trans_snps = bim[!cis,]

wgt.matrix = matrix(0,nrow=nrow(genos$bim),ncol=M)
colnames(wgt.matrix) = models
rownames(wgt.matrix) = genos$bim[,2]


alphas = matrix()
names = c("intercept","Cis","MatrixeQTL","GBAT","transPCO")
for (mod in 1:M) {
        if (models[mod] == 'xtune') {
                results = weights.xtune( genos$bed , as.matrix(pheno[,3]) , alpha=1,0)
                wgt.matrix[,mod] = results$weights
                alphas = results$alphas
        }
        else if (models[mod] == 'lasso') {
                wgt.matrix[,mod] = weights.enet( genos$bed , as.matrix(pheno[,3]) , alpha=1 )
        }
        else if (models[mod] == 'enet') {
                wgt.matrix[,mod] = weights.enet( genos$bed , as.matrix(pheno[,3]) , alpha=0.5 )
        }
        else if ( models[mod] == "blup" ) {
                wgt.matrix[,mod] = weights.bslmm( geno.file , bv_type=2 , snp=genos$bim[,2] , fold=0 )
        }
	else if ( models[mod] == "top1" ) {
                wgt.matrix[,mod] = weights.marginal( genos$bed , as.matrix(pheno[,3]) , beta=F )
        }

}

print('gw')
print(cv.performance)
print('cis')
print(cv.performance_cis)
print('trans')
print(cv.performance_trans)

max_values = c(max(cv.performance[1,], na.rm = TRUE),
                max(cv.performance_cis[1,], na.rm = TRUE),
                max(cv.performance_trans[1,], na.rm = TRUE))

# Find the index of the maximum value
best_model = which.max(max_values)

cv.performance_gw = cv.performance
# Assign wgt.matrix based on the best model
if (best_model == 2) {
    	wgt.matrix = wgt.matrix[cis,,drop = FALSE]
	snps = cis_snps
    print("Best model is cis")
        best_model = "cis"
        cv.performance = cv.performance_cis
} else if (best_model == 3) {
        wgt.matrix = wgt.matrix[!cis,,drop = FALSE]
        snps = trans_snps
        print("Best model is trans")
        best_model = "trans"
        cv.performance = cv.performance_trans
} else {
    print("Best model is gw")
    snps=genos$bim
    best_model = "gw"
}


write.table(alphas,paste0(opt$output_dir,"/",opt$gene,"_alphas.txt"),row.names = F, col.names = F,quote=F,sep='\t')
print(paste0("The weight dir is ",opt$weights_dir ,"/",opt$gene, ".wgt.RDat" ))
save( wgt.matrix , snps , cv.performance,cv.performance_gw,cv.performance_cis, cv.performance_trans,best_model,alphas, file = paste0( opt$weights_dir ,"/",opt$gene, ".wgt.RDat" ) )

cleanup()
