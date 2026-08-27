library('data.table')
library('optparse')

option_list = list(
  make_option("--tissue", action="store", default=NA, type='character',
              help="tissue"),
  make_option("--base_dir", action="store", default=NA, type='character',
              help="base directory prefixed to all file paths"),
  make_option("--MatrixeQTL_bed_dir", action="store", default=NA, type='character',
              help="the dir containing MatrixeQTL bed files"),
  make_option("--GBAT_bed_dir", action="store", default=NA, type='character',
              help="the dir containing GBAT bed files"),
  make_option("--transPCO_bed_dir", action="store", default=NA, type='character',
              help="the dir which contains transPCO bed files"),
  make_option("--exclude_crossmap", action="store", default=TRUE, type='logical',
              help="if True, will exclude regions which are cross mappable to that gene"),
  make_option("--crossmap_dir", action="store", default="0_mismatches", type='character',
              help="subdirectory under bed_files/cross_mapped/ for cross-mappable regions"),
  make_option("--output_dir", action = 'store', default=NA, type = 'character',
	          help = 'directory to store bed files, plink results, and z matrices'),
  make_option("--models", action = 'store', default=NA, type = 'character',
	          help = 'Cis, MatrixeQTL, GBAT, transPCO'),
  make_option("--fold", action = 'store', default=0, type = 'numeric',
              help = 'should be 0:num_folds'),
  make_option("--genotype_prefix", action="store", default="GTEX_v8_genotypes_pruned",type='character',
              help="prefix for genotype files (.bim/.bed/.fam)"),
  make_option("--plink_path", action="store", default="plink2", type='character',
              help="path to plink2 executable"),
  make_option("--expression_file", action="store", default=NA, type='character',
              help="optional expression file to define gene list; if not provided, uses cis bed dir"),
  make_option("--gene_information", action="store", default=NA, type='character',
              help="Path to file containing gene information"),
  make_option("--cis_window_size", action="store", default=1000000, type='numeric',
              help="size of cis window to create")
  )

opt = parse_args(OptionParser(option_list=option_list))
tissue = opt$tissue
fold = opt$fold
base_dir = opt$base_dir

if (!is.na(opt$expression_file)) {
  expression = fread(opt$expression_file, header = T)
  genes = expression[[1]]
} else {
  print("please use valid expression file")
  q()
}
all_gene_info = fread(opt$gene_information, header = T)

allsnps = fread(paste0(base_dir, "/genotype_files/", opt$genotype_prefix, ".bim"), header = F)

output_dir = opt$output_dir

models = unique( c(unlist(strsplit(opt$models,','))) )

crossmap_dir = opt$crossmap_dir
for (gene in genes) {
	bed_file = matrix(ncol = 6,nrow = 0)

	gene_name = gene

	gene_info = all_gene_info[all_gene_info$geneId == gene_name, ]
	chr = strsplit(gene_info$'#chrom', split = 'chr')[[1]][2]
	tss = gene_info$chromStart
	if (chr == 'X' | chr == 'Y') {
		next
	}

	#cis
	cis_bed_path = paste0(base_dir, "/bed_files/",tissue,"/cis/",gene)
	if (file.exists(cis_bed_path)) {
		cis_bed = fread(cis_bed_path,header = F)
	} else {
		bounds_u = max(1, tss - floor(opt$cis_window_size/2))
		bounds_l = tss + floor(opt$cis_window_size/2)
		bed = c(chr, bounds_u, bounds_l, gene_name, 0, "+")
		cis_bed = matrix(bed, nrow = 1, ncol = 6)
	}

	#GBAT
	GBAT_bed_path = paste0(base_dir, "/GBAT/",tissue,"/fold_",fold,"/",opt$GBAT_bed_dir,"/",gene_name,".txt")
	if (file.exists(GBAT_bed_path)) {
		GBAT_snps = fread(GBAT_bed_path,header = F)
	} else {
		GBAT_snps = data.table()
	}
	# Matrix eQTL
	MatrixeQTL_bed_path = paste0(base_dir, "/MatrixeQTL/",tissue,"/fold_",fold,"/",opt$MatrixeQTL_bed_dir,"/",gene_name,".txt")
	if (file.exists(MatrixeQTL_bed_path)) {
		MatrixeQTL_snps = fread(MatrixeQTL_bed_path,header = F)
	} else {
		MatrixeQTL_snps = data.table()
	}

	# trans-PCO
	transPCO_bed_path = paste0(base_dir, "/transPCO/",tissue,"/fold_",fold,"/",opt$transPCO_bed_dir,"/",gene_name,".txt")
	if (file.exists(paste0(transPCO_bed_path)) ) {
		transPCO_snps = fread(transPCO_bed_path,header = F)
	} else {
		transPCO_snps = data.table()
	}

	GBAT_matches = allsnps[match(unlist(GBAT_snps$V1),allsnps$V2),]
	if (nrow(GBAT_matches) == 0) {
		GBAT_bed = c()
	} else {
		GBAT_bed = cbind(GBAT_matches$V1,GBAT_matches$V4,GBAT_matches$V4+1, GBAT_matches$V2,"0","+")
	}
	MatrixeQTL_matches = allsnps[match(unlist(MatrixeQTL_snps$V1),allsnps$V2),]
	if (nrow(MatrixeQTL_matches) == 0) {
		MatrixeQTL_bed = c()
	} else {
		MatrixeQTL_bed = cbind(MatrixeQTL_matches$V1,MatrixeQTL_matches$V4,MatrixeQTL_matches$V4+1, MatrixeQTL_matches$V2,"0","+")
	}
	transPCO_matches = allsnps[match(unlist(transPCO_snps$V1),allsnps$V2),]
	if(nrow(transPCO_matches) == 0) {
		transPCO_bed = c()
	} else {
		transPCO_matches = transPCO_matches[!is.na(transPCO_matches$V2),]
		transPCO_bed = cbind(transPCO_matches$V1,transPCO_matches$V4,transPCO_matches$V4+1,transPCO_matches$V2, "0","+")
	}

	if ("Cis" %in% models) { bed_file = rbind(bed_file,cis_bed) }
	if ("MatrixeQTL" %in% models ) { bed_file = rbind(bed_file,MatrixeQTL_bed) }
	if ("GBAT" %in% models) { bed_file = rbind(bed_file,GBAT_bed) }
	if ("transPCO" %in% models) { bed_file = rbind(bed_file,transPCO_bed) }

	bed_file_dir = paste0(base_dir, "/bed_files/",tissue,"/",output_dir,"/fold_",fold,"/")
	dir.create(bed_file_dir,recursive = T)
	fwrite(bed_file,paste0(bed_file_dir,gene),row.names = F , col.names = F,quote = F, sep = '\t')

	plink_output_dir = paste0(base_dir, "/plink_results/",tissue,"/",output_dir,"/fold_",fold,"/")
	dir.create(plink_output_dir,recursive = T)

	crossmap_bed_path = paste0(base_dir, "/cross_mapped/",crossmap_dir,"/",gene,".bed")
	if(file.exists(crossmap_bed_path) & opt$exclude_crossmap) {
		print("exluding cross mappable genes")
		arg = paste0(opt$plink_path, " --bfile ", base_dir, "/genotype_files/", opt$genotype_prefix, " --extract bed1 ",bed_file_dir,gene, " --make-bed --exclude bed1 ", crossmap_bed_path, " --out ",plink_output_dir,gene_name)
	} else {
		print("no cross mappable genes to exlude")
		arg = paste0(opt$plink_path, " --bfile ", base_dir, "/genotype_files/", opt$genotype_prefix, " --extract bed1 ",bed_file_dir,gene, " --make-bed --out ",plink_output_dir,gene_name)
	}
	system(arg)

	if(file.exists(paste0(plink_output_dir,gene_name,".bim"))) {
		plink_results = fread(paste0(plink_output_dir,gene_name,".bim"), header =  F)
	} else { next }

	MatrixeQTL_matches = match(unlist(MatrixeQTL_snps$V1),plink_results$V2)
	plink_results$MatrixeQTL = 0
	plink_results$MatrixeQTL[MatrixeQTL_matches] = 1

	GBAT_matches = match(unlist(GBAT_snps$V1),plink_results$V2)
	plink_results$GBAT = 0
	plink_results$GBAT[GBAT_matches] = 1

	transPCO_matches = match(unlist(transPCO_snps$V1),plink_results$V2)
	plink_results$transPCO = 0
	plink_results$transPCO[transPCO_matches] = 1

	z_matrix = data.table()
	if ("MatrixeQTL" %in% models ) { z_matrix = cbind(z_matrix,plink_results$MatrixeQTL) }
	if ("GBAT" %in% models) { z_matrix = cbind(z_matrix,plink_results$GBAT) }
	if ("transPCO" %in% models) { z_matrix = cbind(z_matrix,plink_results$transPCO) }

	z_matrix_dir = paste0(base_dir, "/z_matrices/",tissue,"/",output_dir,"/fold_",fold,"/")
	dir.create(z_matrix_dir, recursive = T)
	fwrite(z_matrix,paste0(z_matrix_dir,gene_name,"_z_matrix.txt"),sep = '\t',col.names = F,row.names = F,quote=F)
}
