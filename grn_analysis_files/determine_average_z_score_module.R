library(data.table)
library(enrichR)

# List of 10 tissues
tissues = c(
  "Adipose_Visceral_Omentum", "Artery_Tibial", "Brain_Cortex",
  "Colon_Transverse", "Esophagus_Mucosa", "Heart_Left_Ventricle",
  "Lung", "Muscle_Skeletal", "Skin_Not_Sun_Exposed_Suprapubic",
  "Whole_Blood"
)

tissues = c(
  "Adipose_Subcutaneous","Adipose_Visceral_Omentum","Artery_Tibial","Brain_Cortex","Colon_Transverse",
  "Esophagus_Mucosa","Heart_Left_Ventricle","Lung","Muscle_Skeletal","Skin_Not_Sun_Exposed_Suprapubic",
  "Whole_Blood","Adrenal_Gland","Artery_Aorta","Artery_Coronary","Brain_Amygdala",
  "Brain_Anterior_cingulate_cortex_BA24","Brain_Caudate_basal_ganglia","Brain_Cerebellar_Hemisphere",
  "Brain_Cerebellum","Brain_Frontal_Cortex_BA9","Brain_Hippocampus","Brain_Hypothalamus",
  "Brain_Nucleus_accumbens_basal_ganglia","Brain_Putamen_basal_ganglia",
  "Brain_Spinal_cord_cervical_c-1","Brain_Substantia_nigra","Breast_Mammary_Tissue",
  "Cells_Cultured_fibroblasts","Cells_EBV-transformed_lymphocytes","Colon_Sigmoid",
  "Esophagus_Gastroesophageal_Junction","Esophagus_Muscularis","Heart_Atrial_Appendage",
  "Kidney_Cortex","Liver","Minor_Salivary_Gland","Nerve_Tibial","Ovary","Pancreas","Pituitary",
  "Prostate","Skin_Sun_Exposed_Lower_leg","Small_Intestine_Terminal_Ileum","Spleen",
  "Stomach","Testis","Thyroid","Uterus","Vagina"
)

general_sumstats = c(
  "PASS_IBD_deLange2017.sumstats.gz", "PASS_Rheumatoid_Arthritis.sumstats.gz", 
  "UKB_460K.blood_PLATELET_COUNT.sumstats.gz", "UKB_460K.blood_RBC_DISTRIB_WIDTH.sumstats.gz", 
  "UKB_460K.blood_RED_COUNT.sumstats.gz", "UKB_460K.blood_WHITE_COUNT.sumstats.gz", "UKB_460K.disease_ALLERGY_ECZEMA_DIAGNOSED.sumstats.gz", 
  "UKB_460K.disease_HYPOTHYROIDISM_SELF_REP.sumstats.gz", "PASS_AtrialFibrillation_Nielsen2018.sumstats.gz", 
  "PASS_HDL.sumstats.gz", "PASS_IschemicStroke_Malik2018.sumstats.gz", "PASS_LDL.sumstats.gz", 
  "UKB_460K.biochemistry_Cholesterol.sumstats.gz", "UKB_460K.body_BMIz.sumstats.gz", 
  "UKB_460K.bp_DIASTOLICadjMEDz.sumstats.gz"
)

general_sumstats = c('UKB_460K.body_BMIz.sumstats.gz', 'UKB_460K.cov_EDU_YEARS.sumstats.gz', 'UKB_460K.lung_FVCzSMOKE.sumstats.gz', 'PASS_Intelligence_SavageJansen2018.sumstats.gz', 'UKB_460K.mental_NEUROTICISM.sumstats.gz', 'UKB_460K.bp_DIASTOLICadjMEDz.sumstats.gz', 'UKB_460K.blood_WHITE_COUNT.sumstats.gz', 'UKB_460K.body_HEIGHTz.sumstats.gz', 'UKB_460K.other_MORNINGPERSON.sumstats.gz', 'UKB_460K.repro_MENARCHE_AGE.sumstats.gz', 'UKB_460K.biochemistry_TotalProtein.sumstats.gz', 'UKB_460K.body_WHRadjBMIz.sumstats.gz', 'UKB_460K.lung_FEV1FVCzSMOKE.sumstats.gz', 'PASS_ReactionTime_Davies2018.sumstats.gz', 'PASS_Myopia_Hysi2020.sumstats.gz', 'UKB_460K.biochemistry_Creatinine.sumstats.gz', 'UKB_460K.blood_RED_COUNT.sumstats.gz', 'PASS_Schizophrenia_Pardinas2018.sumstats.gz', 'UKB_460K.blood_PLATELET_COUNT.sumstats.gz', 'UKB_460K.bmd_HEEL_TSCOREz.sumstats.gz', 'PASS_SleepDuration_Dashti2019.sumstats.gz', 'UKB_460K.biochemistry_IGF1.sumstats.gz', 'PASS_GeneralRiskTolerance_KarlssonLinner2019.sumstats.gz', 'PASS_Eosino_Vuckovic2020.sumstats.gz', 'UKB_460K.biochemistry_AspartateAminotransferase.sumstats.gz', 'UKB_460K.blood_RBC_DISTRIB_WIDTH.sumstats.gz', 'PASS_Insomnia_Jansen2019.sumstats.gz', 'PASS_Height1.sumstats.gz', 'PASS_LymP_Vuckovic2020.sumstats.gz', 'UKB_460K.repro_NumberChildrenEverBorn_Pooled.sumstats.gz', 'PASS_AtrialFibrillation_Nielsen2018.sumstats.gz', 'UKB_460K.biochemistry_Testosterone_Male.sumstats.gz', 'PASS_BMI1.sumstats.gz', 'UKB_460K.disease_HYPOTHYROIDISM_SELF_REP.sumstats.gz', 'PASS_AgeFirstBirth.sumstats.gz', 'PASS_Glaucoma_Craig2020.sumstats.gz', 'PASS_RTC_Vuckovic2020.sumstats.gz', 'UKB_460K.body_BALDING1.sumstats.gz', 'UKB_460K.biochemistry_AlkalinePhosphatase.sumstats.gz', 'UKB_460K.disease_ALLERGY_ECZEMA_DIAGNOSED.sumstats.gz', 'PASS_DrinksPerWeek_Liu2019.sumstats.gz', 'UKB_460K.biochemistry_Phosphate.sumstats.gz', 'PASS_MedicationUse_Wu2019.sumstats.gz', 'PASS_IBD_deLange2017.sumstats.gz', 'UKB_460K.biochemistry_VitaminD.sumstats.gz', 'PASS_MDD_Wray2018.sumstats.gz', 'PASS_RD_Zhao2021.sumstats.gz', 'UKB_460K.biochemistry_Cholesterol.sumstats.gz', 'PASS_ADHD_Demontis2018.sumstats.gz', 'PASS_MonoP_Vuckovic2020.sumstats.gz', 'UKB_460K.biochemistry_TotalBilirubin.sumstats.gz', 'UKB_460K.repro_MENOPAUSE_AGE.sumstats.gz', 'PASS_SA_Grasby2020.sumstats.gz', 'PASS_HipOA_Tachmazidou2019.sumstats.gz', 'PASS_AnorexiaNervosa_Watson2019.sumstats.gz', 'UKB_460K.pigment_SUNBURN.sumstats.gz', 'PASS_NumberChildrenEverBorn.sumstats.gz', 'PASS_BipolarDisorder_Ruderfer2018.sumstats.gz', 'PASS_Years_of_Education1.sumstats.gz', 'PASS_PancreasVol_Liu2021.sumstats.gz', 'PASS_HDL.sumstats.gz', 'PASS_CaudateVol_Satizabal2019.sumstats.gz', 'PASS_BrainstemVol_Satizabal2019.sumstats.gz', 'PASS_TH_Grasby2020.sumstats.gz', 'PASS_LDL.sumstats.gz', 'PASS_BasoP_Vuckovic2020.sumstats.gz', 'PASS_IschemicStroke_Malik2018.sumstats.gz', 'PASS_Anorexia.sumstats.gz', 'PASS_Alzheimers_deRojas2021.sumstats.gz', 'PASS_Rheumatoid_Arthritis.sumstats.gz', 'PASS_Ever_Smoked.sumstats.gz', 'PASS_MO_Zhao2021.sumstats.gz', 'PASS_ProstateCancer.sumstats.gz', 'UKB_460K.cancer_PROSTATE.sumstats.gz', 'PASS_AccumbensVol_Satizabal2019.sumstats.gz', 'UKB_460K.cancer_BREAST.sumstats.gz', 'PASS_Type_2_Diabetes.sumstats.gz', 'PASS_Autism.sumstats.gz')


#general_sumstats = c(
#  "PASS_IBD_deLange2017.sumstats.gz", "PASS_Rheumatoid_Arthritis.sumstats.gz", 
#  "UKB_460K.blood_PLATELET_COUNT.sumstats.gz", "UKB_460K.blood_RBC_DISTRIB_WIDTH.sumstats.gz", 
#  "UKB_460K.blood_RED_COUNT.sumstats.gz", "UKB_460K.blood_WHITE_COUNT.sumstats.gz", "PASS_LDL.sumstats.gz", 
#  "UKB_460K.biochemistry_Cholesterol.sumstats.gz", "UKB_460K.body_BMIz.sumstats.gz"
#)

gene_info = fread("../../data/GTEx_V8.txt.gz", header = TRUE)
significant_results = data.frame(
  tissue = character(),
  trait = character(),
  module = character(),
  pathway_name = character(),
  p_value = numeric(),
  overlap_genes = numeric(),
  avg_z_EGRET = numeric(),
  avg_z_cis = numeric(),
  num_trans_regulated = numeric(),
  prop_trans_regulated = numeric()
)

for (tissue in tissues) {
  model_performance = fread(paste0("../results_sumstats/",tissue,"/cis_MatrixeQTL_GBAT_transPCO_FDR_0.1.txt"), header = T)

  for (trait in general_sumstats) {
    trait_file = paste0(substr(trait, 1, nchar(trait) - 12), ".dat")
	  print(trait_file)
    TWAS_path = paste0("../TWAS_association_/results/", tissue, "/cis_MatrixeQTL_GBAT_transPCO_FDR_0.1/", trait_file)
    cis_path = paste0("../TWAS_association_results/", tissue, "/cis/", trait_file)

    if (!file.exists(TWAS_path) || !file.exists(cis_path)) {
      print("Skipping", trait, "TWAS results not present")
      next
    }
    TWAS_results = fread(TWAS_path, header = TRUE)
    cis_TWAS_results = fread(cis_path, header = TRUE)

    TWAS_results[is.na(TWAS_results$TWAS.Z), "TWAS.Z"] = 0
    cis_TWAS_results[is.na(cis_TWAS_results$TWAS.Z), "TWAS.Z"] = 0

    modules = list.files(paste0("coregulation_grns/",tissue,"/"))

    for (module in modules) {
      genes_in_module = fread(paste0("coregulation_grns/",tissue, "/", module), header = FALSE)
      gene_ids = gene_info$name[gene_info$geneId %in% genes_in_module$V1]

      TWAS_module_subset = merge(TWAS_results, genes_in_module, by.x = 'ID', by.y = 'V1', all.y = TRUE)
      TWAS_module_subset[is.na(TWAS_module_subset$TWAS.Z), "TWAS.Z"] = 0

      cis_module_subset = merge(cis_TWAS_results, genes_in_module, by.x = 'ID', by.y = 'V1', all.y = TRUE)
      cis_module_subset[is.na(cis_module_subset$TWAS.Z), "TWAS.Z"] = 0

      avg_z_EGRET = mean(abs(TWAS_module_subset$TWAS.Z), na.rm = TRUE)
      avg_z_cis = mean(abs(cis_module_subset$TWAS.Z), na.rm = TRUE)

      z_egret <- abs(TWAS_module_subset$TWAS.Z)
      z_cis   <- abs(cis_module_subset$TWAS.Z)

      # Ensure pairing is valid
      stopifnot(length(z_egret) == length(z_cis))

      d = z_egret - z_cis

      if (length(d) < 2 || sd(d) == 0) {
          t_test_results = list(
            statistic = NA,
            p.value = 1,
            estimate = mean(d),
            method = "paired t-test (skipped: zero variance)"
          )
        } else {
          t_test_results = t.test(
            z_egret,
            z_cis,
            paired = TRUE,
            alternative = "greater"
          )
        }

      if (length(unique(d)) > 1) {
        wilcoxon_test = wilcox.test(
          z_egret,
          z_cis,
          paired = TRUE,
          alternative = "greater",
          exact = FALSE
        )
        wilcoxon_p = wilcoxon_test$p.value
      } else {
        wilcoxon_p = 1
      }

      module_gene_performance = model_performance[gene %in% unlist(genes_in_module),]
      module_gene_performance$p_trans_part[is.na(module_gene_performance$p_trans_part)] = 1

      num_trans_regulated = sum(module_gene_performance$p_trans_part < 0.01)
      prop_trans_regulated = num_trans_regulated/nrow(genes_in_module)


      #dbs <- c("GO_Biological_Process_2023")
      #result <- enrichr(gene_ids, dbs)
      result = data.frame()
      if (!is.null(result[["GO_Biological_Process_2023"]]) &&
          nrow(result[["GO_Biological_Process_2023"]]) > 0) {

        top_pathway <- result[["GO_Biological_Process_2023"]][
          which.min(result[["GO_Biological_Process_2023"]]$Adjusted.P.value),
        ]

        adj_p <- top_pathway$Adjusted.P.value
        overlap = as.numeric(strsplit(top_pathway$Overlap, "/")[[1]][1])
      } else {
        # Fallback when no enrichment results
        adj_p <- 1
        top_pathway <- data.frame(
          Term = NA,
          Adjusted.P.value = 1,
          overlap = NA
        )
        adj_p <- top_pathway$Adjusted.P.value
        overlap = top_pathway$overlap
      }
      significant_results = rbind(significant_results, 
        data.frame(
          tissue = tissue,
		      trait = trait_file,
          module = module,
          pathway_name = top_pathway$Term,
          p_value = adj_p,
          overlap_genes = overlap,
          avg_z_EGRET = avg_z_EGRET,
          avg_z_cis = avg_z_cis,
          num_trans_regulated = num_trans_regulated,
          prop_trans_regulated = prop_trans_regulated,
          module_size = length(gene_ids),
          t_test_p = t_test_results$p.value,
          wilcoxon_p = wilcoxon_p
        )
      )

    }
  }
}

fwrite(significant_results, "coregulation_significant_enrichment_results_new.txt", sep = "\t", quote = FALSE, row.names = FALSE)
