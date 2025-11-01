#################### Please refer to the Methods section for detailed methodology ####################
#################### This outlines only the general data analysis workflow; for specific data processing details, please contact the author if needed ####################



setwd(" ")
library(Seurat)
library(multtest)
library(dplyr)
library(patchwork)
library(R.utils)
library(ggplot2)
library(ggrepel)
library(TwoSampleMR)
library(plinkbinr)
library(ieugwasr)
library(clusterProfiler)
library(tidyverse)
library(forcats)
library(org.Hs.eg.db)
library(pathview)
library(enrichplot)
library(ggnewscale)
library(DOSE)
library(stringr)
library(forestplot)
library(data.table)



raw_counts <- readRDS("snRNA-seq_human_raw_counts.RDS")
dim(raw_counts)
barcode_meta <- read.csv("snRNA-seq_human_barcode_meta.csv", row.names = 1)
seuratobject <- CreateSeuratObject(counts = raw_counts,
                           meta.data = barcode_meta,
                           project = "humanTG",
                           min.features = 1000,
                           min.cells = 3)
dim(seuratobject)
saveRDS(seuratobject, "snRNA-seq_human_raw_counts_fea1000cell3.RDS")


seuratobject[["percent.mt"]] <- PercentageFeatureSet(seuratobject, pattern = "^MT-") 
head(seuratobject@meta.data,5)

seuratobject <- subset(seuratobject, subset = nFeature_RNA < 15000 & percent.mt < 5)

seuratobject <- NormalizeData(seuratobject, normalization.method = "LogNormalize", scale.factor = 10000)

seuratobject <- FindVariableFeatures(seuratobject, selection.method = "vst", nfeatures = 2000)

hvf_info <- HVFInfo(seuratobject)

hvf_info_filtered <- hvf_info[hvf_info$mean > 0 & hvf_info$variance > 0 & hvf_info$variance.standardized > 0, ]


seuratobject <- ScaleData(object = seuratobject, vars.to.regress = c("nCount_RNA", "percent.mt"))

seuratobject <- RunPCA(seuratobject, features = VariableFeatures(object = seuratobject), verbose=F)

print(seuratobject[["pca"]], dims = 1:5, nfeatures = 5)


seuratobject <- RunUMAP(seuratobject, dims = 1:20)
DimPlot(seuratobject, reduction = "umap", group.by = "subtype", label = TRUE) +
  ggtitle("UMAP of Cell Subtypes")


seuratobject <- RunTSNE(seuratobject, dims = 1:20)
DimPlot(seuratobject, reduction = "tsne", group.by = "subtype", label = TRUE) +
  ggtitle("t-SNE of Cell Subtypes")


subtype <- SetIdent(seuratobject, value = seuratobject@meta.data$subtype)
marker <- FindAllMarkers(subtype, only.pos = F, min.pct = 0.25, logfc.threshold = 0.5)



filtered_allmarker <- marker %>% filter(p_val_adj < 0.05)


brain <- data.table::fread("PsychENCODE", header = TRUE, sep = ",")

eqtl <- merge(filtered_allmarker, brain, by.x = "ENSEMBL", by.y = "gene_id", all.x = T, all.y = F)


disease <- read_outcome_data(filename = " ", 
                         sep = " ", 
                         snp_col = " ", 
                         beta_col = " ", 
                         se_col = " ", 
                         effect_allele_col = " ", 
                         other_allele_col = " ", 
                         pval_col = " ", 
                         eaf_col = " ")


grouped_df <- group_by(eqtl, ENSEMBL)
grouped_list <- group_split(grouped_df)

harmo_allSNP <- data.frame()
mr_results <- data.frame()


for (i in 1 : length(grouped_list)) {
  print("current i:")
  print(i)
  current_gene <- grouped_list[[i]]
  print("SNPsum")
  print(nrow(current_gene))
  write.csv(current_gene, file = "tmp.csv")
  
  exposure_data <- read_exposure_data("tmp.csv", sep = " ", snp_col = " ", beta_col = " ", se_col = " ", effect_allele_col = " ", other_allele_col = " ", pval_col = " ")
  

  p_error <- tryCatch({
    clumpSNP <- ld_clump(
      dplyr::tibble(rsid = exposure_data$SNP, pval = exposure_data$pval.exposure, id=exposure_data$id.exposure),
      plink_bin = " ",
      bfile = " "
    )
  }, error = function(e) e)
  if(inherits(p_error, "error")) next
  merge_clump <- merge(exposure_data, clumpSNP, by.x = "SNP", by.y = "rsid", all.x = FALSE)
  new_clump <- subset(merge_clump, select = -c(pval, id))
  
  
  OUTCOME <- filter(disease, SNP %in% new_clump$SNP)
  
  harmo_data <- harmonise_data(exposure_dat = exposure_data, outcome_dat = OUTCOME)
  harmo_data <- harmo_data[harmo_data$mr_keep, ]
  if (nrow(harmo_data) == 0) next
  
  harmo_data$ENSEMBL_id <- current_gene$ENSEMBL[1]
  harmo_allSNP <- rbind(harmo_allSNP, harmo_data)
  
  if(nrow(harmo_data) > 1) {
    mr_result <- generate_odds_ratios(mr_res = mr(harmo_data, method_list = "mr_ivw"))
    mr_result_Q <- mr_heterogeneity(harmo_data) 
    mr_result$Q_pval <- mr_result_Q$Q_pval[1]
    mr_result$ENSEMBL_id = harmo_data$ENSEMBL_id[1]
    
    if(nrow(harmo_data) > 2) {
      mr_result_EGGER <- mr_pleiotropy_test(harmo_data)
      mr_result$MREgger_intercept <- mr_result_EGGER$pval
      mr_results <- rbind(mr_results, mr_result)
    } else {
      mr_result$MREgger_intercept <- NA
      mr_results <- rbind(mr_results, mr_result)
    }
  } else {
    mr_result <- generate_odds_ratios(mr_res = mr(harmo_data, method_list = "mr_wald_ratio"))
    mr_result$ENSEMBL_id = harmo_data$ENSEMBL_id[1]
    mr_result$Q_pval <- NA
    mr_result$MREgger_intercept <- NA
    mr_results <- rbind(mr_results, mr_result)
  }
}


