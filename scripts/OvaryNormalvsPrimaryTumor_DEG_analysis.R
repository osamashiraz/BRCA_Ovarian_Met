run_DEseq2_OvaryNormalvsPrimaryTumor = function(cts, clinData, excluded_samples){
  library(DESeq2)
  
  mode(cts)="integer"
  icesTAF::mkdir("../output/deseq2")
  
  
  # ILC Met vs IDC Met
  clinData_select = subset(clinData, TissueType %in% c("Normal","Primary") & Study_Code %in% "Ovary" & SampleID %notin% excluded_samples)
  clinData_select$type_site = paste0(clinData_select$TissueType, clinData_select$TissueSite)
  clinData_select = subset(clinData_select, type_site %notin% "NormalBreast")
  clinData_select$TissueSite = factor(clinData_select$TissueSite,levels = c("Ovary","Breast"))
  clinData_select$TissueType = factor(clinData_select$TissueType,levels = c("Normal","Primary"))
  clinData_select$Patient_ID = factor(clinData_select$Patient_ID)
  table(clinData_select$TissueType)
  
  
  
  ### DEG Analysis ----------------- DO NOT RUN AGAIN
  # library(DESeq2)
  # 
  # dds_OvaryVsPri <- DESeqDataSetFromMatrix(countData = cts[, clinData_select$SampleID], colData = clinData_select, design = ~TissueType)
  # dds_OvaryVsPri <- DESeq(dds_OvaryVsPri, parallel = F)
  # 
  # save(dds_OvaryVsPri,file = "../output/deseq2/dds_OvaryVsPri.Rdata")
  # write.dcf(data.frame(design="= ~TissueType"),file = "../output/deseq2/dds_OvaryVsPri Run Info.txt")
  ######## ----------------------
  
  
  print("LOADING DESEQ2 results")
  dds_OvaryVsPri = get(load("../output/DESeq2/dds_OvaryVsPri.Rdata"))
  
  # FILTER OUT OUTLIERS
  OvaryVsPri <- processRes(results(dds_OvaryVsPri, contrast = c("TissueType","Normal","Primary")), 
                            lfcSE_filter = T, lfcSEcut = 2)
  
  print("ANNOTATING DESEQ2 results")
  # ANNOTATE WITH GENES
  OvaryVsPri = as.data.frame(annotate_DESeqOut(OvaryVsPri, geneFamilyMap))
  
  return(OvaryVsPri)
}