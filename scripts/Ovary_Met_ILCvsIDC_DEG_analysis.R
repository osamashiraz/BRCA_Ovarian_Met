run_DEseq2_Ovary_Met_ILCvIDC = function(cts, clinData, excluded_samples){
  library(DESeq2)
  
  mode(cts)="integer"
  icesTAF::mkdir("../output/deseq2")
  
  
  # ILC Met vs IDC Met
  clinData_select = subset(clinData, SampleID %notin% excluded_samples & Subtype %in% c("ILC","IDC") & 
                                 TissueType %in% "Metastasis" & Study_Code %in% "Ovary")
  table(clinData_select$Subtype)
  
  # OM5B, 0001-56M-02 # remove these duplicated
  clinData_select = subset(clinData_select, SampleID %notin% c("OM5B", "0001-56M-02"))
  table(clinData_select$Subtype)
  length(clinData_select$SampleID)
  
  
  #### DEG Analysis ----------------- DO NOT RUN AGAIN
  # library(DESeq2)
  # 
  # dds_ILCvsIDC_mets <- DESeqDataSetFromMatrix(countData = cts[, clinData_select$SampleID],
  #                                                    colData = clinData_select, design = ~Subtype)
  # dds_ILCvsIDC_mets <- DESeq(dds_ILCvsIDC_mets, parallel = T)
  # 
  # save(dds_ILCvsIDC_mets,file = "../output/deseq2/dds_Ovary_Met_ILCvsIDC.Rdata")
  # write.dcf(data.frame(design="= ~Subtype"),file = "../output/deseq2/dds_Ovary_Met_ILCvsIDC Run Info.txt")
  ######### ----------------------
  
  
  print("LOADING DESEQ2 results")
  dds_ILCvsIDC_mets = get(load("../output/deseq2/dds_Ovary_Met_ILCvsIDC.Rdata"))
  
  # FILTER OUT OUTLIERS
  OvILCvsIDC_met <- processRes(results(dds_ILCvsIDC_mets, contrast = c("Subtype","ILC","IDC")),
                               lfcSE_filter = T, lfcSEcut = 2)
  
  print("ANNOTATING DESEQ2 results")
  # ANNOTATE WITH GENES
  OvILCvsIDC_met = as.data.frame(annotate_DESeqOut(OvILCvsIDC_met, geneFamilyMap))
  
  return(OvILCvsIDC_met)
}







