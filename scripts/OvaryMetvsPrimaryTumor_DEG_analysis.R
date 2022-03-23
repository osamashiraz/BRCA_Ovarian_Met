run_DEseq2_OvaryMetvsPrimaryTumor = function(cts, clinData, excluded_samples){
  library(DESeq2)
  
  mode(cts)="integer"
  icesTAF::mkdir("../output/deseq2")
  
  
  # ILC Met vs IDC Met
  clinData_select = subset(clinData, Subtype %in% c("IDC","ILC","mDLC") & 
                             Study_Code %in% "Ovary" &  SampleID %notin% excluded_samples)
  table(clinData_select$TissueType)
  

  #### DEG Analysis ----------------- DO NOT RUN AGAIN
  # library(DESeq2)
  # 
  # dds_OvaMetvsPri <- DESeqDataSetFromMatrix(countData = cts[, clinData_select$SampleID],
  #                                                    colData = clinData_select, design = ~TissueType+Patient_ID)
  # dds_OvaMetvsPri <- DESeq(dds_OvaMetvsPri, parallel = T)
  # 
  # save(dds_OvaMetvsPri,file = "../output/deseq2/dds_OvaMetvsPri.Rdata")
  # write.dcf(data.frame(design="= ~TissueType+Patient_ID"),file = "../output/deseq2/dds_OvaMetvsPri Run Info.txt")
  ######### ----------------------
  
  
  print("LOADING DESEQ2 results")
  dds_OvaMetvsPri = get(load("../output/DESeq2/dds_OvaMetvsPri.Rdata"))
  
  # FILTER OUT OUTLIERS
  OvaMetvsPri <- processRes(results(dds_OvaMetvsPri, contrast = c("TissueType","Metastasis","Primary")), 
                            lfcSE_filter = T, lfcSEcut = 2)
  
  print("ANNOTATING DESEQ2 results")
  # ANNOTATE WITH GENES
  OvaMetvsPri = as.data.frame(annotate_DESeqOut(OvaMetvsPri, geneFamilyMap))
  
  return(OvaMetvsPri)
}