# load clinical data (remove NOT sure) - add PAM50, Tumor purity, Immune
#  Fig 1A - Clinical data pie chart
#  Sup Fig 1 - Met site and histology association


prepareAnnotations = function(){
  
  # Load clinical annotations (merged annotations for ovary metastasis cohort and breastMET from L. Zhu et al 2019, D. Vareslija et al 2019 and K. M. Levine et al 2019)
  clinData = read.csv("../input/Clinical_annotations.csv",stringsAsFactors = F)
  rownames(clinData) = clinData$SampleID
  table(clinData$Center)
  
  load("../input/log2_tpm_genes.Rdata")
  
  
  print("LOADING TUMOR PURITY")
  
  # 1. ESTIMATE
  source("./preparation/estimate_tumorPurity.R")
  clinData$Estimate_Purity = tumorPurityScores[4, clinData$SampleID] %>% as.numeric()
  
  print("LOADING IMMUNE SCORES")
  
  # 2. ConsensusTME
  source("./preparation/consensusTME.R")
  clinData$Immune_Scores = as.numeric(ConsensusTME_out["Immune_Score", clinData$SampleID])
  
  #plot(ConsensusTME_out["Immune_Score", colnames(tumorPurityScores)] %>% as.numeric(), tumorPurityScores[4, ])
  
  print("INFERRING ER STATUS FROM mRNA")
  
  # 3. ER inferred from mRNA
  source("./preparation/inferER.R")
  # er_status_df = inferER(log2_tpm_genes, clinData) # do not run again
  load("../output/annotations/ER_infered/er_status_df.Rdata")
  table(er_status_df$ER, er_status_df$ER_infered)
  dim(er_status_df)
  
  clinData$ER_infered = er_status_df[clinData$SampleID, "ER_infered"]
  
  print("LOADING PAM50 CALLS")
  
  ## 4. PAM50
  source("./preparation/PAM50.R")
  # makePAM50calls(log2_tpm_genes, clinData, er_status_df) # do not run again
  
  ## add PAM50 to clinical data
  pam50.output.df = read.delim("../output/annotations/PAM50/pam50_calls.txt")
  rownames(pam50.output.df)=pam50.output.df$X
  
  clinData$PAM50 = pam50.output.df[clinData$SampleID,"PAM50.subtype"]
  
  return(clinData)
}
