load("../output/annotations/ESTIMATE/ESTIMATE_scores_final.RData")
sort(tumorPurityScores[4,])

### DO NOT RUN AGAIN -----------------------

# 
# # # LOAD ANNOTATIONS
# clinData = read.csv("../input/Clinical_annotations.csv",stringsAsFactors = F)
# rownames(clinData) = clinData$SampleID
# 
# # LOAD TPM VALUES
# load("../input/log2_tpm_genes.Rdata")
# rownames(log2_tpm_genes) = make.unique(rownames(log2_tpm_genes))
# 
# 
# library(estimate)
# 
# ##### Estimate purity
# icesTAF::mkdir("../output/annotations/ESTIMATE/")
# write.table(log2_tpm_genes, file="../output/annotations/ESTIMATE/input.txt", quote=FALSE,	sep="\t")
# filterCommonGenes(input.f="../output/annotations/ESTIMATE/input.txt",
#                   output.f="../output/annotations/ESTIMATE/input.gct",
#                   id="GeneSymbol")
# 
# # no purity estimate for illumina paltform
# estimateScore("../output/annotations/ESTIMATE/input.gct",
#               "../output/annotations/ESTIMATE/ESTIMATE_scores.gct",
#               platform="illumina")
# 
# 
# #### read in tumor purity
# library("CePa")
# estimate_out_illumina = read.gct("../output/annotations/ESTIMATE/ESTIMATE_scores.gct")
# ids_illumina = gsub(x=colnames(estimate_out_illumina),pattern = "X","") # CORRECT COLNAMES
# ids_illumina = gsub(x=ids_illumina, pattern = "[.]","-")
# colnames(estimate_out_illumina)=ids_illumina
# 
# # COMPUTE PURITY
# get_purity <- function(est_score){
#   # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3826632/
#   return(cos(0.6049872018 + 0.0001467884 * est_score))
# }
# 
# rwNames=rownames(estimate_out_illumina)
# estimate_out_illumina=rbind(estimate_out_illumina,get_purity(estimate_out_illumina[3,]))
# rownames(estimate_out_illumina)=c(rwNames,"Tumor Purity")
# 
# # SAVE ESTIMATE OUT
# tumorPurityScores = data.matrix(estimate_out_illumina)
# save(tumorPurityScores, file= "../output/annotations/ESTIMATE/ESTIMATE_scores_final.RData")



###  -----------------------
