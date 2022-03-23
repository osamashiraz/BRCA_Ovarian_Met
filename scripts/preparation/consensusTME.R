load("../output/annotations/ConsensusTME/ConsensusTME_immunescores.Rdata")

### DO NOT RUN AGAIN -----------------------

## LOAD TPM VALUES
# load("../input/log2_tpm_genes.Rdata")
# rownames(log2_tpm_genes) = make.unique(rownames(log2_tpm_genes))
# 
# # devtools::install_github("cansysbio/ConsensusTME")
# library(ConsensusTME)
# ConsensusTME_out = consensusTMEAnalysis(log2_tpm_genes, cancer = "BRCA", statMethod = "ssgsea",
#                                         parallel.sz = 2) # GSVA is slow ~ 5-10min
# ConsensusTME_out = as.data.frame(ConsensusTME_out)
# sort(ConsensusTME_out["Immune_Score",])
# 
# icesTAF::mkdir("../output/annotations/ConsensusTME")
# save(ConsensusTME_out, file = "../output/annotations/ConsensusTME/ConsensusTME_immunescores.Rdata",quote = F)

### -----------------------