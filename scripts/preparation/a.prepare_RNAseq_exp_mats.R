load("../input/cts_ens.Rdata")
load("../input/cts_genes.Rdata")
load("../input/cts_genes_log2CPM.Rdata")
load("../input/log2_tpm_genes.Rdata")
load("../input/log2_tpm_ens.Rdata")


# LOAD GENEMAP
geneMap = get(load(file = "../input/bioMart_geneMap.Rdata")) # hugo symbols, gene names downloaded from biomart on 2021-Feb-04
pcGeneMap = subset(geneMap, gene_biotype == "protein_coding" & ensembl_gene_id %in% unique(c(rownames(cts))))

# LOAD GENE FAMILIES
geneFamilyMap = get(load(file = "../input/geneFamilyMap.Rdata")) # genename_DB.tsv downloaded from genenames website on 2021-Feb-04


### DO NOT RUN AGAIN -----------------------

# ## LOAD CLINICAL DATA
# clinData = read.csv("../input/Clinical_annotations.csv",stringsAsFactors = F)
# rownames(clinData) = clinData$SampleID
# 
# #### Load expression data - tximport
# txiOut = get(load("../input/txiOutput_bo_br_ov_gi_normal.rda"))
# cts = do.call("cbind", sapply(txiOut,`[[`,2)) # raw counts for DESeq
# dim(cts)
# tpm = do.call("cbind", sapply(txiOut,`[[`,1)) # length & library normalized TPM for visualization and other purpose
# dim(tpm)
# 
# # Convert Ensemble IDs to gene symbols
# pcGeneMap = subset(geneMap, gene_biotype == "protein_coding" &
#                     ensembl_gene_id %in% unique(c(rownames(cts), rownames(tpm))))
# 
# cts = cts[pcGeneMap$ensembl_gene_id,]; #rownames(cts)=pcGeneMap$external_gene_name
# dim(cts)
# tpm = tpm[pcGeneMap$ensembl_gene_id,]; #rownames(tpm)=pcGeneMap$external_gene_name
# dim(tpm)
# 
# # filter counts == 0
# cts = cts[rowSums(cts)>0,]
# dim(cts)
# tpm = tpm[rowSums(tpm)>0,]
# dim(tpm)
# 
# # filter TPM matrix - remove lowly exp genes
# tpmThreshold <- apply(tpm, MARG = 1, function(x) length(which(x > 0.5))) # tag genes > 0.5 tpm
# keep.index <- as.vector(which(tpmThreshold > 1)); tpm = tpm[keep.index,] # remove genes <= 0.5 tpm and not expressed at least 1 sample
# tpm = lib_size_scaling(tpm, 1e6) # renormalize transcripts per million
# colSums(tpm)[1:5]
# dim(tpm)
# 
# #### Subset for common samples
# commonSamples = as.character(na.omit(colnames(tpm)[match(rownames(clinData),colnames(tpm))])) # common samples
# 
# clinData = clinData[commonSamples,] # clinical data
# cts = cts[,commonSamples]
# dim(cts)
# tpm = tpm[,commonSamples]
# log2_tpm = as.matrix(log2(tpm+1)) # log2 tpm
# dim(tpm)
# 
# out_ = gene_name_conversion(rownames(log2_tpm), fromType = "ENSEMBL", toType = "SYMBOL")
# log2_tpm_genes = log2_tpm[out_$ENSEMBL,]
# rownames(log2_tpm_genes)=out_$SYMBOL
# dim(log2_tpm_genes)
# 
# save(log2_tpm, file = "../input/log2_tpm_ens.Rdata")
# save(log2_tpm_genes, file = "../input/log2_tpm_genes.Rdata")
# 
# out_ = gene_name_conversion(rownames(cts), fromType = "ENSEMBL", toType = "SYMBOL")
# cts_genes = cts[out_$ENSEMBL,]
# rownames(cts_genes) = out_$SYMBOL
# dim(cts_genes)
# 
# cts_CPM = getCPM(cts, logged = F)
# cts_log2CPM = getCPM(cts, logged = T)
# cts_genes_log2CPM = getCPM(cts_genes, logged = T)
# 
# dim(log2_tpm_genes)
# dim(cts_genes_log2CPM)
# 
# save(cts, file = "../input/cts_ens.Rdata")
# save(cts_genes, file = "../input/cts_genes.Rdata")
# save(cts_genes_log2CPM, file = "../input/cts_genes_log2CPM.Rdata")



#### --------------------



