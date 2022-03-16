makePAM50calls = function(log2_tpm_genes, clinData, er_status_df){
  
  # ER pos and neg samples
  er_pos_samples = subset(er_status_df, ER_infered %in% "Pos")$SampleID
  er_neg_samples = subset(er_status_df, ER_infered %in% "Neg")$SampleID

  # PAM50 calls
  library(genefu)
  data("pam50")
  data("pam50.robust")
  data("nkis")
  
  
  # pam50 genes 
  pam50.genes = rownames(pam50$centroids)
  pam50.genes[pam50.genes %notin% rownames(log2_tpm_genes)]
  pam50.genes[pam50.genes %in% "CDCA1"]="NUF2"
  pam50.genes[pam50.genes %in% "KNTC2"]="NDC80"
  pam50.genes[pam50.genes %in% "ORC6L"]="ORC6"
  
  # pam50 tpm
  pam50_tpm = t(log2_tpm_genes[pam50.genes, subset(clinData, TissueType %notin% "Normal")$SampleID])
  
  pam50.calls = as.data.frame(matrix(NA, nrow = nrow(pam50_tpm), ncol = 20))
  rownames(pam50.calls)=rownames(pam50_tpm)
  pam50.prob = data.frame(Basal = rep(NA,20), Her2 = rep(NA,20), LumA = rep(NA,20), LumB = rep(NA,20), Normal = rep(NA,20),  row.names = c(1:20))
  pam50.prob.df <- data.frame(Basal = rep(NA,nrow(pam50_tpm)), Her2 = rep(NA,nrow(pam50_tpm)), LumA = rep(NA,nrow(pam50_tpm)), LumB = rep(NA,nrow(pam50_tpm)),
                              Normal = rep(NA,nrow(pam50_tpm)),  row.names = rownames(pam50_tpm))
  
  for (i in row.names(pam50_tpm)) {
    for (j in c(1:20)) {
      er_pos_samples_test <- er_pos_samples[!er_pos_samples %in% i]
      er_neg_samples_test <- er_neg_samples[!er_neg_samples %in% i]
      if(identical(er_pos_samples_test, er_pos_samples)) {
        test.set <- c(sample(er_pos_samples_test,7),sample(er_neg_samples_test,6))
      }
      if(identical(er_neg_samples_test, er_neg_samples)) {
        test.set <- c(sample(er_pos_samples_test,6),sample(er_neg_samples_test,7))
      }
      sample.set <- c(test.set,i)
      pam50.sample.set.final <- pam50_tpm[sample.set,]
      PAM50.subtype <- molecular.subtyping(sbt.model="pam50", data = pam50.sample.set.final, annot = annot.nkis, do.mapping = FALSE)
      pam50.calls[i,j] <- as.character(PAM50.subtype$subtype[i])
      pam50.prob[j,] <- PAM50.subtype$subtype.proba[14,]
      print(paste(i, "PAM50 Test", j, "=>",as.character(PAM50.subtype$subtype[i])))
    }
    pam50.prob.df[i,] <- colMeans(pam50.prob)
  }
  
  Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  pam50.output.df <- as.data.frame(apply(X = pam50.calls, FUN = Mode, MARGIN = 1))
  colnames(pam50.output.df) <- "PAM50.subtype"
  
  icesTAF::mkdir("../output/annotations/PAM50")
  write.table(pam50.calls,
              file = "../output/annotations/PAM50/pam50_20foldTest.txt",
              quote = FALSE, col.names = NA, sep = "\t")
  write.table(pam50.output.df,
              file = "../output/annotations/PAM50/pam50_calls.txt",
              quote = FALSE, col.names = NA, sep = "\t")
  write.table(signif(pam50.prob.df,3),
              file = "../output/annotations/PAM50/pam50_probs.txt",
              quote = FALSE, col.names = NA, sep = "\t")
  return(list(pam50.output.df = pam50.output.df, pam50.calls = pam50.calls, pam50.prob.df = pam50.prob.df))
}