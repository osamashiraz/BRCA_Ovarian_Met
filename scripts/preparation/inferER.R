inferER = function(log2_tpm_genes, clinData){
  # ER status inferrence
  samples = subset(clinData, TissueType %notin% "Normal")$SampleID
  length(samples)
  er.status = clinData[colnames(log2_tpm_genes[, samples]), "ER"] # no normals
  
  er_status_df = data.frame(ER = er.status, esr1 = log2_tpm_genes["ESR1", samples], SampleID = clinData[colnames(log2_tpm_genes[, samples]),"SampleID"],
                            patientID = clinData[colnames(log2_tpm_genes[, samples]),"Patient_ID"],tissueType = clinData[colnames(log2_tpm_genes[, samples]),"TissueType"],
                            site = clinData[colnames(log2_tpm_genes[, samples]),"TissueSite"], study = clinData[colnames(log2_tpm_genes[, samples]),"Study_Code"])
  er_status_df$q3 = er_status_df$esr1 > quantile(er_status_df$esr1)[2]
  
  # ggplot(er_status_df, aes(x = reorder(SampleID,esr1), y = esr1, color = ER)) + geom_point(size = 5, aes(label = SampleID)) + myTheme +
  #   theme(axis.text.x = element_text(angle = 90, size = 6))+xlab("")+
  #   geom_hline(yintercept = quantile(er_status_df$esr1)[2])+facet_wrap(.~site, nrow = 3)
  
  ## second quantile --> median as cutoff for ER+ and ER-
  er_status_df$ER_infered = ifelse(er_status_df$esr1 > quantile(er_status_df$esr1)[2], "Pos", "Neg")
  
  # ggplot(er_status_df, aes(x = reorder(SampleID, esr1), y = esr1, color = ER_infered))+geom_point(size = 5)+myTheme+
  #   theme(axis.text.x = element_text(angle = 90))+xlab("")+
  #   geom_hline(yintercept = quantile(er_status_df$esr1)[2]) +
  #   ggrepel::geom_label_repel(nudge_x = 2, nudge_y = 3, data = subset(er_status_df, SampleID %in% c("OM6","0001-56M-01","0014-56M","OM7","0008-56M","0006-56M")),
  #                             aes(label = patientID))
  # 
# subset(er_status_df, patientID %in% c("Ov13","Ov12"))
  return(er_status_df)
}