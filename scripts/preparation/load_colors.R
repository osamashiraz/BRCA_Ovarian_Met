## COLOR PALLETS
print("LOAD COLORS")

pam50_cols = RColorBrewer::brewer.pal(n = 5, name = "Set1")[c(1,2,5,3,4)]
names(pam50_cols)=c("LumA","LumB","Basal","Her2","Normal")

# pgr_stat = summary(log2_tpm_genes["PGR",selected_clinData$SampleID])
#esr1_stat = summary(as.numeric((log2_tpm_genes["ESR1",])))
#erbb2_stat = summary(as.numeric((log2_tpm_genes["ERBB2",])))
# ar_stat = summary(log2_tpm_genes["AR",selected_clinData$SampleID])
annot_col = list(Subtype = c("IDC"="blue","ILC"="red",'mDLC'="green","Normal"="gray","Adeno"="#FF9800"),
                 ER = c("Pos"="#F57F17","Neg"="#212121"),
                 `Tissue Type` = c("Primary"="#00C853","Metastasis"="#304FFE","Normal"="gray"),
                 Study_Code = c("Brain"="#9C27B0","Bone"="#2196F3","GI"="#43A047","Ovary"="#FFEB3B","Breast"="#8BC34A"),
                 Site = c("Brain"="#9C27B0","Bone"="#2196F3","GI"="#43A047","Ovary"="#FFEB3B","Breast"="#8BC34A","NOT SURE"="gray"),
                 PAM50 = pam50_cols,
                 `Metastasis Type` = c("Orphan"="#BDBDBD","Paired"="#E91E63"), #c("FALSE"="#212121","TRUE"="#F9A825"),
                 Purity = circlize::colorRamp2(c(0,0.4, 0.8), c("#1B5E20", "black", "#f44336")),
                 # ESR1 = circlize::colorRamp2(c(esr1_stat[1],esr1_stat[2],esr1_stat[6]), c("#1B5E20", "black", "#f44336")),
                 # HER2 = circlize::colorRamp2(c(erbb2_stat[1],erbb2_stat[2],erbb2_stat[6]), c("#1B5E20", "black", "#f44336")),
                 Immune = circlize::colorRamp2(c(-0.2, 0, 0.2), c("#1B5E20", "black", "#f44336"))
)

