# LOAD PACKAGES
library(magrittr)
library(ggplot2)
library(cowplot)
library(dplyr)
library(EnvStats)
library(ggpubr)
library(gridExtra)
library(ggthemes)
library(reshape2)

`%notin%`=Negate(`%in%`)
options("scipen"=100, "digits"=4)

## THEMEING
blankTheme = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), 
                   panel.background = element_blank(), axis.ticks.x = element_blank(),plot.margin = grid::unit(c(0,0,0,0), "mm")) 
myTheme = theme_bw()+theme(# TEXT
  text = element_text(size = 15), 
  axis.text = element_text(size = 12, face = "bold", color = "#1A237E"), 
  legend.text = element_text(size = 12, face = "bold"),
  plot.subtitle = element_text(size = 8,face="bold"),
  # BACKGROUND
  legend.background = element_rect(fill = "white", size = 4, colour = "white"),
  panel.grid.major = element_blank(),#(color = "#9E9E9E",linetype = "dotted", size = 0.1),
  panel.grid.minor = element_line(color="#BDBDBD",linetype = "dotted"),
  # TICKS
  #LEGEND POS
)
myTheme_barplot = theme_bw()+theme(# TEXT
  text = element_text(size = 15, family = "Helvetica",face = "bold"), 
  legend.text = element_text(size = 10,face="bold"),
  legend.title = element_text(size = 12,face="bold"),
  rect = element_blank(),
  #rect = element_rect(linetype = "dotted",colour = "grey"),
  # axis.title = element_text(size = 20),
  axis.line.y = element_line(colour = "#424242",size = 0.5),
  #axis.line.x = element_blank(),
  plot.title = element_text(hjust = 0.5),
  # line = element_blank(),
  panel.background = element_rect(fill =  "white"),
  panel.grid = element_blank(),
  axis.title = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank()
  # TICKS
  #LEGEND POS
)
myTheme_umap = theme_bw()+theme(# TEXT
  text = element_text(size = 15), 
  legend.text = element_text(size = 10,face="bold"),
  legend.title = element_text(size = 12,face="bold"),
  rect = element_blank(),
  # rect = element_rect(linetype = "dotted",colour = "grey"),
  # axis.title = element_text(size = 20),
  axis.line = element_line(linetype = "dashed",colour = "gray",#"#0091EA",
                           size = 0.5,
                           arrow = ggplot2::arrow(angle = 30,length = unit(0.24, "cm"), type = "closed")),
  plot.title = element_text(hjust = 0.5),
  # line = element_blank(),
  # panel.background = element_rect(colour = "white",color = "white"),
  panel.grid = element_blank(),
  axis.title = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank()
  # TICKS
  #LEGEND POS
)
myEnrich_barplot = function(df, order = T, order_by="logpval", plot_top_n= 10, plot_fdr = T, 
                            abr = F, abr_len = 50, title="", abline_val = 0.05,
                            title_color="#0097A7", font_color="black", bar_color="#80DEEA", add_cutoff = F){
  # browser()
  if (sum(names(df) %in% "label") == 0){ # cluster profiler output does not have label - change names for cluster profiler compatibility
    df$label = df$Description
    df$pval = df$pvalue
    df$fdr = df$p.adjust
  }
  if (abr){
    df$abr_label = abbreviate(df[,"label"], minlength = abr_len)
  }else{
    df[,"abr_label"]=df[,"label"]
  }
  if (plot_fdr){
    df[,"logpval"]=-log10(df[,"fdr"])
  }else{
    df[,"logpval"]=-log10(df[,"pval"])
  }
  if (order){
    df = df[order(df[,order_by], decreasing = T),]
  }
  # browser()
  df_sub = na.omit(df[1:plot_top_n,])
  hide_NA = NULL
  if (nrow(df_sub) < plot_top_n){
    pad_length = plot_top_n-nrow(df_sub)
    padding = as.data.frame(matrix(data = NA, nrow = pad_length, ncol = ncol(df_sub)))
    colnames(padding)=colnames(df_sub)
    padding$logpval = 0
    padding$abr_label = as.character(1:pad_length)
    hide_NA = theme(axis.text.x = element_text(colour = c(rep("#424242", nrow(df_sub)), rep("#FAFAFA", pad_length))))
    df_sub = rbind(df_sub, padding)
  }
  if(add_cutoff == T){
    geomVline = geom_vline(xintercept = -log10(abline_val), linetype='dashed', color="red", alpha = 0.6, size = 1)
  }else{
    geomVline = NULL
  }
  plt = ggplot(df_sub, aes(y = reorder(abr_label, logpval), x = logpval, label = label)) +geomVline+
    geom_bar(stat = "identity", fill = bar_color, alpha = 0.5,width = 0.75) + myTheme_barplot + ggtitle(title) +
    geom_text(aes(x = 0.1), size = 4.4, face="bold", hjust = 0, vjust = 0.35, alpha = 0.6, color = font_color) + 
    scale_x_continuous(expand = c(0, 0.2))+theme(title = element_text(size = 12, color = title_color))+
    scale_x_continuous(labels = function(x) sprintf("%.1f", x))+hide_NA #scale_y_discrete(limits = rev(levels(reorder(df$abr_label, dflogpval))))
  return(list(plt = plt,df = df))
}


# gene name conversion
gene_name_conversion = function(geneNames, fromType = "SYMBOL", toType = "ENTREZID"){
  load("Y:/02_Datasets/Resources/bioMart_geneMap.Rdata")
  library(org.Hs.eg.db)
  library(clusterProfiler)
  map = bitr(geneID = geneNames, fromType = fromType, toType = toType, OrgDb = org.Hs.eg.db)
  map = map[!duplicated(map[,1]),]
  rownames(map)=map[,1]
  geneNames = intersect(map[,1],geneNames)
  return(map)
}

# Process DESeq results
processRes = function(res, lfcSE_filter = F, lfcSEcut = 3){
  res = na.omit(res)
  res$absLogFC = abs(res$log2FoldChange)
  res = res[order(res$padj),]  
  res$padj[res$padj == 0]=min(res$padj[res$padj != 0])
  res$log10Padj=-log10(res$padj)
  res$log10Pvalue=-log10(res$pvalue)
  if(lfcSE_filter){
    res = subset(res, lfcSE < lfcSEcut) 
  }
  return(res)
}


## COMPUTE
getCPM = function(counts, logged = F){
  y <- edgeR::DGEList(counts)
  y <- edgeR::calcNormFactors(y)
  CPM <- edgeR::cpm(y, normalized.lib.sizes = TRUE)
  if(logged){CPM = log2(CPM+1)}
  return(CPM)
}

lib_size_scaling = function(mat, lib.scale = 1e6){
  as.data.frame(apply(mat, 2, function(x) lib.scale*x/sum(x)))
}




## Qorts QC functions
getColors = function(x){
  c(RColorBrewer::brewer.pal(9, name = "Set1"),RColorBrewer::brewer.pal(8, name = "Set2"))[1:x]
}
myBasicPlotter = function(res, samples){
  resSub = res
  resSub@lanebam.list = resSub@lanebam.list[resSub@lanebam.list %in% samples]
  resSub@sample.list = resSub@sample.list[resSub@sample.list %in% samples]
  resSub@decoder = resSub@decoder[samples,]
  resSub@summaryTable = resSub@summaryTable[,samples]
  resSub@calc.data = lapply(resSub@calc.data, function(x) {
    if (is.data.frame(x)){
      return(x[samples,])
    }
    if (is.list(x)){
      return(x[samples])
    }
  })
  resSub@qc.data = lapply(resSub@qc.data, function(x) {
    if (is.data.frame(x)){
      return(x[samples,])
    }
    if (is.list(x)){
      return(x[samples])
    }
  })
  basic.plotter = build.plotter.colorBySample(resSub, plotter.params = list("std.lines.alpha"=255))
  # n = length(basic.plotter$lanebam.params$lines.col)
  # basic.plotter$lanebam.params$lines.col = getColors(n)
  # basic.plotter$lanebam.params$points.col = getColors(n)
  # basic.plotter$lanebam.params$points.tcol = getColors(n)
  # basic.plotter$lanebam.params$lines.tcol = getColors(n)
  # basic.plotter$legend.params$lines.col = getColors(n)
  return(basic.plotter)
}



## Score Data frame
mkDF_scores = function(scoreMat, signature, clinData){
  # browser()
  IDs = intersect(colnames(scoreMat),clinData$SampleID)
  clinData = clinData[IDs,]
  data.frame(
    score = as.numeric(scoreMat[signature, clinData$SampleID]), 
    SampleIDs= clinData$SampleID, 
    Subtype = clinData$Subtype, 
    Group = clinData$Patient_ID,
    TissueType = clinData$TissueType,#relevel(factor(clinData$TissueType),ref = "Primary"), 
    Site = clinData$TissueSite, 
    ER_infered = clinData$ER_infered, 
    ER = clinData$ER, 
    PR = clinData$PR,
    HER2 = clinData$HER2,
    Study = clinData$Study_Code,
    PAM50 = clinData$PAM50,
    Paired = clinData$Paired,
    Purity = clinData$Estimate_Purity,
    stringsAsFactors = F
  )
}


