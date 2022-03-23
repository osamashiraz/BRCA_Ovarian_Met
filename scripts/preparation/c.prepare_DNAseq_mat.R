
load_filter_SNVs = function(clinData){
  load("../input/MammaSeq_PON_variants.Rdata")
  load("../input/MammaSeq_SNVs_Annotated.Rdata")
  
  # Profiled Samples
  commonSamples = intersect(clinData$SampleID, Mamma$Sample.ID %>% unique())
  
  
  # All profiled samples
  Mamma_all = subset(Mamma, Sample.ID %in% clinData[commonSamples,]$SampleID)
  length(table(Mamma_all$Sample.ID))
  
  # Filtering
  dim(Mamma_all) #2158   
  Mamma_sub = Mamma_all #Mamma_paired # Mamma_all #Mamma_paired
  
  # FILTER 0: REMOVE VARIANTS IN BUFFY CONTROLS
  Mamma_sub$alteration = paste0(Mamma_sub$HUGO.symbol,".",Mamma_sub$Protein.sequence.change)
  Mamma_sub = subset(Mamma_sub, alteration %notin% exclude_variants)
  dim(Mamma_sub) # 1145   
  
  # FILTER 1: REMOVE GENOMIC VARIANTS THAT ARE PRESENT IN 30% OR MORE SAMPLES (LIKELY TO BE SEQUENCING ARTIFACTS OR COMMON VARIATIONS)
  mamma.filter.1 = Mamma_sub[Mamma_sub$Position %in% names(which(table(Mamma_sub$Position) < length(unique(Mamma_sub$Sample.ID))*20/100)),]
  dim(mamma.filter.1)# 757
  
  # FILTER 2: REMOVE ALL SYNONYMOUS VARIANTS (SY)
  mamma.filter.2 = mamma.filter.1[-grep('SY',mamma.filter.1$Sequence.ontology),]
  dim(mamma.filter.2) # 475
  
  # FILTER 3: REMOVE LOW CONFIDENCE VARIANTS (SELECT FOR VARIANTS WITH STRAND BIAS BETWEEN 0.5 and 0.6)
  # STB FORMULA: READS IN FORWARD STRAND / (READS IN FORWARD STRAND + READS IN REVERSE STRAND), 0.6 >= x >= 0.5 
  # MEANS WE DISREGARD REGIONS WHERE REVERSE STRAND READS ARE HIGHER THAN FORWARD 
  # BUT WE INCLUDE VARIANTS IF THE FOWARD STRAND READS ARE LITTLE BITT MORE THAN REVERSE ONES, REASONING
  mamma.filter.3 = mamma.filter.2[mamma.filter.2$STB >= 0.4 & mamma.filter.2$STB <= 0.6, ]
  dim(mamma.filter.3) # 356
  
  # FILTER 4: REMOVE VARIANTS ANNOTATED IN GNOMAD AND 1000GENOME WITH AF >= 1%
  mamma.filter.3$gnomAD.AF.Total = as.numeric(mamma.filter.3$gnomAD.AF.Total); mamma.filter.3$X1000.Genomes.AF = as.numeric(mamma.filter.3$X1000.Genomes.AF)
  mamma.filter.4 = mamma.filter.3[mamma.filter.3$gnomAD.AF.Total < 0.01 | mamma.filter.3$X1000.Genomes.AF < 0.01, ]
  dim(mamma.filter.4) # 349
  
  # FILTER 5: REMOVE VARIANTS WITH AF > 90%, THESE ARE LIKELY GERMLINE
  mamma.filter.4$AF = as.numeric(mamma.filter.4$AF)
  mamma.filter.5 = mamma.filter.4[mamma.filter.4$AF < 0.9, ]
  dim(mamma.filter.5) # 344
  
  # FILTER 6: READ DEPTH
  mamma.final = subset(mamma.filter.5, Total.reads > 3) # Sequencing depth filter
  dim(mamma.final) # 344
  
  return(list(unfiltered_SNVs = Mamma_all, filtered_SNVs = mamma.final))
}


prepare_SNV_matrix = function(mamma.final, clinData){
  MUT = mamma.final[, c("Sample.ID", "HUGO.symbol", "Sequence.ontology")]
  MUT$HUGO.symbol <- factor(MUT$HUGO.symbol)
  
  MUT = subset(MUT, Sample.ID %in% unique(clinData$SampleID))
  table(clinData[unique(MUT$Sample.ID), "TissueType"])
  sort(unique(MUT$Sample.ID))
  
  MUT_SS = subset(MUT,Sequence.ontology %in% "SS")
  MUT_MS = subset(MUT,Sequence.ontology %in% "MS")
  MUT_SG = subset(MUT,Sequence.ontology %in% "SG")
  
  MUT_SS$Sequence.ontology = rep(1,nrow(MUT_SS))
  MUT_SG$Sequence.ontology = rep(1,nrow(MUT_SG))
  MUT_MS$Sequence.ontology = rep(1,nrow(MUT_MS))
  
  if (nrow(MUT_SS) > 0){
    MUT_SS = dcast(MUT_SS, HUGO.symbol~Sample.ID, value.var = "Sequence.ontology",fun.aggregate = sum)
  }else{
    MUT_SS = NULL
  }
  if (nrow(MUT_SG) > 0){
    MUT_SG = dcast(MUT_SG, HUGO.symbol~Sample.ID, value.var = "Sequence.ontology",fun.aggregate = sum)
  }else{
    MUT_SG = NULL
  }
  MUT_MS = dcast(MUT_MS, HUGO.symbol~Sample.ID, value.var = "Sequence.ontology",fun.aggregate = sum)
  
  mat.list <- list(MUT_SS = MUT_SS, MUT_SG = MUT_SG, MUT_MS = MUT_MS)#, AMP = NULL, DEL = NULL)
  head(mat.list)
  
  mat.list = lapply(mat.list, function(x){
    if (!is.null(nrow(x))){
      rownames(x)=x[,1]
      x = x[,-1]
      x = as.matrix(x)
    }else{
      return(NULL)
    }
  })
  library(ComplexHeatmap)
  mat.list.uni <- unify_mat_list(mat.list)
  
  names(mat.list.uni) = c("Splice Site","Stop Gain","Mis-sense")
  return(mat.list.uni)
}

prepare_SNV_matrix_proteinChange = function(mamma.final, clinData){
  MUT = mamma.final[, c("Sample.ID", "alteration", "Sequence.ontology")]
  MUT$alteration <- factor(MUT$alteration)
  
  MUT = subset(MUT, Sample.ID %in% unique(clinData$SampleID))
  table(clinData[unique(MUT$Sample.ID), "TissueType"])
  sort(unique(MUT$Sample.ID))
  
  MUT_SS = subset(MUT,Sequence.ontology %in% "SS")
  MUT_MS = subset(MUT,Sequence.ontology %in% "MS")
  MUT_SG = subset(MUT,Sequence.ontology %in% "SG")
  
  MUT_SS$Sequence.ontology = rep(1,nrow(MUT_SS))
  MUT_SG$Sequence.ontology = rep(1,nrow(MUT_SG))
  MUT_MS$Sequence.ontology = rep(1,nrow(MUT_MS))
  
  if (nrow(MUT_SS) > 0){
    MUT_SS = dcast(MUT_SS, alteration~Sample.ID, value.var = "Sequence.ontology",fun.aggregate = sum)
  }else{
    MUT_SS = NULL
  }
  if (nrow(MUT_SG) > 0){
    MUT_SG = dcast(MUT_SG, alteration~Sample.ID, value.var = "Sequence.ontology",fun.aggregate = sum)
  }else{
    MUT_SG = NULL
  }
  MUT_MS = dcast(MUT_MS, alteration~Sample.ID, value.var = "Sequence.ontology",fun.aggregate = sum)
  
  mat.list <- list(MUT_SS = MUT_SS, MUT_SG = MUT_SG, MUT_MS = MUT_MS)#, AMP = NULL, DEL = NULL)
  head(mat.list)
  
  mat.list = lapply(mat.list, function(x){
    if (!is.null(nrow(x))){
      rownames(x)=x[,1]
      x = x[,-1]
      x = as.matrix(x)
    }else{
      return(NULL)
    }
  })
  library(ComplexHeatmap)
  mat.list.uni <- unify_mat_list(mat.list)
  
  names(mat.list.uni) = c("Splice Site","Stop Gain","Mis-sense")
  return(mat.list.uni)
}


generate_oncoprint = function(split_by = c("Subtype","Patient_ID"), clinData, mat.list.uni, n_genes =5){
  library("ComplexHeatmap")

  # PREPARING VISUALS
  col = c("Stop Gain" = "#b71c1c", "Splice Site" = "#03A9F4", "Mis-sense" = "#8BC34A")
  alter_fun = list(
    background = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = 'NA')) },
    `Mis-sense` = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#8BC34A", col = NA)) },
    `Splice Site` = function(x, y, w, h) {
      grid.rect(x, y, w, h, gp = gpar(fill = NA,col="#03A9F4", lwd = 2))},
    `Stop Gain` = function(x, y, w, h) {
      grid.points(x, y, pch = 16, gp = gpar(col="#b71c1c"))}
  )

  sampleIDS = intersect(colnames(mat.list.uni$`Mis-sense`), clinData %>% rownames())
  clinData_sub = clinData[sampleIDS,]
  
  htanno = HeatmapAnnotation(Subtype = factor(clinData_sub[,"Subtype"]),
                             ER = clinData_sub[,"ER_infered"],
                             PAM50 = clinData_sub[,"PAM50"],
                             annotation_name_gp = gpar(fontsize = 10, fontface="bold", col="#616161"),
                             na_col = "white", annotation_height = unit(0.05, "mm"), col = annot_col)
  
  column_ordeR = order(clinData_sub$TissueType, decreasing = T)
  
  # Select top N genes
  select_genes = head(names(rev(rowSums(Reduce(`+`,mat.list.uni) > 0) %>% sort())), n_genes) # genes altered in many samples
  matt = lapply(mat.list.uni, function(x) x[rownames(x) %in% select_genes, clinData_sub$SampleID])
  
  clinData_sub$SimpleID = gsub(x = clinData_sub$SimpleID, pattern = "M_1", replacement = "M1")
  clinData_sub$SimpleID = gsub(x = clinData_sub$SimpleID, pattern = "M_2", replacement = "M2")
  clinData_sub$SimpleID = gsub(x = clinData_sub$SimpleID, pattern = "_", replacement = " ")
  
  columnSplit = clinData_sub[,split_by]
  
  tumor.ht = oncoPrint(matt, alter_fun = alter_fun, col = col, remove_empty_columns = F,bottom_annotation = htanno,
                       row_names_gp = gpar(fontsize = 10, fontface="bold"),
                       column_names_gp = gpar(fontsize = 12, fontface="bold"), column_labels = clinData_sub$SimpleID,
                       column_title = NULL, show_column_names = T,remove_empty_rows = F, column_order = column_ordeR,#),
                       column_split = columnSplit)
  return(tumor.ht)
}


mutational_changes = function(snv_list, clinData){
  mamma.final = snv_list$filtered_SNVs
  
  sampleSpace = intersect(clinData$SampleID, unique(snv_list$unfiltered_SNVs$Sample.ID))
  priSamples = subset(clinData[sampleSpace,], TissueType %in% "Primary")$SampleID
  metSamples = subset(clinData[sampleSpace,], TissueType %in% "Metastasis")$SampleID
  
  MUT = mamma.final[, c("Sample.ID", "HUGO.symbol", "Sequence.ontology", "Protein.sequence.change")]
  MUT_sub = subset(MUT, Sample.ID %in% sampleSpace)
  
  #subset(MUT_sub, Sample.ID %in% priSamples)$Sample.ID %>% unique() %>% length()
  #subset(MUT_sub, Sample.ID %in% metSamples)$Sample.ID %>% unique() %>% length()
  
  library(dplyr)
  library(magrittr)
  
  pri = subset(MUT_sub, Sample.ID %in% priSamples)
  pri$patientID = clinData[pri$Sample.ID, "Patient_ID"]
  
  #unique(pri$patientID) %>% length() ### multiple tissue samples from same patient combined
  
  pri$ID = paste0(pri$patientID,pri$HUGO.symbol)
  pri = pri[!duplicated(pri$ID),]
  pri = pri %>% dplyr::group_by(HUGO.symbol) %>% dplyr::summarise(n = n()) %>% mutate(Percentage = 100*n/length(clinData[priSamples, "Patient_ID"] %>% unique())) %>% as.data.frame()
  rownames(pri)=pri$HUGO.symbol
  
  met = subset(MUT_sub, Sample.ID %in% metSamples)
  met$patientID = clinData[met$Sample.ID, "Patient_ID"]
  
  #unique(met$patientID) %>% length()
  
  met$ID = paste0(met$patientID, met$HUGO.symbol)
  met = met[!duplicated(met$ID),]
  met = met %>% dplyr::group_by(HUGO.symbol) %>% dplyr::summarise(n = n()) %>% mutate(Percentage = 100*n/length(clinData[metSamples, "Patient_ID"] %>% unique())) %>% as.data.frame()
  rownames(met)=met$HUGO.symbol
  
  allGenes = unique(c(met$HUGO.symbol %>% as.character(), pri$HUGO.symbol %>% as.character()))
  
  df = data.frame(genes = allGenes, 
                  pri_pct = pri[allGenes,"Percentage"], 
                  met_pct = met[allGenes,"Percentage"])
  df[is.na(df$pri_pct),"pri_pct"]=0
  df[is.na(df$met_pct),"met_pct"]=0
  
  names(df)[2:3]=c("Primary","Metastasis")
  df$change = df$Metastasis-df$Primary
  
  df_change = df[order(abs(df$change), decreasing = T),]
  rownames(df_change)=df_change$genes
  
  dn_mutations = na.omit(subset(df_change, change < 1)[1:5,1:3])
  up_mutations = na.omit(subset(df_change, change > -1)[1:5,1:3])
  df_melt = reshape2::melt(rbind(dn_mutations,up_mutations))
  df_melt$change_type = ifelse(df_change[df_melt$genes,"change"] > 0, "Met Up","Pri Up")
  df_melt$change = df_change[df_melt$genes,"change"]
  
  df_melt$genes = factor(df_melt$genes, levels = unique(df_melt$genes)[order(df_melt$change)])
  df_melt$value = round(df_melt$value, digits = 1)
  
  return(df_melt)
}



clinically_significant_variants = function(mamma.final) {
  load("../input/2020-02-20-PMKB.Rdata")
  load("../input/2020-02-20-oncokb.Rdata")
  
  # ADD IDENTIFIER
  oncokb$Alteration = gsub(x = oncokb$Alteration, pattern = "Amplification", replacement = "Gain")
  oncokb$Alteration = gsub(x = oncokb$Alteration, pattern = "Deletion", replacement = "Loss")
  oncokb$combined = as.character(interaction(oncokb$Gene, oncokb$Alteration))
  mamma.final$combined = as.character(interaction(mamma.final$HUGO.symbol,mamma.final$Protein.sequence.change))
  length(table(mamma.final$Sample.ID))
  
  cnvCombined = NULL
  # CREATE DATA FRAME
  combinedAlterations = data.frame(Sample.ID = rbind(cbind(as.character(mamma.final$Sample.ID)),
                                                     cbind(cnvCombined$Sample.ID)), 
                                   HUGO.symbol = rbind(cbind(as.character(mamma.final$HUGO.symbol)), 
                                                       cbind(as.character(cnvCombined$gene))),
                                   AF_CN = rbind(cbind(mamma.final$AF),                           
                                                 cbind(cnvCombined$copy.number)),
                                   combined = rbind(cbind(mamma.final$combined), cbind(cnvCombined$alteration)),stringsAsFactors = F)
  
  # OVERLAP BETWEEN FOUND ALTERATIONS AND ONCOKB DB
  overlapOnco = combinedAlterations[complete.cases(match(combinedAlterations$combined, oncokb$combined)),]
  overlapOnco = merge(overlapOnco,oncokb,by="combined")
  
  # UPDATE DATA FRAME
  clinicallySig.OncoKB = as.data.frame(cbind( as.character(overlapOnco$Sample.ID),
                                              as.character(overlapOnco$HUGO.symbol),
                                              as.character(overlapOnco$combined),
                                              as.character(overlapOnco$Level),
                                              as.character(overlapOnco$AF_CN)),stringsAsFactors = F)
  names(clinicallySig.OncoKB) = c("Sample.ID", "HUGO.symbol","Gene.Variation","Level","AF/CN")
  clinicallySig.OncoKB$Source = rep("OncoKB", length(clinicallySig.OncoKB$Sample.ID))
  clinicallySig.OncoKB$ID.Variant = interaction(clinicallySig.OncoKB$Sample.ID, clinicallySig.OncoKB$Gene.Variation)
  
  ### 2. FIND ALTERATIONS LISTED IN PMKB
  
  # CREATE DATA FRAME
  PMKBAll$Variants = gsub(x = PMKBAll$Variants, pattern = "copy.number.gain", replacement = "Gain")
  PMKBAll$Variants = gsub(x = PMKBAll$Variants, pattern = "copy.number.loss", replacement = "Loss")
  PMKBAll$combined = PMKBAll$Variants
  combinedAlterations = data.frame(Sample.ID = rbind(cbind(as.character(mamma.final$Sample.ID)),
                                                     cbind(cnvCombined$Sample.ID)), 
                                   HUGO.symbol = rbind(cbind(as.character(mamma.final$HUGO.symbol)), 
                                                       cbind(as.character(cnvCombined$gene))),
                                   AF_CN = rbind(cbind(mamma.final$AF),                           
                                                 cbind(cnvCombined$copy.number)),
                                   combined = rbind(cbind(mamma.final$combined), cbind(cnvCombined$alteration)),stringsAsFactors = F)
  
  # OVERLAP BETWEEN FOUND ALTERATIONS AND PMKB DB
  overlapPMKB = combinedAlterations[complete.cases(match(combinedAlterations$combined,PMKBAll$Variants)),]
  overlapPMKB = merge(overlapPMKB, PMKBAll,by="combined")
  
  # UPDATE DATAFRAME
  clinicallySig.PMKB = as.data.frame(cbind( as.character(overlapPMKB$Sample.ID),
                                            as.character(overlapPMKB$HUGO.symbol),
                                            as.character(overlapPMKB$combined),
                                            as.character(overlapPMKB$Tissue.Type.s.),
                                            as.character(overlapPMKB$Tier),
                                            as.character(overlapPMKB$AF_CN)),stringsAsFactors = F)
  names(clinicallySig.PMKB) = c("Sample.ID", "HUGO.symbol","Gene.Variation","Tissue","PMKB.Tier","AF/CN")
  clinicallySig.PMKB$ID.Variant = interaction(clinicallySig.PMKB$Sample.ID, clinicallySig.PMKB$Gene.Variation)
  clinicallySig.PMKB$Source = rep("PMKB", length(clinicallySig.PMKB$ID.Variant))
  
  # ADD TISSUE INFORMATION TO DATAFRAME
  aggTissues = aggregate(Tissue~ID.Variant,data = clinicallySig.PMKB,FUN= paste)$Tissue
  collapsedTissues = sapply(X = aggTissues,FUN = function(x) paste(x, collapse =  "", sep = ", "))
  
  # REMOVE DUPLICATES
  clinicallySig.PMKB = clinicallySig.PMKB[!duplicated(clinicallySig.PMKB$ID.Variant),]
  clinicallySig.PMKB$Tissue = collapsedTissues
  clinicallySig.PMKB$TissueSummary = rep("Breast", length(clinicallySig.PMKB$Sample.ID))
  clinicallySig.PMKB$TissueSummary[grep(pattern = "reast", x = collapsedTissues,invert = T)] = "Other"
  
  
  
  ### 3. COMBINE ALL CLINICALLY SIGNIFICANT ALTERATIONS - OVERLAP
  
  # FIND COMMON/OVERLAPPING
  clinicalSig.Merged = merge(clinicallySig.PMKB,clinicallySig.OncoKB,by="ID.Variant",all = T)
  detected.By.PMKB = complete.cases(clinicalSig.Merged$Source.x)
  detected.By.OncoKB = complete.cases(clinicalSig.Merged$Source.y)
  detected.By.Both = detected.By.PMKB == detected.By.OncoKB
  
  # ADD ANNOATIONS 
  clinicalSig.Merged$Detectedby = rep("", length(detected.By.PMKB))
  clinicalSig.Merged$Detectedby[which(detected.By.PMKB)] = "PMKB"
  clinicalSig.Merged$Detectedby[which(detected.By.OncoKB)] = "OncoKB"
  clinicalSig.Merged$Detectedby[which(detected.By.Both)] = "Both"
  
  # IDENTIFIER
  extract.gene.variant = function(x) {
    splitted = stringr::str_split(string = x, pattern = "[.]")[[1]]
    return (paste(splitted[2], paste(splitted[3:length(splitted)],collapse = "."), sep = "."))
  }
  clinicalSig.Merged$gene.variants.combined = sapply(clinicalSig.Merged$ID.Variant, FUN = extract.gene.variant)
  clinicalSig.Merged$SampleID.combined = sapply(clinicalSig.Merged$ID.Variant, FUN = function(x) stringr::str_split(x, pattern="[.]")[[1]][1])
  
  return(list(clinicalSig.Merged = clinicalSig.Merged, overlapOnco = overlapOnco))
}



