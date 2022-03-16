[Github.io page](https://osamashiraz.github.io/LO_Lab_CASR_2021/)

main.RMD contains the code used to generate data presented in following
work: (Identifying Genomic Alterations in Stage IV Breast Cancer
Patients using MammaSeq<sup>TM</sup>: An International Collaborative
Study)\[<a href="https://www.sciencedirect.com/science/article/abs/pii/S1526820920302135" class="uri">https://www.sciencedirect.com/science/article/abs/pii/S1526820920302135</a>\]

# PREPARATION
```r

# Install all these packages
#>.....

setwd("../scripts/")

# PREPARE ENVIRONMENT
source("helperFunctions.R")

# LOAD INPUT RNASEQ DATA AND ANNOTATIONS
source("./preparation/a.prepare_RNAseq_exp_mats.R")
source("./preparation/b.prepare_annotations.R")
clinData = prepareAnnotations()

source("./preparation/load_colors.R")
#source("./preparation/c.prepare_DNAseq_mat.R")

# Supplementary data 2 
icesTAF::mkdir("../output/supplementary_data")
writexl::write_xlsx(subset(clinData, Study_Code %in% "Ovary" & TissueType %notin% "Normal"), 
                    path = "../output/supplementary_data/Supplementary data 2 - PAM50 and tumor purity.xlsx")

```
