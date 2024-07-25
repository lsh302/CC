library(maftools)
library(R.utils)
library(tidyverse)

cl <- read.csv("clinicalData.csv", check.names=F)
cl <- rename(cl,Tumor_Sample_Barcode=bcr_patient_barcode)
mafFilePath <- dir(path = "CESC",pattern = "masked.maf.gz$",full.names = T,recursive=T)
mafdata <- lapply(mafFilePath,  
                  function(x){
                    read.maf(x,
                             clinicalData = cl, 
                             isTCGA=TRUE)
                  })
maf <- merge_mafs(mafdata)
maf
cl_young <- subset(cl, cluster%in%c("young"))$Tumor_Sample_Barcode
cl_medium <- subset(cl, cluster%in%c("medium"))$Tumor_Sample_Barcode
cl_old <- subset(cl, cluster%in%c("old"))$Tumor_Sample_Barcode
maf_young <- subsetMaf(maf = maf, tsb = cl_young, isTCGA = TRUE)
maf_medium <- subsetMaf(maf = maf, tsb = cl_medium, isTCGA = TRUE)
maf_old <- subsetMaf(maf = maf, tsb = cl_old, isTCGA = TRUE)
plotmafSummary(maf = maf,
                     rmOutlier = TRUE,
                     addStat = "median",
                     dashboard = TRUE,
                     titvRaw = FALSE)
titv <- titv(maf = maf,
             plot = F,
             useSyn = TRUE)
plotTiTv(res = titv)
somaticInteractions(maf = maf,
                    top = 25,
                    pvalue = c(0.05, 0.1))
vc_cols = c("#1F6CB6","red3","#70BC6B","#F4B9C5","#784194","#B81A7B","#A65628","#9E1F63") 
names(vc_cols) = c('Missense_Mutation','Multi_Hit','Frame_Shift_Del','Nonsense_Mutation','Frame_Shift_Ins','In_Frame_Ins',
  'Splice_Site','In_Frame_Del')
col1 = c("#00468B","#925E9F","#759EDD")
names(col1) = c('young','medium','old')
oncoplot(maf = maf,
         top = 25, 
         colors = vc_cols, 
         clinicalFeatures = 'cluster',
         sortByAnnotation = T, 
         annotationColor = list(cluster=col1),
         anno_height = 0.5, 
         legend_height = 4, 
         drawRowBar = T, 
         drawColBar = T, 
         draw_titv = F) 
fab.ce <- clinicalEnrichment(maf = maf,clinicalFeature = 'cluster')
plotEnrichmentResults(enrich_res = fab.ce,
                      pVal = 0.05,
                      geneFontSize = 0.5,
                      annoFontSize = 0.35,
                      legendFontSize=0.6)
mafSurvival(maf = maf,
            genes = 'SMC4', 
            time = 'OS.time',
            Status = 'OS', 
            isTCGA = TRUE)