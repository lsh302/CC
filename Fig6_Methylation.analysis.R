library(TCGAbiolinks)
library(tidyverse)
library(ChAMP)
require(GEOquery) 

cesc_methy <- GDCquery(
  project = "TCGA-CESC", 
  data.category = "DNA Methylation", 
  data.type = "Methylation Beta Value",
  platform = "Illumina Human Methylation 450")
GDCdownload(cesc_methy)
txts <- dir(
  path = 'TCGA-CESC-Methylation/',
  pattern = "*.methylation_array.sesame.level3betas.txt$", 
  recursive = TRUE)
beta_value <- function(x){
  dff <- data.table::fread(file.path('TCGA-CESC-Methylation/',x),select = 2)
  return(dff)
}
value_list <- lapply(X = txts,
                     FUN = beta_value)
beta <- do.call(cbind, value_list)
beta <- as.data.frame(beta)
rownames(beta) <- test1$V1
sample_sheet <- data.table::fread('gdc_sample_sheet.tsv',header = T)
ids <- data.frame(file_name = sample_sheet$`File Name`,sample_id = sample_sheet$`Sample ID`)
txts2 <- str_split(txts, "/", simplify = T) 
txts3 <- txts2[,2]
ids2 <- ids[match(txts3, ids$file_name),]
colnames(beta) <- ids2$sample_id
tumor_col <-  colnames(beta)[substr(colnames(beta),14,15) =="01"] 
beta_tumor <- beta[,tumor_col]
colnames(beta_tumor) <- substr(colnames(beta_tumor), 1, 12)
beta_tumor <- as.matrix(beta_tumor)
beta_tumor <- impute.knn(beta_tumor) 
betaData <- beta_tumor$data
betaData <- betaData+0.0001
beta_tumor <- betaData
cl<- read.table("cluster.txt",header=T)
rownames(cl) <- cl[,1]
cl <- cl[colnames(beta_tumor),,drop=F]  
myLoad=champ.filter(beta = beta_tumor,pd = cl) 
myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=5) 
pD <- myLoad$pd
group_list=pD$cluster
myDMP <- champ.DMP(beta = myNorm,pheno=group_list) 