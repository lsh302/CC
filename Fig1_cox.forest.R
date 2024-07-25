library(ggplot2)
library(survival)
library(survminer)
library(tidyverse)

meta<-read.csv("clinicalData.csv",header = T,row.names = 1) 
meta$OS.time <- meta$OS.time/30
meta$DSS.time <- meta$DSS.time/30 
meta$DFI.time <- meta$DFI.time/30 
meta$PFI.time <- meta$PFI.time/30
surType <- c("OS","DSS","PFI","DFI")
for(i in surType){
  j <- paste0(i,".time")
  meta1 <- meta
  meta1 <- rename(meta1, survival.time = j, survival = i) 
  if (T) {
    res.cox <- coxph(Surv(survival.time, survival) ~ Age + Race + Pathology + Grade + Stage,data = meta1) 
    summary(res.cox)
  }
  p1 <- ggforest(res.cox,  
                 data = meta1, 
                 main = i,  
                 cpositions = c(0.01, 0.15, 0.35,0.55), 
                 fontsize = 1,
                 refLabel = 'reference', 
                 noDigits = 4) 
  ggsave(paste0(i,".pdf"),p1,width = 10, height = 6)
}