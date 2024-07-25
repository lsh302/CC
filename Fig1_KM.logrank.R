library(survival)
library(survminer)
library(RColorBrewer)
library(tibble) 
library(ggpp) 
library(tidyverse)

dat<-read.csv("clinicalData.csv",header = T)
dat$OS.time <- dat$OS.time/30 
dat$DSS.time <- dat$DSS.time/30 
dat$PFI.time <- dat$PFI.time/30 
dat$DFI.time <- dat$DFI.time/30 
surType <- c("OS","DSS","PFI","DFI")
for(i in surType){
  j <- paste0(i,".time")
  dat1 <- dat
  dat1 <- rename(dat1, survival.time = j, survival = i)
  fitd <- survdiff(Surv(survival.time, survival) ~ group,
                   data      = dat1,
                   na.action = na.exclude)
  p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  fit <- survfit(Surv(survival.time, survival)~ group,
                 data      = dat1,
                 type      = "kaplan-meier",
                 error     = "greenwood",
                 conf.type = "plain",
                 na.action = na.exclude)
  ps <- pairwise_survdiff(Surv(survival.time, survival)~ group,
                          data            = dat1,
                          p.adjust.method = "fdr") 
  mycol <- c("#00468B","#925E9F","#759EDD")
  names(fit$strata) <- gsub("group=", "", names(fit$strata)) 
  p <- ggsurvplot(fit               = fit,
                  conf.int          = FALSE, 
                  risk.table        = TRUE, 
                  risk.table.col    = "strata",
                  palette           = mycol, 
                  data              = dat1,
                  xlim              = c(0,120), 
                  size              = 1,
                  break.time.by     = 12, 
                  legend.title      = "",
                  xlab              = "Time (months)",
                  ylab              = i,
                  risk.table.y.text = FALSE,
                  tables.height     = 0.3) 
  p.lab <- paste0("log-rank test P",
                  ifelse(p.val < 0.001, " < 0.001", paste0(" = ",round(p.val, 3)))) 
  p$plot <- p$plot + annotate("text",
                              x = 0, y = 0.35, 
                              hjust = 0,
                              fontface = 4,
                              label = p.lab)
  addTab <- as.data.frame(as.matrix(ifelse(round(ps$p.value, 3) < 0.001, "<0.001",round(ps$p.value, 3))))
  addTab[is.na(addTab)] <- "-"
  df <- tibble(x = 0, y = 0, tb = list(addTab))
  p$plot <- p$plot + 
    geom_table(data = df, 
               aes(x = x, y = y, label = tb),
               table.rownames = TRUE)
  ggsave(paste0(i,".pdf"),p$plot,width = 6, height = 6)
}