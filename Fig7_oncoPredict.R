library(oncoPredict)
library(ggplot2)
library(tidyverse)

rt <- read.csv("expdataTumor.csv",check.names = F,row.names = 1)
data=as.matrix(rt)
calcPhenotype(trainingExprData = Expr,    
              trainingPtype = Res,        
              testExprData = data,            
              batchCorrect = 'eb',              
              powerTransformPhenotype = F,      
              removeLowVaryingGenes = 0.2,      
              minNumSamples = 20,               
              printOutput = TRUE,               
              removeLowVaringGenesFrom = 'homogenizeData') 
senstivity <- read.csv("calcPhenotype_Output/DrugPredictions.csv", header=T, sep=",", check.names=F, row.names=1)
cl=read.table("cluster.txt", header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(cl), row.names(senstivity))
cl=cl[sameSample,,drop=F]
senstivity=senstivity[sameSample,,drop=F]
all(rownames(Type)==rownames(senstivity)) 
rt=cbind(cl, senstivity)
rt$cluster=factor(rt$cluster, levels= c('young',"medium",'old'))
rt2 <- data.frame(drug=colnames(senstivity),Pvalue=NA)
rownames(rt2) <- rt2$drug
obj <- rownames(rt2)
for(i in obj){
  result1 <- rt 
  result1 <- rename(result1, marker = i)
  test=kruskal.test(marker ~ cluster, data=result1) 
  diffPvalue=test$p.value
  rt2[i,2]=diffPvalue
  if(diffPvalue<0.05){ 
    p<- ggstatsplot::ggbetweenstats(
      data = result1, 
      x = cluster,
      y = marker,  
      fill = cluster,
      xlab = NULL,
      ylab = "Imputed Sensitivity Score",
      type = "np",     
      pairwise.display = "all",  
      p.adjust.method = "fdr",  
      title = i, 
      k = 3L, 
      centrality.plotting = F,
      point.args = list(position = ggplot2::position_jitterdodge(jitter.width = 0.3,dodge.width = 0.1),alpha =1, 
                        size = 2, stroke = 0, na.rm = TRUE),
      boxplot.args = list(aes(fill=cluster),notch=F,width = 0.25, size=0.4,alpha = 0.4, na.rm = TRUE),
      violin.args = list(aes(fill=cluster),width = 0.95, size=0.4,alpha = 0.2, na.rm = TRUE),
      ggsignif.args = list(textsize = 2, tip_length = 0.01,size=0.4, na.rm = TRUE),
      ggtheme = theme_classic())+
      theme(axis.line = element_line(linewidth = 0.4),
            axis.title.y.left  = element_text(size = 12),
            axis.title.y.right  = element_text(size = 5),
            plot.title = element_text(hjust = 0.5),
            axis.text = element_text(size = 12),plot.subtitle = element_text(size = 5))+ 
      ggplot2::scale_color_manual(values = c("#00468B","#925E9F","#759EDD"))+ 
      ggplot2::scale_fill_manual(values = c("#00468B","#925E9F","#759EDD")) 
    ggsave(paste0("chemosensitivity",i,".pdf"),p,height = 4,width = 4)
  }
}