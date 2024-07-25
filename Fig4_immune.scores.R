library(ggplot2)
library(tidyverse)

exp=read.csv("expdataTumor.csv" , header=T, check.names=F,row.names = 1)
target_genelist <- c("LAG3","ITGA5","ESM1","DES","CXCL2") 
target_exp <- exp[target_genelist,]
for (i in c(1:304)){
  target_exp[j+1,i] =  -0.297*target_exp[1,i] + 0.334*target_exp[2,i]+0.19*target_exp[3,i]-0.214*target_exp[4,i] + 0.115*target_exp[5,i]
}
rownames(target_exp)[j+1] <- target
target_exp <- as.data.frame(t(target_exp))
cluster=read.table("cluster.txt", header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(target_exp), row.names(cluster))
rt1=cbind(target_exp[sameSample,,drop=F], cluster[sameSample,,drop=F])
result1 <- rt1[,c(j+1,j+2)]
result1$cluster <- factor(result1$cluster,levels =c("young","medium","old"))
colnames(result1) <- c("score","cluster")
test=kruskal.test(score ~ cluster, data=result1) 
p <- ggstatsplot::ggbetweenstats(
      data = result1, 
      x = cluster,
      y = score,  
      type = "np",     
      pairwise.display = "all",  
      p.adjust.method = "fdr",  
      title = target,
      xlab = "", 
      ylab = "Scores", 
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