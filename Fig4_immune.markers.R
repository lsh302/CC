library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(tidyverse)

data=read.csv("expdataTumor.csv", header=T, check.names=F,row.names = 1)
data=log2(data+1)
gene=read.table("markerlist", header=F, sep="\t", check.names=F)
sameGene=intersect(row.names(data),as.vector(gene[,1]))
data=t(data[sameGene,])
cluster=read.table("cluster.txt", header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(data), row.names(cluster))
rt1=cbind(data[sameSample,,drop=F], cluster[sameSample,,drop=F])
data=melt(rt1,id.vars=c("cluster")) 
colnames(data)=c("cluster", "Gene", "Expression")
data$cluster=factor(data$cluster, levels=c('young',"medium",'old')) 
boxplot=ggboxplot(data, x="Gene", y="Expression", fill="cluster",
				  orientation="horizontal",
				  xlab="",
				  ylab="Gene expression",
				  legend.title="Cluster",
				  width=0.8,
				  outlier.shape=NA, 
				  size = 0.2)+
  stat_compare_means(aes(group=cluster),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")),label="p.signif", method="kruskal.test",label.y= 12)+ 
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1), 
        axis.line = element_line(linewidth = 0),
        axis.title.x = element_text(face = "plain", size = 10),
        axis.text.x = element_text(face = "plain", size = 8),
        axis.text.y = element_text(face = "plain", size = 8))