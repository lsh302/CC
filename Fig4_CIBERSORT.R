library(CIBERSORT)
library(reshape2)
library(ggpubr) 
library(dplyr)

sig_matrix <- system.file("extdata", "LM22.txt", package = "CIBERSORT")
data(LM22)
mixed_expr <- read.csv("expdataTumor.csv",check.names = F,row.names = 1)
results <- cibersort(sig_matrix = LM22, mixture_file = mixed_expr,perm = 1000,QN = F)
data=filter(results,P-value < 0.05)
data=as.matrix(data[,1:(ncol(data)-3)])
cluster=read.table("cluster.txt", header=T, check.names=F, row.names=1)
sameSample=intersect(row.names(data), row.names(cluster))
data=cbind(data[sameSample,,drop=F], cluster[sameSample,,drop=F])
data$cluster=factor(data$cluster, levels=c('young',"medium",'old')) 
data1=melt(data,id.vars=c("cluster"))
colnames(data1)=c("cluster", "celltype", "value")
bioCol=c("#00468B","#925E9F","#759EDD")
ggboxplot(data1, x="celltype", y="value", fill="cluster",
                  orientation="horizontal", 
                  outlier.shape=NA,
                  xlab="",
                  ylab="Fraction",
                  legend.title="cluster",
                  width=0.8,
                  size = 0.2,
                  palette=bioCol)+
  theme_test()+
  stat_compare_means(aes(group=cluster),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")),label="p.signif", method="kruskal.test",label.y= 0.45) 
