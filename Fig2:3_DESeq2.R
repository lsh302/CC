library(DESeq2)
library(ggplot2)
library(tidyverse)

data<-read.csv("expdataTumor_count.csv",header = T,check.names = F,row.names = 1) 
target <- "Young" 
Non_target <- "Old" 
list <- read.table("cluster.txt",header = T,check.names = F,row.names = 1)
list.group_target <- filter(list,cluster=="young") 
group_target <- intersect(rownames(list.group_target),colnames(data))
target_number <- length(group_target) 
list.group_Nontarget <- filter(list,cluster =="old") 
group_Nontarget <- intersect(rownames(list.group_Nontarget),colnames(data))
Nontarget_number <- length(group_Nontarget)
data.target <- data[,group_target] 
data.Nontarget<- data[,group_Nontarget]
all(row.names(data.target)==row.names(data.Nontarget))
data <- cbind(data.target,data.Nontarget)
condition <- as.factor(c(rep(target,target_number),rep(Non_target,Nontarget_number)))
condition <- relevel(condition,ref=Non_target) 
condition
colData <- data.frame(condition=condition,row.names = colnames(data))
dds <- DESeqDataSetFromMatrix(countData = data, colData = colData, design = ~ condition) 
vst <- vst(dds)
plotPCA(vst)
dds <- DESeq(dds)
res <- results(dds,contrast = c("condition",target,Non_target)) 
res
res <- as.data.frame(res)
res$name <- rownames(res)
cpm <- fpm(dds,robust = F)
cpm <- as.data.frame(cpm)
cpm$name <- rownames(cpm)
all(rownames(cpm)==rownames(res))
res.cpm <- merge(cpm,res,by="name")
res.cpm$change <- as.factor(ifelse(res.cpm$padj<0.05 & abs(res.cpm$log2FoldChange)>log2(2),  
                                       ifelse(res.cpm$log2FoldChange>log2(2), "Up", "Down"), "NoDiff"))
volcano.data <- res.cpm
rownames(volcano.data) <- volcano.data$name 
volcano.data <- volcano.data[complete.cases(volcano.data$padj),] 
suppressMessages(library(ggrepel))
up <- volcano.data[volcano.data$log2FoldChange>1 & volcano.data$padj<0.05,]
down <- volcano.data[volcano.data$log2FoldChange< -1 & volcano.data$padj<0.05,]
up.top <- rownames(up[order(up$padj,-up$log2FoldChange),])[1:5]
down.top <- rownames(down[order(down$padj,down$log2FoldChange),])[1:5]
top_genes <- c(up.top,down.top)
top_genes_data <- volcano.data[top_genes,]
valcano <- ggplot(data=volcano.data, aes(x=log2FoldChange, y=-log10(padj), color=change)) + 
  geom_point(alpha=0.8, size=1.5) + 
  theme_bw(base_size=15) + 
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank(),
    axis.title=element_text(size=15),
    axis.text=element_text(size=15)) +
  ggtitle(paste0(target,"VS",Non_target)) + 
  scale_color_manual(name="", values=c("red", "blue", "black"), limits=c("Up", "Down", "NoDiff")) + 
  geom_vline(xintercept=c(-log2(2), log2(2)), lty=2, col="gray", lwd=0.5) + 
  geom_hline(yintercept=-log10(0.05), lty=2, col="gray", lwd=0.5)+
  geom_text_repel(data = top_genes_data, aes(label = name), size = 2.5,color = "black")
deg.data <- res.cpm[res.cpm$change %in% c("Up", "Down"), ] 
rownames(deg.data) <- deg.data[,1] 
deg.data <- arrange(deg.data,desc(log2FoldChange)) 
deg.data <- deg.data[,2:(target_number + Nontarget_number + 1)] 
deg.data1=log2(deg.data+1) 
deg.data1 <- as.matrix(deg.data1) 
deg.data2 <- t(scale(t(deg.data1)))
group=data.frame(type = c(rep(target,target_number),rep(Non_target,Nontarget_number))) 
row.names(group)=colnames(deg.data2)
group_colors=list(type=c(Young = "#8c61ff", Old = "#36c3fe")) 
pheatmap::pheatmap(deg.data2,
                   cluster_cols = F,
                   cluster_rows = F,
                   angle_col=315, display_numbers = FALSE,
                   fontsize_number = 10,
                   legend = T, show_rownames = F, show_colnames = F,
                   annotation_col = group,
                   annotation_colors = group_colors,
                   filename = paste0(target,"VS",Non_target,"_heatmap.pdf")) 