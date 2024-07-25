library(Seurat) 
library(tidyverse) 
library(patchwork) 
library(ggplot2)
library(pheatmap)
library(reshape2)
library(RColorBrewer)
library(scRNAtoolVis)
library(GSVA)
library(org.Hs.eg.db)
library(msigdbr)

Cells.sub<- subset(sce@meta.data,celltype=="Epithelial") 
scRNAsub <- subset(sce,cells=row.names(Cells.sub))
scRNAsub <- FindVariableFeatures(scRNAsub,selection.method = "vst",nfeatures=2000) 
scale.genes <- rownames(scRNAsub)
scRNAsub <- ScaleData(scRNAsub,features = scale.genes)
scRNAsub <- RunPCA(scRNAsub,features = VariableFeatures(scRNAsub))
ElbowPlot(scRNAsub,ndims=20,reduction="pca")
pc.num=1:10
scRNAsub <- FindNeighbors(scRNAsub,dims = pc.num) 
scRNAsub <- FindClusters(scRNAsub,resolution =0.15) 
cell_cluster <- data.frame(cell_ID=rownames(metadata),cluster_ID=metadata$seurat_clusters) 
scRNAsub <- RunUMAP(scRNAsub,dims = pc.num) 
DimPlot(scRNAsub,reduction = "umap",label = T) 
markers <- FindAllMarkers(object = scRNAsub,test.use="wilcox",
                          only.pos = F, 
                          logfc.threshold =0.25) 
all.markers =markers %>% dplyr::select(gene,everything()) %>% subset(p_val<0.05) 
jjVolcano(diffData = all.markers,
          log2FC.cutoff = 0.25,
          size =3.5,
          fontface = 'italic',
          topGeneN =3
)
scRNAsub@meta.data$group <- ifelse(scRNAsub@meta.data$seurat_clusters %in% c("0","1","2","3","4","6","7","8","10","11","12","13","14"),"rest",ifelse(scRNAsub@meta.data$seurat_clusters =="5","5","9"))
scRNAsub@meta.data$group <- factor(scRNAsub@meta.data$group,levels=c("5","9","rest"))
Idents(scRNAsub) <- scRNAsub@meta.data$group
expr <- AverageExpression(scRNAsub, assays = "RNA", layer = "data")[[1]]
expr <- expr[rowSums(expr)>0,]  
expr <- as.matrix(expr)
human_KEGG = msigdbr(species = "Homo sapiens",
                     category = "H")
human_KEGG_Set = human_KEGG %>% split(x = .$gene_symbol, f = .$gs_name)
gsva.kegg <- gsva(expr, gset.idx.list = human_KEGG_Set, 
                  kcdf="Gaussian", 
                  method = "gsva",
                  parallel.sz=1)
df1 <- as.data.frame(t(scale(t(gsva.kegg))))
df1 <- df1[order(df1[,1], decreasing = TRUE), ] 
colnames(df1) <- c("Epi_5","Epi_9","Epi_rest")
pheatmap(df1,
         show_colnames = T, 
         angle_col = "0",
         fontsize_row = 7,
         fontsize_col = 10,
         cluster_row = F,
         cluster_col = F,
         color = colorRampPalette(c("#4B5B95", "white", "#CD322C"))(50))
filelist<-paste0("CancerSEA/",dir("CancerSEA/")) 
CancerSEA_Set<-list()
for(i in 1:length(filelist)){   
  CancerSEA_Set[[i]]<-read.table(filelist[i],header = T,sep = "\t")[,2]
}
names(CancerSEA_Set) <-  str_extract(filelist, "(?<=/).*?(?=\\.)")
gsva.cancersea <- gsva(expr, gset.idx.list = CancerSEA_Set, 
                  kcdf="Gaussian", 
                  method = "gsva",
                  parallel.sz=1)
gsva.cancersea1 <- as.data.frame(t(scale(t(gsva.cancersea))))
df2 <- gsva.cancersea1
colnames(df2) <- paste("Epi", c(0:14), sep = "_")
pheatmap(df2, 
         show_colnames = T,
         angle_col = "45",
         fontsize_row = 8.5, 
         fontsize_col = 8.5,
         cluster_row = F,
         cluster_col = F,
         color = colorRampPalette(c("#4B5B95", "white", "#CD322C"))(50))
dff <- as.data.frame(table(scRNAsub@meta.data$seurat_clusters,scRNAsub@meta.data$Phase))
dff$cluster <- paste0("Epi_",dff$cluster)
dff$cluster <- factor(dff$cluster,levels=paste0("Epi_",c(0:14)))
dff$stage <- factor(dff$stage, levels = c("G2M","S","G1"))
ggplot(data = dff,aes(x = cluster,y = cells, fill = stage))+
  geom_bar(stat = 'identity',position = 'fill',width = 0.9)+ 
  labs(x = "",y = "Relative abundance")+
  scale_fill_manual(values = color)+
  guides(fill = guide_legend(ncol = 1,bycol = T,override.aes = list(size = 5)))+
  theme_test()