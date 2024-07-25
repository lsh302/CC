library(Seurat) 
library(harmony) 
library(tidyverse) 
library(ggplot2)
library(patchwork) 

dir_name=list.files("GSE197461_GSE208653/")
scRNAlist <- list()
for(i in 1:length(dir_name)){
  counts <- Read10X(data.dir = paste("GSE197461_GSE208653/",dir_name[i],sep =""))
  scRNAlist[[i]] <- CreateSeuratObject(counts,project = dir_name[i],min.cells = 3, min.features=200)
}
for(i in 1:length(scRNAlist)){
  sc <- scRNAlist[[i]]
  sc[["mt_percent"]] <- PercentageFeatureSet(sc,pattern ="^MT-") 
  HB_genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ") 
  HB_m <- match(HB_genes,rownames(sc@assays$RNA)) 
  HB_genes <- rownames(sc@assays$RNA)[HB_m]
  HB_genes <- HB_genes[!is.na(HB_genes)]
  sc[["HB_percent"]] <- PercentageFeatureSet(sc,features=HB_genes) 
  scRNAlist[[i]] <- sc
  rm(sc)
}
violin_before <-list()
for(i in 1:length(scRNAlist)){
  violin_before[[i]] <- VlnPlot(scRNAlist[[i]],
                                features = c("nFeature_RNA","nCount_RNA","mt_percent","HB_percent"),
                                pt.size =0.01, 
                                ncol =4)
}
violin_before_filter <- CombinePlots(plots = violin_before, nrow=length(scRNAlist),legend='none')
scRNAlist <- lapply(X = scRNAlist,FUN = function(x){ 
  x<- subset(x,
            subset = nFeature_RNA > 200 & nFeature_RNA <5000 &
            mt_percent < 20 & 
            HB_percent < 3 &
            nCount_RNA < quantile(nCount_RNA,0.97) & 
            nCount_RNA > 1000)})
violin_after <-list()
for(i in 1:length(scRNAlist)){
  violin_after[[i]] <- VlnPlot(scRNAlist[[i]],
                                features = c("nFeature_RNA","nCount_RNA","mt_percent","HB_percent"),
                                pt.size =0.01, 
                                ncol =4)
}
violin_after_filter <- CombinePlots(plots = violin_after, nrow=length(scRNAlist),legend='none')
scRNAlist <- merge(x=scRNAlist[[1]],y=scRNAlist[-1])
violin_after_filter_combined <- VlnPlot(scRNAlist,
                        features = c("nFeature_RNA","nCount_RNA","mt_percent","HB_percent"), pt.size = 0.01, ncol=4)
scRNAlist <- NormalizeData(scRNAlist) %>% 
  FindVariableFeatures(selection.method = "vst",nfeatures =2000) %>% 
  ScaleData() %>% 
  RunPCA(npcs = 20,verbose = T) 
DimPlot(scRNAlist,reduction = "pca",group.by ="orig.ident") 
g2m_genes = cc.genes$g2m.genes
g2m_genes = CaseMatch(search = g2m_genes,match = rownames(scRNAlist))
s_genes = cc.genes$s.genes
s_genes = CaseMatch(search = s_genes,match =rownames(scRNAlist)) 
scRNAlist <- CellCycleScoring(object=scRNAlist,g2m.features=g2m_genes,s.features=s_genes) 
scRNAlist=CellCycleScoring(object =scRNAlist,
                           s.features = s_genes,
                           g2m.features = g2m_genes,
                           set.ident =TRUE) 
scRNAlist <- CellCycleScoring(object=scRNAlist,g2m.features=g2m_genes,s.features=s_genes) 
VlnPlot(scRNAlist,features = c("S.Score","G2M.Score"), group.by = "orig.ident",
                                  ncol = 2, pt.size =0.1)
cellcy <- scRNAlist@meta.data %>% ggplot(aes(S.Score,G2M.Score))+geom_point(aes(color=Phase))+ theme_minimal()
scRNA_harmony <- RunHarmony(scRNAlist,group.by.vars = "orig.ident") 
DimPlot(scRNA_harmony,reduction = "harmony",group.by = "orig.ident")
ElbowPlot(scRNA_harmony,ndims=50,reduction="harmony")
scRNA_harmony <- FindNeighbors(scRNA_harmony,reduction ="harmony", dims = 1:10) %>% FindClusters(resolution = 0.5)
scRNA_harmony <- RunUMAP(scRNA_harmony,reduction = "harmony",dims = 1:10)
DimPlot(scRNA_harmony,reduction ="umap",group.by = "orig.ident")
DimPlot(scRNA_harmony,reduction = "umap",label =TRUE) 
markers <- FindAllMarkers(object =scRNA_harmony,test.use="wilcox",
                          only.pos =TRUE,
                          logfc.threshold =0.25)
all.markers =markers %>% dplyr::select(gene,everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10,wt = avg_log2FC)
#Epithelial
genes_to_check = c("EPCAM","KRT18","CDKN2A","CDH1","KLF5")
DotPlot(scRNA_harmony,features = genes_to_check, assay="RNA",group.by = "seurat_clusters") + coord_flip() +ggtitle("Epithelial")
#Endothelial                  
genes_to_check = c("PECAM1","VWF","EMCN","ENG", "CDH5")
DotPlot(scRNA_harmony,features = genes_to_check, assay="RNA",group.by = "seurat_clusters") + coord_flip() +ggtitle("Endothelial")
#Fibroblast
genes_to_check = c("DCN","COL1A1","COL3A1","ACTA2","TAGLN")
DotPlot(scRNA_harmony,features = genes_to_check, assay="RNA",group.by = "seurat_clusters") + coord_flip() +ggtitle("Fibroblast")
#Smooth muscle
genes_to_check = c("TAGLN","ACTA2","RGS5")
DotPlot(scRNA_harmony,features = genes_to_check, assay="RNA",group.by = "seurat_clusters") + coord_flip() +ggtitle("Smooth muscle")
#T/NK
genes_to_check = c("CD3D","NKG7","CD8A","CCL5","GZMA","CD3E","CD3G")
DotPlot(scRNA_harmony,features = genes_to_check, assay="RNA",group.by = "seurat_clusters") + coord_flip() +ggtitle("NK/T")
#B cell
genes_to_check = c("CD19","CD79A","MS4A1","BANK1")
DotPlot(scRNA_harmony,features = genes_to_check, assay="RNA",group.by = "seurat_clusters") + coord_flip() +ggtitle("B cell")
#Plasma
genes_to_check = c("IGHG1","TNFRSF17","MZB1","IGKC","IGHG3","XBP1","JCHAIN")
DotPlot(scRNA_harmony,features = genes_to_check, assay="RNA",group.by = "seurat_clusters") + coord_flip() +ggtitle("Plasma")
#Myeloid
genes_to_check = c("CD14","C1QA","CD163","CD68","CSF1R","LYZ")
DotPlot(scRNA_harmony,features = genes_to_check, assay="RNA",group.by = "seurat_clusters") + coord_flip() +ggtitle("Myeloid")
#Neutrophils
genes_to_check = c("CSF3R","CPA3","FCGR3B","NCF1","SORL1")
DotPlot(scRNA_harmony,features = genes_to_check, assay="RNA",group.by = "seurat_clusters") + coord_flip() +ggtitle("Neutrophils")
#Mast
genes_to_check = c("MS4A2","KIT","TPSAB1","CPA3")
DotPlot(scRNA_harmony,features = genes_to_check, assay="RNA",group.by = "seurat_clusters") + coord_flip() +ggtitle("Mast")
celltype=data.frame(ClusterID=0:18,celltype='unkown')
celltype[celltype$ClusterID %in% c(3,7,10,13),2]='Epithelial cell' 
celltype[celltype$ClusterID %in% c(14),2]='Endothelial cell' 
celltype[celltype$ClusterID %in% c(2),2]='Fibroblast'
celltype[celltype$ClusterID %in% c(12),2]='Smooth muscle cell' 
celltype[celltype$ClusterID %in% c(0,1,11,17),2]='T/NK cell' 
celltype[celltype$ClusterID %in% c(5,8),2]='B cell' 
celltype[celltype$ClusterID %in% c(5),2]='Plasma cell' 
celltype[celltype$ClusterID %in% c(4,9,16,18),2]='Myeloid cell'
celltype[celltype$ClusterID %in% c(6),2]='Neutrophil' 
celltype[celltype$ClusterID %in% c(15),2]='Mast cell' 
sce=scRNA_harmony
sce@meta.data$celltype = NA 
for(i in 1:nrow(celltype)){
  sce@meta.data[which(sce@meta.data$seurat_clusters == celltype$ClusterID[i]),"celltype"] <- celltype$celltype[i]} 
genes_to_check = c("EPCAM","KRT18","CDH1","KLF5",
                   "PECAM1","VWF","EMCN","ENG", "CDH5",
                   "DCN","COL1A1","COL3A1",
                   "TAGLN","ACTA2","RGS5",
                   "CD3D","NKG7","CD8A","CCL5","GZMA","CD3E","CD3G",
                   "CD19","MS4A1","BANK1",
                   "IGHG1","TNFRSF17","MZB1","IGKC","IGHG3","JCHAIN",
                   "C1QA","CD68","LYZ",
                   "CSF3R","FCGR3B","NCF1","SORL1",
                   "MS4A2","KIT","TPSAB1","CPA3")
DotPlot(sce,features = genes_to_check, assay="RNA",group.by = "celltype") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+coord_flip()
DimPlot(sce,reduction = "umap",group.by = "celltype",label =T, pt.size=1)