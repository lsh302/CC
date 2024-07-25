library(CellChat)
library(patchwork)
library(Seurat)
library(ggplot2)

data.input = sce@assays$RNA@data
meta.data =  sce@meta.data
data.input = data.input[,row.names(meta.data)]
meta.data$celltype = factor(meta.data$celltype,levels = c("Epithelial", "Endothelial", "Fibroblast", "Smooth muscle cell",
                                                          "NK/T","B cell", "Plasma", "Myeloid", "Neutrophils","Mast"))
meta.data1 <- meta.data
meta.data1$celltype1 <- meta.data1$celltype
meta.data1_epi <- filter(meta.data1,celltype1=="Epithelial")
meta.data1_Non_epi<- filter(meta.data1,celltype1!="Epithelial") 
epi_subgroup <- read.csv("epithelial_cellTOcluster.csv",header = T,row.names = 1)
epi_subgroup$cluster_ID = factor(epi_subgroup$cluster_ID,levels = c("Epi_5", "Epi_9", "Epi_rest"))
epi_subgroup <- epi_subgroup[row.names(meta.data1_epi),,drop=F]
meta.data1_epi$celltype1 <- epi_subgroup$cluster_ID 
meta.data1 <- rbind(meta.data1_epi,meta.data1_Non_epi)
meta.data1$celltype1 = factor(meta.data1$celltype1,levels = c("Epi_5", "Epi_9", "Epi_rest",
                                       "Endothelial", "Fibroblast", "Smooth muscle cell", "NK/T", 
                                       "B cell", "Plasma", "Myeloid", "Neutrophils","Mast"))
data.input = data.input[,row.names(meta.data1)] 
cellchat <- createCellChat(object = data.input, 
                           meta = meta.data1, 
                           group.by = "celltype1")
cellchat <- setIdent(cellchat, ident.use = "celltype1") 
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB 
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) 
future::plan("multiprocess", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat) 
cellchat <- identifyOverExpressedInteractions(cellchat) 
cellchat <- computeCommunProb(cellchat,population.size = F) 
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
df.pathway = subsetCommunication(cellchat,slot.name = "netP")
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F,
                 title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F,
                 title.name = "Interaction weights/strength") 
netVisual_heatmap(cellchat,color.heatmap = c("white", "#b2182b"))
netVisual_heatmap(cellchat, measure = "weight",color.heatmap = c("white", "#b2182b"))
netAnalysis_contribution(cellchat, signaling = pathways.show)
pathways.show.all <- cellchat@netP$pathways
netVisual_bubble(cellchat, sources.use = c(1:3), targets.use = c(7:12), remove.isolate = FALSE,sort.by.target = T) 
netVisual_bubble(cellchat, sources.use = c(7:12), targets.use = c(1:3), remove.isolate = FALSE,sort.by.target = T) 
pathways.show <- c("MIF") 
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle",label.edge= T)