library(Seurat)
library(monocle)
library(CytoTRACE)

exp <- scRNAsub@assays[["RNA"]]@counts
fdata <- data.frame(gene_short_name = row.names(scRNAsub), row.names = row.names(scRNAsub))
pdata <- scRNAsub@meta.data
fd <- new("AnnotatedDataFrame", data = fdata) 
pd <- new("AnnotatedDataFrame", data = pdata)
CDS <- newCellDataSet(cellData = exp,phenoData = pd,featureData = fd)
CDS <- estimateSizeFactors(CDS)
CDS <- estimateDispersions(CDS)
CDS <- detectGenes(CDS, min_expr = 0.1)
expressed_genes <-  row.names(subset(fData(CDS),num_cells_expressed >= 20));length(expressed_genes)
clustering_DEG_genes <-differentialGeneTest(CDS[expressed_genes,],fullModelFormulaStr = '~seurat_clusters')
ordering_genes <-row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:2000]
CDS <- setOrderingFilter(CDS,ordering_genes = ordering_genes)
plot_ordering_genes(CDS)
CDS <- reduceDimension(CDS, method = 'DDRTree')
CDS <-orderCells(CDS)
p1 <- plot_cell_trajectory(CDS, color_by = "Pseudotime")
p2 <- plot_cell_trajectory(CDS, color_by = "State")
p3 <- plot_cell_trajectory(CDS, color_by = "seurat_clusters")
monocle_meta <- data.frame(t(CDS@reducedDimS), 
                           CDS$Pseudotime, 
                           CDS$State, 
                           CDS$seurat_clusters)
colnames(monocle_meta) <- c("C1", "C2", "Pseudotime", "State", "group")
phenot1 <- monocle_meta$group
phenot1 <- as.character(phenot1)
names(phenot1) <- rownames(monocle_meta)
emb_monocle <- monocle_meta[,1:2]
mat<-as.matrix(scRNAsub@assays$RNA@counts)
results <- CytoTRACE(mat = mat)
plotCytoTRACE(results, phenotype = phenot1, emb = emb_monocle)