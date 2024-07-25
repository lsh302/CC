library(clusterProfiler)
library(ggplot2)
library(org.Hs.eg.db)

target <- "Young"
Non_target <- "Old"
data <- read.csv("cpmResultDiff.csv")
data_sgni <- data[c("name","log2FoldChange","pvalue","padj")]
entrezid_all = mapIds(x = org.Hs.eg.db,  
                      keys = data_sgni$name, 
                      keytype = "SYMBOL", 
                      column = "ENTREZID") 
entrezid_all = na.omit(entrezid_all)  
entrezid_all = data.frame(entrezid_all) 
go_enrich = enrichGO(gene = entrezid_all[,1],
                     OrgDb = org.Hs.eg.db, 
                     keyType = "ENTREZID", 
                     ont = "ALL", 
                     pAdjustMethod = "fdr", 
                     pvalueCutoff = 1, 
                     qvalueCutoff = 1, 
                     readable = T) 
go_enrich  = data.frame(go_enrich) 
KEGG_enrich = enrichKEGG(gene = entrezid_all[,1], 
                         keyType = "kegg",
                         organism= "human",  
                         pAdjustMethod = "fdr", 
                         pvalueCutoff = 1,  
                         qvalueCutoff =1) 
KEGG_enrich  = data.frame(KEGG_enrich)
go_enrich.BP <- arrange(arrange(filter(go_enrich,ONTOLOGY=="BP"),pvalue)[1:6,]) 
go_enrich.CC <- arrange(filter(go_enrich,ONTOLOGY=="CC"),pvalue)[1:6,]
go_enrich.MF <- arrange(filter(go_enrich,ONTOLOGY=="MF"),pvalue)[1:6,]
go_enrich.BP <- arrange(go_enrich.BP,Count) 
go_enrich.MF <- arrange(go_enrich.MF,Count)
go_enrich <- rbind(go_enrich.BP,go_enrich.CC,go_enrich.MF)
go_enrich$Description <- factor(go_enrich$Description,levels=go_enrich$Description)
p1 <- ggplot(go_enrich,
       aes(x=Description,y=Count, fill=ONTOLOGY)) + 
  geom_bar(stat="identity", width=0.8) + 
  scale_fill_nejm(alpha = 0.7) +
  facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')+
  coord_flip() +  
  xlab("GO term") +  
  ylab("Gene_Number") +  
  labs(title = paste0(target," VS ",Non_target))+  
  theme_bw()
go_enrich.BP <- arrange(go_enrich.BP,GeneRatio) 
go_enrich.CC <- arrange(go_enrich.CC,GeneRatio)
go_enrich.MF <- arrange(go_enrich.MF,GeneRatio)
go_enrich <- rbind(go_enrich.BP,go_enrich.CC,go_enrich.MF)
go_enrich$Description <- factor(go_enrich$Description,levels=go_enrich$Description)
p2 <- ggplot(go_enrich,
       aes(y=Description,x=GeneRatio))+
  geom_point(aes(size=Count,color=pvalue))+ 
  facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')+
  scale_color_gradient(low = "red",high ="blue")+
  labs(color=expression(pvalue,size="Count"),  
       x="GeneRatio",y="GO term",title=paste0(target," VS ",Non_target))+ 
  theme_bw()
kegg_enrich <- arrange(kegg_enrich,GeneRatio)
kegg_enrich$Description <- factor(kegg_enrich$Description,levels=kegg_enrich$Description)
p3 <- ggplot(kegg_enrich,
             aes(y=Description,x=GeneRatio))+
  geom_point(aes(size=Count,color=pvalue))+ 
  scale_color_gradient(low = "red",high ="blue")+
  labs(color=expression(pvalue,size="Count"), 
       x="GeneRatio",y="",title=paste0(target," VS ",Non_target))+ 
  scale_y_discrete(labels=function(x) str_wrap(kegg_enrich$Description, width = 50)) +
  theme_bw()
