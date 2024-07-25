library(ReactomePA)
library(reactome.db)
library(org.Hs.eg.db)
library(ggplot2)

data <- read.csv("cpmResultDiff.csv")
sig <- data$name
sig_ID <- bitr(sig,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = "org.Hs.eg.db")
eReac <- enrichPathway(gene = sig_ID$ENTREZID,organism = 'human',pvalueCutoff = 1,qvalueCutoff = 1)
bardata <- eReac@result 
rownames(bardata) <- bardata$Description
bardata <- bardata[sort(bardata$Description,decreasing = TRUE),]
bardata$Description <- factor(bardata$Description,levels=bardata$Description)
barplot <- ggplot(bardata,aes(y=Description,x=Count))+
  geom_bar(stat = "identity",aes(fill=pvalue))+
  scale_fill_gradientn(colours =  c("#4D6381", "#BAC2CC"))+
  theme_bw()
eReacx <- setReadable(eReac, 'org.Hs.eg.db', 'ENTREZID')
cnetplot <-  cnetplot(eReacx, 
                      circular = TRUE, 
                      colorEdge = TRUE,
                      categorySize="pvalue",
                      showCategory = bardata$Description,
                      node_label = "gene",
                      layout = 'kk') 