library(TCGAbiolinks)
library(SummarizedExperiment)
library(data.table)
library(dplyr)
library(stringr)
library(maftools)
library(ggplot2)
library(ggprism)
library(BSgenome.Hsapiens.UCSC.hg38)

query <- GDCquery(
  project = "TCGA-CESC",
  data.category = "Copy Number Variation",
  data.type = "Masked Copy Number Segment",              
  access = "open"
)
GDCdownload(query, files.per.chunk = 100)
data <- GDCprepare(query,save=T,save.filename = "01.CNV_TCGA_CESC.rda")
cc_cnv <- data[,-1]
cc_cnv <- cc_cnv[,c(6,1:5)]
tumor_seg <- cc_cnv[substr(cc_cnv$Sample,14,15)=="01",] 
names(tumor_seg) <- c("Sample","Chromosome","Start Position","End Position","Num markers","Seg.CN") 
hg_marker_file <- read.delim("snp6.na35.remap.hg38.subset.txt.gz")
Hg_marker_file <- hg_marker_file[hg_marker_file$freqcnv=="FALSE",]
Hg_marker_file <- Hg_marker_file[,1:3]
names(Hg_marker_file) <- c("Marker Name","Chromosome","Marker Position")
all.lesions <- "gistic_output/all_lesions.conf_90.txt"
amp.genes <- "gistic_output/amp_genes.conf_90.txt"
del.genes <- "gistic_output/del_genes.conf_90.txt"
scores.gis <- "gistic_output/scores.gistic"
coad.gistic = readGistic(gisticAllLesionsFile = all.lesions, 
                         gisticAmpGenesFile = amp.genes, 
                         gisticDelGenesFile = del.genes, 
                         gisticScoresFile = scores.gis, isTCGA = TRUE)
gisticChromPlot(gistic = coad.gistic,
                markBands = "all",
                ref.build = "hg38")
df <- data.frame(chromName = seqnames(BSgenome.Hsapiens.UCSC.hg38), 
                 chromlength = seqlengths(BSgenome.Hsapiens.UCSC.hg38))
df$chromNum <- 1:length(df$chromName) 
df <- df[1:22,] 
df$chromlengthCumsum <- cumsum(as.numeric(df$chromlength)) 
df$chormStartPosFrom0 <- c(0,df$chromlengthCumsum[-nrow(df)])
tmp_middle <- diff(c(0,df$chromlengthCumsum)) / 2
df$chromMidelePosFrom0 <- df$chormStartPosFrom0 + tmp_middle
scores <- read.table("gistic_output/scores.gistic",sep="\t",header=T,stringsAsFactors = F)
all_lesions <- read.table("gistic_output/all_lesions.conf_90.txt",sep="\t",header=T,stringsAsFactors = F)
chromID <- scores$Chromosome
scores$StartPos <- scores$Start + df$chormStartPosFrom0[chromID]
scores$EndPos <- scores$End + df$chormStartPosFrom0[chromID]
scores[scores$Type == "Del", "frequency"] <- scores[scores$Type == "Del", "frequency"] * -1
df$ypos <- rep(c(0.85,0.95),11) 
p1=ggplot(scores, aes(x = StartPos,y = frequency))+
  geom_area(aes(group=Type, fill=factor(Type,levels = c("Del","Amp"))))+
  scale_fill_manual(values = c("#00468b","#ed0000"),guide=guide_legend(reverse = T),name="Type")+
  geom_vline(data = df ,mapping=aes(xintercept=chromlengthCumsum),linetype=2)+
  geom_text(data = df,aes(x=chromMidelePosFrom0,y=ypos,label=chromName))+
  scale_x_continuous(expand = c(0,0),limits = c(0,2.9e9),name = NULL,labels = NULL)+
  scale_y_continuous(expand = c(0.02,0.02),limits = c(-1,1),guide = "prism_offset")+ 
  theme_prism()+
  theme(axis.ticks.x = element_blank(),
        axis.line.x = element_blank())