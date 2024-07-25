exp <- read.csv("expdataTumor.csv", header=T, check.names=F, row.names=1)
RSI_gene <- read.table("RSI_genelist.txt")
RSI_gene <- RSI_gene[,1]
RSI_exp <- exp[RSI_gene,]
RSI_exp <- as.data.frame(t(RSI_exp))
RSI_exp$score = -0.0098009*RSI_exp$AR + 0.0128283*RSI_exp$JUN + 0.0254552*RSI_exp$STAT1 - 0.0017589*RSI_exp$PRKCB - 0.0038171*RSI_exp$RELA + 0.1070213*RSI_exp$ABL1- 0.0002509*RSI_exp$SUMO1 - 0.0092431*RSI_exp$PAK2 - 0.0204469*RSI_exp$HDAC1 - 0.0441683*RSI_exp$IRF1
group <- read.table("cluster.txt", header = T, row.names = 1)
sameSample=intersect(row.names(RSI_exp), row.names(group))
group=group[sameSample,,drop=F]
RSI_exp=RSI_exp[sameSample,,drop=F]
all(rownames(group)==rownames(RSI_exp)) 
rt=cbind(group, RSI_exp)
rt$cluster=factor(rt$cluster, levels=c("young","medium", "old"))
p <- ggstatsplot::ggbetweenstats(
  data = rt, 
  x = cluster, 
  y = score,
  type = "np",    
  mean.ci = TRUE,    
  pairwise.display = "all",  
  p.adjust.method = "fdr",  
  title = "RSI",
  caption = "",
  k = 3L, 
  centrality.plotting = F,
  point.args = list(position = ggplot2::position_jitterdodge(jitter.width = 0.3,dodge.width = 0.1),alpha =1, 
                    size = 2, stroke = 0, na.rm = TRUE),
  boxplot.args = list(aes(fill=cluster),notch=F,width = 0.25, size=0.4,alpha = 0.4, na.rm = TRUE),
  violin.args = list(aes(fill=cluster),width = 0.95, size=0.4,alpha = 0.2, na.rm = TRUE),
  ggsignif.args = list(textsize = 2, tip_length = 0.01,size=0.4, na.rm = TRUE),
  ggtheme = theme_classic())+
  theme(axis.line = element_line(linewidth = 0.4),
        axis.title.y.left  = element_text(size = 12),
        axis.title.y.right  = element_text(size = 5),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 12),plot.subtitle = element_text(size = 5))+ 
  ggplot2::scale_color_manual(values = c("#00468B","#925E9F","#759EDD"))+
  ggplot2::scale_fill_manual(values = c("#00468B","#925E9F","#759EDD"))