library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(getopt)
library(extrafont)
loadfonts(device = "pdf")

command = matrix(c('diff','d',1,"character",
		   'file','f',1,"character",
                   'output','o',1,"character")
                 ,byrow=TRUE, ncol=4)
#传参
args=getopt(command)
options(stringsAsFactors = FALSE)
outpath = args$output
diffile=args$diff
filename = args$file
result1 = paste0(outpath,filename,".png")
result2 = paste0(outpath,filename,".pdf")
Data<-read.table(diffile,header = T,check.names = F,sep='\t')
Data$threshold = factor(ifelse(Data$p_val < 0.05 & abs(Data$avg_log2FC) >= 0.25, ifelse(Data$avg_log2FC >= 0.25 ,'Up','Down'),'NoSignificant'),levels=c('Up','Down','NoSignificant'))
datap = filter(Data, p_val<0.05)
ranka<- datap[order(datap$avg_log2FC),]
ranka <- filter(ranka,avg_log2FC< -0.25)
rankc<- datap[order(-datap$avg_log2FC),]
rankc <- filter(rankc,avg_log2FC>0.25)

p=ggplot(Data,aes(x=avg_log2FC,y=-log10(p_val),color=threshold))+
geom_point()+
scale_color_manual(values=c("#DC143C","#00008B","#808080"))+
geom_text_repel(
data = rbind(ranka[1:min(nrow(ranka),5),],rankc[1:min(nrow(rankc),5),]),
aes(label = Gene),
size = 7,
family = "serif",
fontface = "italic",
segment.color = "black", show.legend = FALSE )+
theme_classic()+
theme(legend.title = element_blank())+
theme(legend.text=element_text(size=30, family = "serif"))+
theme(axis.title=element_text(size=30, family = "serif"),axis.text=element_text(size=30,color='black'))+
labs(x="log2(FoldChange)",y="-log10(p-Value)")+
#ylab('-log10(p-Value)')+
#xlab('log2(FoldChange)')+
geom_vline(xintercept=c(-0.25,0.25),lty=3,col="black",lwd=0.5) +
geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)
ggsave(p,file=result1,height = 8,width = 9)
ggsave(p,file=result2,height = 8,width = 9)
