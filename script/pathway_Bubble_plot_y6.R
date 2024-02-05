library(ggplot2)
library("getopt")
library(extrafont)
loadfonts(device = "pdf")
##################Rscript GO_KEGG_Bubble_plot.R -i infile -n name -t GO/KEGG/Reactome
command=matrix(c("file","i",1,"character",
                                "name","n",1,"character",
				"type","t",1,"character",
                                "help","h",0,"logical"),byrow=T,ncol=4)
args=getopt(command)
if (!is.null(args$help) || is.null(args$file)|| is.null(args$name)) {
    cat(paste(getopt(command, usage = T), "\n"))
    q()
}

pathway = read.table(args$file,sep="\t",header=T)
diff_gene<-as.numeric(unlist(strsplit(as.character(pathway$GeneRatio[1]),"/"))[2])
pathway$GeneRatio<-pathway$Count/diff_gene

p<-ggplot(pathway, aes(x = pathway$GeneRatio,y = reorder(Term,-pathway$PValue)))+theme_bw()+
theme(panel.background = element_blank(),panel.border=element_rect(fill="transparent", color="black"))+
geom_point(aes(colour=PValue,size=Count))+scale_colour_gradientn(colours=c("red","green")) +
expand_limits(color=seq(0, 0.05, by=0.01))+
ggtitle(paste(paste(args$name,args$type,sep="."),"_Enrich",sep="")) + xlab("GeneRatio") +ylab("")+
theme(title=element_text(size=20, family = "serif"),axis.text.y=element_text( color="black", size=20, family = "serif"),axis.text.x=element_text( color="black", size=16, family = "serif"),legend.title = element_text(size=18,face="bold", family = "serif"),legend.text=element_text(size = 15, family = "serif"),plot.title=element_text(hjust=0.5, family = "serif"),legend.key=element_blank())

ggsave(p,file=paste(paste(args$name,args$type,sep="."),".png",sep=""),  device="png",  dpi=700, height=7,width=13)
ggsave(p,file=paste(paste(args$name,args$type,sep="."),".pdf",sep=""),  device="pdf", height=7,width=13)

