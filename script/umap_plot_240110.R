library(tidyverse)
library(ggplot2)
library(Seurat)
library(dplyr)
library(getopt)

command = matrix(c('rds','r',1,"character",
                   'genes','g',1,"character",
                  'output','o',1,"character")
                 ,byrow=TRUE, ncol=4)
args=getopt(command)


seuset = readRDS(args$rds)

output = args$output
if(!dir.exists(output)){
  dir.create(output)
}

#gene = read.table(args$genes,header = F,sep = "\t")
#genelist =gene$V1
genelist = strsplit(args$genes,split=",")[[1]]
for(i in genelist){
  if(i %in% rownames(seuset)){
    print(paste0("Plot ",i))
    p = FeaturePlot(seuset,i,label=T)+theme(plot.title = element_text(size=30))
    ggsave(paste0(output,"/",i,".UMAP.png"),p,width = 8.5,height = 8)
    p = VlnPlot(object = seuset, features = i,pt.size=0)
    ggsave(paste0(output,"/",i,".Violin.png"),p,width = 8.5,height = 8)
  }else{
    print(paste0("Cannot find gene ",i))
  }
}
