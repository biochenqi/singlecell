library(Seurat)
library(ggplot2)
library(dplyr)
library(reshape2)
library(getopt)
library(extrafont)
loadfonts(device = "pdf")
# par(family="serif")
command = matrix(c('rds','r',1,"character",
                   'celltype','t',1,"character",
                   'title','s',1,"character",
		                'color','c',1,"character",
                   'output','o',1,"character")
                 ,byrow=TRUE, ncol=4)
args=getopt(command)

mycolors <-c("#d6604d","#5a8ca8","#d8869d","#689e45","#dfc27d","#93c1b6")
# colors = strsplit(args$color,split=",")[[1]]
# mycolors <- c("CD4+TN" = "blue", "CD4+TRM" = "green", "TFH" = "red", "TH17" = "purple", "Treg" = "orange","Undefined"="gray")
# mycolors <- c("CD4+TN" = "#dfc27d", "CD4+TRM" = "#689e45", "TFH" = "#d8869d", "TH17" = "#5a8ca8", "Treg" = "#d6604d","Undefined"="gray")

if(!is.null(args$color)){
  # 读取文件
  color_df <- read.table(args$color, header = TRUE, stringsAsFactors = FALSE,sep='\t',comment.char="")

  # 将数据框的列转换为向量
  mycolors <- setNames(as.vector(color_df$col), color_df$type)

  # 添加"Undefined"="gray"到向量
  mycolors["Undefined"] = "gray"
}


output = args$output
if(!dir.exists(output)){
  dir.create(output)
}

type = read.table(args$celltype,header = T,sep = "\t",row.names = 1)
seuset = readRDS(args$rds)
seuset = AddMetaData(seuset,type,col.name = "Type")
# levels(seuset$Type) <- mycolors



p = DimPlot(seuset,label = F,group.by="Type",raster = F,label.size=7)+
  labs(title = args$title)+
  coord_cartesian(clip="off")+
  theme(text=element_text(size=30, family = "serif"),axis.text=element_text(size=30, family = "serif"),plot.title = element_text(hjust=0.5,size=40, family = "serif"))+
  scale_color_manual(values=mycolors)

ggsave(paste0(output,"/UMAPPlot_Type.png"),p,width = 9,height = 8)
ggsave(paste0(output,"/UMAPPlot_Type.pdf"),p,width = 9,height = 8)
