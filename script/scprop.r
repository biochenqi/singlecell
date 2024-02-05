library(dplyr)
library(getopt)
library(ggplot2)

command = matrix(c('summerycell','t',1,"character", 
                   'cellinfo','s',1,"character", 
                   'sample_color','S',2,"character", 
                   'type_color','T',2,"character",
                   'sheight','h',2,"numeric",
                   'swidth','w',2,"numeric",
                   'theight','H',2,"numeric",
                   'twidth','W',2,"numeric",
                   'output', 'o', 1, "character")
                 ,byrow = TRUE, ncol = 4)

args = getopt(command)
options(stringsAsFactors = FALSE)

#读取细胞类型文件
summerycell <- read.table(args$summerycell, sep = "\t",header = T,comment.char = "")
colnames(summerycell) <- c("Cell","Celltypes")

#读取样本类型文件
cellinfo <- read.table(args$cellinfo, sep = "\t",header = T)
colnames(cellinfo) <- c("Cell", "Sample")

#合并数据
dat <- merge(summerycell, cellinfo, by = "Cell")

#读取样本颜色信息
if (!is.null(args$sample_color)) {
  samplecolor <- read.table(args$sample_color, sep = "\t",header = T,comment.char = "")
  if(ncol(samplecolor)==1){
    colnames(samplecolor) <- "Sample"
  }else if(ncol(samplecolor)==2){
  colnames(samplecolor) <- c("Sample", "SampleColor")
  dat <- left_join(dat, samplecolor, by = "Sample")
  }else{
    warning("Check columns of sample color file")
  }
}

#读取细胞类型颜色信息
if (!is.null(args$type_color)) {
  typecolor <- read.table(args$type_color, sep = "\t",header = T,comment.char = "")
  if(ncol(typecolor)==1){
    colnames(typecolor) <- "Celltypes"
  }else if(ncol(typecolor)==2){
    colnames(typecolor) <- c("Celltypes", "TypeColor")
    dat <- left_join(dat, typecolor, by = "Celltypes")
  }else{
    warning("Check columns of type color file")
  }

}

sample_num <- cellinfo %>% group_by(Sample) %>% count()
colnames(sample_num) <- c("Sample", "Count")
typesample_num <- dat %>% group_by(Celltypes,Sample) %>% count()
df = left_join(typesample_num,sample_num,by="Sample") %>% mutate(prop=n/Count)
write.table(df,"prop.data.txt",quote = F,sep = "\t",col.names = T,row.names = F)
if(!is.null(args$sample_color)){
  if(all(samplecolor$Sample %in% unique(df$Sample))){
    df$Sample=factor(df$Sample,levels = rev(samplecolor$Sample))
  }else{
    warning("Check sample names of sample color file")
  }
}
if(!is.null(args$type_color)){
  if(all(typecolor$Celltypes %in% unique(df$Celltypes))){
    df$Celltypes=factor(df$Celltypes,levels = rev(typecolor$Celltypes))
  }else{
    warning("Check Celltypes of Type color file")
  }
}
#设置颜色
mycolors <-c("#d6604d","#5a8ca8","#d8869d","#689e45","#dfc27d","#93c1b6","#fdae61","#b15928",
             "#e7298a","#cab2d6","#6a339a","#003c30","#d53e4f","#3288bd","#abdda4","#66c2a5",
             '#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
              '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
              '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
              '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
              '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
              '#968175')
sampleheight=min(1.5*length(unique(df$Sample)),8)
typeheight=min(1.5*length(unique(df$Celltypes)),8)
samplewidth=2*sampleheight
typewidth=2*typeheight

if(!is.null(args$sheight)){sampleheight=args$sheight}
if(!is.null(args$swidth)){samplewidth=args$swidth}
if(!is.null(args$theight)){typeheight=args$theight}
if(!is.null(args$twidth)){typewidth=args$twidth}

#设置字体大小
font_size=20

p_sample = ggplot(df, aes(x = Sample, y = n, fill = Celltypes)) +
  geom_col(position = "fill",width=0.7) + labs( x = "", y = "") +
  #coord_flip() + 
  theme_classic() + guides(fill = guide_legend(reverse = T)) + 
  theme(axis.line = element_line(size = 0.8), axis.ticks = element_line(size = 0.8), 
        axis.ticks.length = unit(1.5, "mm"), axis.text = element_text(size = font_size, colour = "black",family = "serif"), 
        legend.text = element_text(colour = "black", size = font_size ,family = "serif"),
        legend.title = element_text(colour = "black", size = font_size, family = "serif"))

if(!is.null(args$type_color)){
  if(ncol(typecolor)==2){
    p_sample = p_sample + scale_fill_manual(values = rev(typecolor$TypeColor))
  }else{
    p_sample = p_sample + scale_fill_manual(values = rev(mycolors[1:length(unique(df$Celltypes))]))
  }
}else{
  p_sample = p_sample + scale_fill_manual(values = rev(mycolors[1:length(unique(df$Celltypes))]))
}

ggsave(paste0(args$output,"/PropPlot_bySample.pdf"),p_sample,height = sampleheight,width = samplewidth)
ggsave(paste0(args$output,"/PropPlot_bySample.png"),p_sample,height = sampleheight,width = samplewidth)


p_type = ggplot(df, aes(x = Celltypes, y = n, fill = Sample)) +
  geom_col(position = "fill") + labs( x = "", y = "") +
  #coord_flip() + 
  theme_classic() + guides(fill = guide_legend(reverse = T)) + 
  theme(axis.line = element_line(size = 0.8), axis.ticks = element_line(size = 0.8), 
        axis.ticks.length = unit(1.5, "mm"), axis.text = element_text(size = font_size+3, colour = "black"), 
        legend.text = element_text(colour = "black", size = font_size+3),
        legend.title = element_text(colour = "black", size = font_size+3))

if(!is.null(args$sample_color)){
  if(ncol(samplecolor)==2){
    p_type = p_type + scale_fill_manual(values = rev(samplecolor$SampleColor))
  }else{
    p_type = p_type + scale_fill_manual(values = rev(mycolors[1:length(unique(df$Sample))]))
  }
}else{
  p_type = p_type + scale_fill_manual(values = rev(mycolors[1:length(unique(df$Sample))]))
}

ggsave(paste0(args$output,"/PropPlot_byCelltypes.pdf"),p_type,height = typeheight,width = typewidth)
ggsave(paste0(args$output,"/PropPlot_byCelltypes.png"),p_type,height = typeheight,width = typewidth)
