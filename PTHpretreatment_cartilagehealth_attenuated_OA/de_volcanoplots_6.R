rm(list=ls())
#BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
library(ggplot2)

De_direc = getwd()
file_list <- list.files(path=DE_direc)


for (i in file_list){
  dat <- read.delim(i)
  ttl <- strsplit(i,"_CPM3")
  dat <- dat[!duplicated(dat$Gene),]
  row.names(dat) <- dat[,1]
  dat <- dat[,-1]
  sortdat <- dat[order(dat$logFC),]
  #top 10 logFC (5 pos, 5 neg) genes
  #topgenes= rownames(sortdat)[c(1:10,14498:145081)]
  a = nrow(dat)
  b = a-3
  topgenes= rownames(sortdat)[c(1:3,a:b)]
  #tiff("test.tiff", units="in", width=5, height=5, res=300)
  tiff(paste(ttl[[1]][[1]],"_volcano.tiff"), units="in", width=5, height=5, res=300)
  # insert ggplot code
  
  plot(EnhancedVolcano(dat,
                       lab= rownames(dat),
                       selectLab= topgenes,
                       xlab = bquote(~Log[2] ~ "fold change"),
                       ylab = bquote(~-Log[10] ~ "FDR"),
                       x= 'logFC',
                       y= 'FDR',
                       title= ttl[[1]][1],
                       axisLabSize = 20,
                       legendPosition = "none",
                       cutoffLineCol = "red2",
                       vline = 0,
                       vlineCol = ("black"), 
                       vlineType = "solid",
                       xlim = c(-10,10),
                       ylim = c(0,8),
                       
                       # hline = 0,
                       # hlineCol = ("black"), 
                       # hlineType = "solid",
                       #legendPosition = "right",
                       subtitle= "Top 10 DEGs labeled",
                       pCutoff= 0.05,
                       shape= 20,
                       colCustom= NULL,
                       colAlpha= 1,
                       FCcutoff= 2,
                       border = "full",
                       borderColour = "white",
                       col = c("grey30", "grey30", "grey30", "red2"),
                       pointSize= c(ifelse(rownames(dat) %in% topgenes,5,2)),
                       labSize= 4,
                       labCol = "royalblue",
                       labFace= 'bold',
                       drawConnectors= TRUE,
                       widthConnectors= 0.75,
                       endsConnectors = 'last'))
  dev.off()
}
