
rm(list=ls())

library(dplyr)

#BiocManager::install("edgeR")
library(edgeR)
#install.packages("DESeq2")

#BiocManager::install("DeSeq2")
library(DESeq2)
if (!require("NMF")) {
  install.packages("NMF", dependencies = TRUE)
  library(NMF)
}

if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}

if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}


#DE_file <- read.delim("DE_VV1loadedvsctrlcartilage_CPM3.txt")
count_file <- read.delim("VV1_VA1_PV1_PA1_cart_norm.txt")

#tagged_genes <- c("Sox9","Col2a1","Acan","Runx2","Col10a1","Mmp13","Spp1","Epas1","Sost","Sox5","Sox6","Fgfr3","Gdf5","Col27a1","Col9a1","Fto","Klhdc5","Il1rn",
#                 "Bmp2","H2-Ab1","Il1b","Mmp13","Il6","Il1a","Il1b","Osm", "Ifngr", "Mc1r")

tagged_genes <- c("Sox9","Col2a1","Acan","Runx2","Col10a1","Mmp13","Spp1","Epas1","Sost","Sox5","Sox6","Fgfr3","Gdf5","Col27a1","Col9a1") 

count_file <- count_file[!duplicated(count_file$Gene),]


grp <- c("VV1","VA1","PV1","PA1")
for (i in grp){
  
  datvv <- count_file %>% dplyr::select(matches(i))
  
  datvv <- datvv[,order(colnames(datvv),decreasing=TRUE)]
  
  rownames(datvv) <- count_file$Gene
  
  datvv <- datvv %>% filter(rownames(datvv) %in% tagged_genes)
  
  datvv_matrix1 <- as.matrix(datvv)
  
  tiff(paste(i,"_heatmap_cart.tiff"), units="in", width=10, height=12, res=300)
  #aheatmap(datvv_matrix1, col=colorRampPalette(c("#FDE725FF", "#27AD81FF", "#440154FF"))(50), scale="row", Colv = NA, fontsize = 18)
  aheatmap(datvv_matrix1, col=brewer.pal(11,"RdBu"), scale="row", Colv = NA, fontsize = 32)
  dev.off()
  }










