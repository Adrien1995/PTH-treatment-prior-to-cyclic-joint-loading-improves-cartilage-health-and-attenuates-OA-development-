#comparing differential expression 

rm(list=ls())

#install.packages("VennDiagram")
library(VennDiagram)
library(dplyr)

vv1_cart <- read.delim("vv1cart_detable.txt")
vv2_cart <- read.delim("vv2cart_detable.txt")

vv1_cort <- read.delim("vv1_cort_detable.txt")
vv2_cort <- read.delim("vv2_cort_detable.txt")

vv1_canc <- read.delim("vv1_canc_detable.txt")
vv2_canc <- read.delim("vv2_canc_detable.txt")

vv1_LN <- read.delim("vv1_LN_detable.txt")
vv2_LN <- read.delim("vv2_LN_detable.txt")


vv1_cart <- vv1_cart %>% filter(TotalgenesFCover2 == 1)
vv2_cart <- vv2_cart %>% filter(TotalgenesFCover2 == 1)

vv1_cort <- vv1_cort %>% filter(TotalgenesFCover2 == 1)
vv2_cort <- vv2_cort %>% filter(TotalgenesFCover2 == 1)

vv1_canc <- vv1_canc %>% filter(TotalgenesFCover2 == 1)
vv2_canc <- vv2_canc %>% filter(TotalgenesFCover2 == 1)

vv1_LN <- vv1_LN %>% filter(TotalgenesFCover2 == 1)
vv2_LN <- vv2_LN %>% filter(TotalgenesFCover2 == 1)


# Chart

vv1_cart_genes <- vv1_cart$Gene
vv2_cart_genes <- vv2_cart$Gene

vv1_cort_genes <- vv1_cort$Gene
vv2_cort_genes <- vv2_cort$Gene

vv1_canc_genes <- vv1_canc$Gene
vv2_canc_genes <- vv2_canc$Gene

vv1_LN_genes <- vv1_LN$Gene
vv2_LN_genes <- vv2_LN$Gene



color4 <- c("purple") 
ct4 <- adjustcolor(color4, alpha.f = 0.2)
ct4b <- adjustcolor(color4, alpha.f = 0.6)
ct4a <- adjustcolor(color4, alpha.f = 0.9)


color3 <- c("blue") 
ct3 <- adjustcolor(color3, alpha.f = 0.2) 
ct3b <- adjustcolor(color3, alpha.f = 0.6) 
ct3a <- adjustcolor(color3, alpha.f = 0.9) 


color2 <- c("red3") 
ct2 <- adjustcolor(color2, alpha.f = 0.2)
ct2b <- adjustcolor(color2, alpha.f = 0.6)
ct2a <- adjustcolor(color2, alpha.f = 0.85) 


color1 <- c("black") 
ct1 <- adjustcolor(color1, alpha.f = 0.2) 
ct1b <- adjustcolor(color1, alpha.f = 0.6)
ct1a <- adjustcolor(color1, alpha.f = 0.9) 

#("lightgreen", 'turquoise', 'orange', 'purple')

color5 <- c("lightgreen") 
ct5 <- adjustcolor(color5, alpha.f = 0.2) 
ct5b <- adjustcolor(color5, alpha.f = 0.6)
ct5a <- adjustcolor(color5, alpha.f = 0.9) 

color6 <- c("turquoise") 
ct6 <- adjustcolor(color6, alpha.f = 0.2) 
ct6b <- adjustcolor(color6, alpha.f = 0.6)
ct6a <- adjustcolor(color6, alpha.f = 0.9) 

color7 <- c("orange") 
ct7 <- adjustcolor(color7, alpha.f = 0.2) 
ct7b <- adjustcolor(color7, alpha.f = 0.6)
ct7a <- adjustcolor(color7, alpha.f = 0.9) 

color8 <- c("purple") 
ct8 <- adjustcolor(color8, alpha.f = 0.2) 
ct8b <- adjustcolor(color8, alpha.f = 0.6)
ct8a <- adjustcolor(color8, alpha.f = 0.9) 




#### 1 and 2 wk comaprison 

venn.diagram(
  x = list(vv1_cart_genes,vv2_cart_genes),
  category.names = c("wk1" , "wk2"),
  filename = 'cart1&2wks_venn.png',
  output=TRUE,
  col=c("black", "black"),
  fill = c(ct4, ct4a),
  #font.size = 15,
  cex =2.0,
  cat.cex = 2.0,
  cat.default.pos = "outer",
  cat.fontfamily="sans",
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055)
  
)

venn.diagram(
  x = list(vv1_cort_genes,vv2_cort_genes),
  category.names = c("wk1" , "wk2"),
  filename = 'cort1&2wks_venn.png',
  output=TRUE,
  col=c("black", "black"),
  fill = c(ct4,ct4a),
  #font.size = 15,
  cex =2.0,
  cat.cex = 2.0,
  cat.default.pos = "outer",
  cat.fontfamily="sans",
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055)
)

venn.diagram(
  x = list(vv1_canc_genes,vv2_canc_genes),
  category.names = c("wk1" , "wk2"),
  filename = 'canc1&2wks_venn.png',
  output=TRUE,
  col=c("black", "black"),
  fill = c(ct4,ct4a),
  #font.size = 15,
  cex =2.0,
  cat.cex = 2.0,
  cat.default.pos = "outer",
  cat.fontfamily="sans",
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055)
)

venn.diagram(
  x = list(vv1_LN_genes,vv2_LN_genes),
  category.names = c("wk1" , "wk2"),
  filename = 'LN1&2wks_venn.png',
  output=TRUE,
  col=c("black", "black"),
  fill = c(ct4,ct4a),
  #font.size = 15,
  cex =2.0,
  cat.cex = 2.0,
  cat.default.pos = "outer",
  cat.fontfamily="sans",
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055)
)


###cross tissue at 1wk 

venn.diagram(
  x = list(vv1_cart_genes,vv1_cort_genes, vv1_canc_genes, vv1_LN_genes),
  category.names = c("cart" , "cort", "canc", "LN"),
  filename = 'crosstissue_1wk_venn.png',
  output=TRUE,
  col=c("lightgreen", 'turquoise', 'orange', 'purple'),
  fill = c(ct5,ct6, ct7, ct8),
  #font.size = 15,
  cex =2.0,
  cat.cex = 2.0
)

###cross tissue at 1wk 

venn.diagram(
  x = list(vv2_cart_genes, vv2_cort_genes, vv2_canc_genes, vv2_LN_genes),
  category.names = c("cart" , "cort", "canc", "LN"),
  filename = 'crosstissue_2wks_venn.png',
  output=TRUE,
  col=c("lightgreen", 'turquoise', 'orange', 'purple'),
  fill = c(ct5,ct6, ct7, ct8),
  #font.size = 15,
  cex =2.0,
  
  cat.cex = 2.0
)



