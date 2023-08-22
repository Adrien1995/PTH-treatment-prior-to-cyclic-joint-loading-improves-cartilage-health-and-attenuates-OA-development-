#computing differential expression
rm(list=ls())

library(dplyr)


DE_direc = getwd()
file_list <- list.files(path=DE_direc)
for (i in file_list){
  pa2cart <- read.delim(i)
  ttl <- strsplit(i,"_CPM")
  
  pa2cart$positivefoldchange <- 2**(pa2cart$logFC)
  pa2cart$negativefoldchange <- -1/(pa2cart$positivefoldchange)
  pa2cart <- pa2cart %>% mutate(foldchange = if_else(positivefoldchange > 1, positivefoldchange,negativefoldchange))
  pa2cart <- pa2cart %>% mutate(FDRsig = if_else(FDR < 0.05, 1, 0))
  pa2cart <- pa2cart %>% mutate(absFC = if_else(FDRsig > 0, abs(foldchange), 0))
  
  pa2cart <- pa2cart %>% mutate(TotalgenesFCover2 = if_else(absFC > 2, 1, 0))
  pa2cart$totalsiggenes <- sum(pa2cart$TotalgenesFCover2)
  pa2cart <- pa2cart %>% mutate(ActualFCofsiggenes = if_else(TotalgenesFCover2> 0,foldchange, 0))
  pa2cart <- pa2cart %>% mutate(Column4 = if_else(ActualFCofsiggenes> 0, 1, 0))
  
  pa2cart$totalupgenes <- sum(pa2cart$Column4)
  pa2cart$totaldowngenes <-pa2cart$totalsiggenes - pa2cart$totalupgenes
  
  write.table(pa2cart, sep = "\t", paste(ttl[[1]][[1]],"_detable.txt"))
  
}





