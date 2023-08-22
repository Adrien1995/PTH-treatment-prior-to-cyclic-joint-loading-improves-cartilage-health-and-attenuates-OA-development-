rm(list=ls())

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("enrichplot")
#BiocManager::install("clusterProfiler")
#BiocManager::install("DOSE")
#BiocManager::install("ggnewscale")
#BiocManager::install("pathview")


library(DOSE)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggnewscale)
library(dplyr)
#BiocManager::install(organism, character.only = TRUE)
library("org.Mm.eg.db", character.only = TRUE)
#BiocManager::install("BiocParallel")
library("BiocParallel")


# Setting desired species
organism = "org.Mm.eg.db"

i = "file_CPM3.txt"

df <- read.delim(i)
ttl <- strsplit(i,"_CPM")

original_gene_list <- df$logFC

names(original_gene_list) <- df$Gene

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

#CC, MF, BP
#ALL

gse2 <- gseGO(geneList=gene_list, 
              ont ="BP", 
              keyType = "SYMBOL", 
              nPerm = 10000, 
              minGSSize = 3, 
              maxGSSize = 800, 
              pvalueCutoff = 0.05, 
              verbose = TRUE, 
              OrgDb = org.Mm.eg.db, 
              pAdjustMethod = "fdr")


write.csv(gse2@result,"gsearesults_....csv")

