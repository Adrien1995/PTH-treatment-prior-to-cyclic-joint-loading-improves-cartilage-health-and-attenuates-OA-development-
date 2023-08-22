#fitting glm models

rm(list=ls())

#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("edgeR")
#BiocManager::install("DESeq2")

library(edgeR)
library(DESeq2)
library(readxl)


x <- read.delim("Gene_count_compiled.txt")
x1 <- x[-c(1:4),]
rownames(x1) <- x1[,1] 
x1 <- x1[,-1]

#group <- factor(c(2,1,2,1,2,1,2,1,2,1,2,1,2,1,4,3,4,3,4,3,4,3,4,3,4,3,4,3,4,3,6,5,6,5,6,5,6,5,6,5,6,5,6,5,6,5,8,7,8,7,8,7,8,7,8,7,7,8,7)) 
#plate <- factor(c(2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)) 

#group <- factor(c(2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,4,3,4,3,4,3,4,3,4,3,4,3,4,3,4,3,5,6,5,6,5,6,5,5,6,6,6,7,7,7,8,7,8,7,8,7,8,7,8)) 
#plate <- factor(c(2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3)) 

group <- factor(c( 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2)) 
plate <- factor(c( 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2)) 
#plate factor start naming at 2


dds <- DESeqDataSetFromMatrix(x1, DataFrame(group), ~group)
rld <- vst(dds)
plotPCA(rld, intgroup = "group")


#sample<- factor(c(7,14,21,28,35,42,49,56,63,70,77,84,91,98,231,238,245,252,259,266,273,280,287,294,301,308,315,322,329,336,119,126,133,140,147,154,161,168,175,182,189,196,203,210,217,224,343,350,357,364,371,378,385,392,399,420,434,441,448)) 

#sample <- factor(c(455,462,469,476,483,490,497,504,511,518,525,532,539,546,553,560,567,574,581,588,595,602,609,616,623,630,637,644,651,658,665,672,B62,B63,B64,B65,B66,B67,B68,B70,B71,B73,B75,B78,B80,B82,B83,B84,B85,B86,B87,B88,B89,B90,B91)) 



#y <-DGEList(counts=x1, group=group)
#y <- calcNormFactors(y)
#pdf("PCA_V0_PO_cart.pdf")
#plotMDS(y, labels=sample)
#dev.off()

#normalizing 
dds <- DESeqDataSetFromMatrix(x1, DataFrame(group), ~group)
ddsEstimate <- estimateSizeFactors(dds)
table_counts_normalized <- counts(ddsEstimate, normalized=TRUE)
write.table(table_counts_normalized, "norm_file.txt", sep="\t")



#####good to output, but not best normalization technique
y<-DGEList(counts=x1, group=group)
y <- calcNormFactors(y)
keep <-rowSums(cpm(y)>=3) >=2
y<-y[keep,] 
d <- cpm(y)
write.table(d, "file_CPM3.txt", sep="\t")
################



####DE 
y<-DGEList(counts=x1, group=group)
y <- calcNormFactors(y)
design<-model.matrix(~0+group)
y <- estimateGLMCommonDisp(y,design)
y <- estimateGLMTrendedDisp(y,design)
y <- estimateGLMTagwiseDisp(y,design)
fit<-glmFit(y,design)


#checking treatment effect 
lrt.2vs1 <- glmLRT(fit, contrast=c(-1,1))
top2v1 <- topTags(lrt.2vs1, n=25000)
write.table(top2v1, "DE_POvsVO_cartilage_CPM3.txt", sep="\t")




##############################################################

lrt.2vs1 <- glmLRT(fit, contrast=c(0,0,1,-1,0,0,0,0))
top2v1 <- topTags(lrt.2vs1, n=25000)
write.table(top2v1, "DE_PV2loadedvsctrlcartilage_CPM3.txt", sep="\t")



lrt.2vs1 <- glmLRT(fit, contrast=c(0,0,0,0,1,-1,0,0))
top2v1 <- topTags(lrt.2vs1, n=25000)
write.table(top2v1, "DE_VA2loadedvsctrlcartilage_CPM3.txt", sep="\t")

lrt.2vs1 <- glmLRT(fit, contrast=c(0,0,0,0,0,0,1,-1))
top2v1 <- topTags(lrt.2vs1, n=25000)
write.table(top2v1, "DE_PA2loadedvsctrlcartilage_CPM3.txt", sep="\t")

#comparing treatment effect (terms used not comptely accurate, jsut to help me remember)

#pth effect during laoding 

lrt.2vs1 <- glmLRT(fit, contrast=c(1,0,-1,0,0,0,0,0))
top2v1 <- topTags(lrt.2vs1, n=25000)
write.table(top2v1, "DE_VV2loaded_PV2loaded_cartilage_CPM3.txt", sep="\t")

#alendronate effect during loading

lrt.2vs1 <- glmLRT(fit, contrast=c(1,0,0,0,-1,0,0,0))
top2v1 <- topTags(lrt.2vs1, n=25000)
write.table(top2v1, "DE_VV2loaded_VA2loaded_cartilage_CPM3.txt", sep="\t")


#combined PTH ALN effect during loading

lrt.2vs1 <- glmLRT(fit, contrast=c(1,0,0,0,0,0,-1,0))
top2v1 <- topTags(lrt.2vs1, n=25000)
write.table(top2v1, "DE_VV2loaded_PA2loaded_cartilage_CPM3.txt", sep="\t")

#is alendroante + PTH additive as compared to VV1 compared to VA1 or is it different?
lrt.2vs1 <- glmLRT(fit, contrast=c(0,0,1,0,0,0,-1,0))
top2v1 <- topTags(lrt.2vs1, n=25000)
write.table(top2v1, "DE_PV2loaded_PA2loaded_cartilage_CPM3.txt", sep="\t")

#is alendroante better after PTH pretreatment as compared to VA1 or is it different?
lrt.2vs1 <- glmLRT(fit, contrast=c(0,0,0,0,1,0,-1,0))
top2v1 <- topTags(lrt.2vs1, n=25000)
write.table(top2v1, "DE_VA2loaded_PA2loaded_cartilage_CPM3.txt", sep="\t")

#PV1 vs VA1
lrt.2vs1 <- glmLRT(fit, contrast=c(0,0,1,0,-1,0,0,0))
top2v1 <- topTags(lrt.2vs1, n=25000)
write.table(top2v1, "DE_PV2loaded_VA2loaded_cartilage_CPM3.txt", sep="\t")




#controls
#pth effect ctrl

lrt.2vs1 <- glmLRT(fit, contrast=c(0,1,0,-1,0,0,0,0))
top2v1 <- topTags(lrt.2vs1, n=25000)
write.table(top2v1, "DE_VV2ctrl_PV2ctrl_cartilage_CPM3.txt", sep="\t")

#alendronate effect during loading

lrt.2vs1 <- glmLRT(fit, contrast=c(0,1,0,0,0,-1,0,0))
top2v1 <- topTags(lrt.2vs1, n=25000)
write.table(top2v1, "DE_VV2ctrl_VA2ctrl_cartilage_CPM3.txt", sep="\t")

#combined PTH ALN effect during loading

lrt.2vs1 <- glmLRT(fit, contrast=c(0,1,0,0,0,0,0,-1))
top2v1 <- topTags(lrt.2vs1, n=25000)
write.table(top2v1, "DE_VV2ctrl_PA2ctrl_cartilage_CPM3.txt", sep="\t")

#is alendroante + PTH additive as compared to VV1 compared to VA1 or is it different?
lrt.2vs1 <- glmLRT(fit, contrast=c(0,0,0,1,0,0,0,-1))
top2v1 <- topTags(lrt.2vs1, n=25000)
write.table(top2v1, "DE_PV2ctrl_PA2ctrl_cartilage_CPM3.txt", sep="\t")

#is alendroante better after PTH pretreatment as compared to VA1 or is it different?
lrt.2vs1 <- glmLRT(fit, contrast=c(0,0,0,0,0,1,0,-1))
top2v1 <- topTags(lrt.2vs1, n=25000)
write.table(top2v1, "DE_VA2ctrl_PA2ctrl_cartilage_CPM3.txt", sep="\t")

#PV1 vs VA1
lrt.2vs1 <- glmLRT(fit, contrast=c(0,0,0,1,0,-1,0,0))
top2v1 <- topTags(lrt.2vs1, n=25000)
write.table(top2v1, "DE_PV2ctrl_VA2ctrl_cartilage_CPM3.txt", sep="\t")





#alternative permutations
lrt.2vs1 <- glmLRT(fit, contrast=c(1,-1,0,0,0,0,0,0))
top2v1 <- topTags(lrt.2vs1, n=25000)
write.table(top2v1, "DE_VV1loadedvsctrlcartilage_CPM3.txt", sep="\t")

lrt.2vs1 <- glmLRT(fit, contrast=c(1,-1,0,0,0,0,0,0))
top2v1 <- topTags(lrt.2vs1, n=25000)
write.table(top2v1, "DE_VV1loadedvsctrlcartilage_CPM3.txt", sep="\t")

lrt.2vs1 <- glmLRT(fit, contrast=c(1,0,0,0)) #wrong thinking
top2v1 <- topTags(lrt.2vs1, n=25000)
write.table(top2v1, "DE_VA1loaded_CPM3.txt", sep="\t")

lrt.2vs1 <- glmLRT(fit, contrast=c(0,0,1,0)) #wrong thinking
top2v1 <- topTags(lrt.2vs1, n=25000)
write.table(top2v1, "DE_PA1loaded_CPM3.txt", sep="\t")

lrt.2vs1 <- glmLRT(fit, contrast=c(1,-1,0,0))
top2v1 <- topTags(lrt.2vs1, n=25000)
write.table(top2v1, "DE_VA1loaded_VA1ctrl_CPM3.txt", sep="\t")


lrt.2vs1 <- glmLRT(fit, contrast=c(1,0,-1,0))
top2v1 <- topTags(lrt.2vs1, n=25000)
write.table(top2v1, "DE_VA1loaded_PA1loaded_CPM3.txt", sep="\t")

lrt.2vs1 <- glmLRT(fit, contrast=c(0,0,1,-1))
top2v1 <- topTags(lrt.2vs1, n=25000)
write.table(top2v1, "DE_PA1loaded_PA1ctrl_CPM3.txt", sep="\t")

lrt.2vs1 <- glmLRT(fit, contrast=c(0,1,0,-1))
top2v1 <- topTags(lrt.2vs1, n=25000)
write.table(top2v1, "DE_VA1ctrl_PA1ctrl_CPM3.txt", sep="\t")


