rm(list=ls())

#install.packages("BiocManager")
#BiocManager::install("WGCNA")


library(WGCNA)

options(stringsAsFactors = FALSE)

#Read in the data set
Data = read.csv("cartilage_norm.csv")

traitData = read.csv("cartilagetrait.csv");

# Take a quick look at what is in the data set:
dim(Data);
names(Data);

datExpr0 = as.data.frame(t(Data[, -c(1)]));


names(datExpr0) = Data$ID 

rownames(datExpr0) = names(Data)[-c(1)];

gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = hclust(dist(datExpr0), method = "average");

sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut
abline(h = 200000, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 200000, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)



# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

softPower = 6;
adjacency = adjacency(datExpr, power = softPower);

TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

############comparison to trait data 


dim(traitData)
names(traitData)

# remove columns that hold information we do not need.
#allTraits = traitData[, -c(1)];
#allTraits = traitData[, c(1, 2, 3, 4, 5, 6) ];
#allTraits = traitData[, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10) ];
allTraits = traitData


dim(allTraits)
names(allTraits)
# Form a data frame analogous to expression data that will hold the clinical traits.
Samples = rownames(datExpr);
traitRows = match(Samples, allTraits$ID);

datTraits = allTraits[traitRows, -1];
rownames(datTraits) = allTraits[traitRows, 1];
collectGarbage();

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")



#memory.limit(size = 56000 )
# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);


# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")


# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")


#cutheight of 0.25 corresponds to a correlation of 0.75
MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off(

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;

# Save module colors and labels for use in subsequent parts - I
save(MEs, moduleLabels, moduleColors, geneTree, file = "networkConstruction-stepByStep.RData")



#################relating modules to exxternal traits 
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");

#moduleTraitCor = cor(MEs, datTraits, method="kendall", use = "everything");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

###Finding significant genes related to a give trait - in a  given module ################# 

# Define variable weight containing the weight column of datTrait
load = as.data.frame(datTraits$Load);
names(load) = "load"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, load, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(load), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");

############Identifying gens with high gene significance and module membership for a given trait 
module = "brown4"   ###########picking brown4 coz it's highest regulated for loading 
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for load",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

######writing results to a dataframe 

probes1 = names(datExpr)
geneInfo1 = data.frame(geneSymbol = probes1,
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)

# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, load, use = "p")));

for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo1)
  geneInfo1 = data.frame(geneInfo1, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo1) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo1$moduleColor, -abs(geneInfo1$load));
geneInfo1 = geneInfo1[geneOrder, ]

write.csv(geneInfo1, file = "geneInfo.csv")



###########  Comparing cyan module and load 

###########Finding significant genes related to a give trait - in a  given module ######################### 

# Define variable weight containing the weight column of datTrait
load = as.data.frame(datTraits$Load);
names(load) = "load"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, load, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(load), sep="");
names(GSPvalue) = paste("p.GS.", names(load), sep="");

###########Identifying genes with high gene significance and module membership for a given trait 
module = "cyan"   
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for load",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

######writing results to a dataframe 

probes1 = names(datExpr)
geneInfo2 = data.frame(geneSymbol = probes1,
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)

# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, load, use = "p")));

for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo2)
  geneInfo2 = data.frame(geneInfo2, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo2) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo2$moduleColor, -abs(geneInfo2$GS.load));

geneInfo2 = geneInfo2[geneOrder, ]

write.csv(geneInfo2, file = "cyan_load.csv")


########## vv1 with grey module 



# Define variable weight containing the weight column of datTrait
vv1 = as.data.frame(datTraits$VV1);
names(vv1) = "vv1"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, vv1, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(vv1), sep="");
names(GSPvalue) = paste("p.GS.", names(vv1), sep="");

###########Identifying gens with high gene significance and module membership for a given trait 
module = "grey"   
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));

verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for vv1",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

######writing results to a dataframe 

probes1 = names(datExpr)
geneInfo3 = data.frame(geneSymbol = probes1,
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)

# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, vv1, use = "p")));

for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo3)
  geneInfo3 = data.frame(geneInfo3, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo3) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo3$moduleColor, -abs(geneInfo3$GS.vv1));
geneInfo3 = geneInfo3[geneOrder, ]

write.csv(geneInfo3, file = "grey_vv1.csv")




#######################Goenichment analysis

source("https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/GeneAnnotation/installAnRichment.R");
installAnRichment();

install.packages("anRichment", repos = NULL, type = "source");

#################Creating starting data frame 

probes2 = names(datExpr)

geneInfo2 = data.frame(geneSymbol = probes2,
                       moduleColor = moduleColors)


write.csv(geneInfo2, file = "geneInfo_Goenrichment.csv")


options(stringsAsFactors = FALSE);
library("anRichment");

symbol.0 = geneInfo2$geneSymbol;
moduleColor = geneInfo2$moduleColor;
table(moduleColor)

GOcollection = buildGOcollection(organism = "mouse")

split = strsplit(symbol.0, split = " /// ", fixed = TRUE);

symbol = sapply(split, function(x) x[1]);
# Convert symbols to Entrez IDs
entrez = convert2entrez(organism = "mouse", symbol = symbol);
# How many conversions were successful?
table(is.finite(entrez))

GOenrichment = enrichmentAnalysis(
  classLabels = moduleColor, identifiers = entrez,
  refCollection = GOcollection,
  useBackground = "given",
  threshold = 1e-4,
  thresholdType = "Bonferroni",
  getOverlapEntrez = TRUE,
  getOverlapSymbols = TRUE,
  ignoreLabels = "grey");

collectGarbage();


#names(GOenrichment)

#table = GOenrichment$enrichmentTable

#table$overlapGenes = shortenStrings(table$overlapGenes, maxLength = 500,
 #                                   split = "|");

save(MEs, moduleLabels, moduleColors, geneTree, GOenrichment, file = "networkConstruction&GOenrichment-stepByStep.RData")

write.csv(GOenrichment$enrichmentTable, file = "GOenrichment-enrichmentTable.csv",
          row.names = FALSE);


###########Example of how to Extract gene names for pathways with more than 50 overlapping genes 

a <- GOenrichment$dataSetDetails$brown4[[1]]$commonGeneEntrez
a
b = as.character(a)
b


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")

BiocManager::install("org.Mm.eg.db")

library(org.Mm.eg.db)

library(annotate)

#d= getSYMBOL(c('236266','68396','21802'), data='org.Mm.eg.db')
#d
e= getSYMBOL(b, data='org.Mm.eg.db')


##########Visualizing eigengenes network 

######Use the one below this block of code 

# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function
sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")



nSelect = 400
# For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
# Open a graphical window
sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")


###########Plotting eigengenes dendrogram
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MEs, "Eigengene dendrogram and adjacency heatmap", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)

# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MEs, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)


# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)





###########################
library(WGCNA);library(ggplot2); library(reshape); library(igraph); library(RColorBrewer); library(WGCNA)

#eigmat = MEs$eigengenes
#eigmat = consMEsC$data
eigmat = MEs
#eigmat = ConsMEC

#eigmat = eigmat[c(1:2),]
#t <- diag(ConsMEC)
#t<- lowerTri2matrix(ConsMEC, diag=1)

#eigmat = ConsME
#eigmat = eigmat[,paste0("ME", labels2colors(1:13))]
adj = bicor(eigmat)
mds = cmdscale(dist(t(eigmat)), eig = T);   
mds$points[,1] = -1* mds$points[,1]
g1 <- graph.adjacency(as.matrix(adj),mode="undirected",weighted=T,diag=FALSE)
layoutFR <- mds$points
#c = paste0("CD", 1:13, ".",labels2colors(1:13))
#c = paste0("ME",labels2colors(1:14))

c <- colnames(MEs)

colors <- gsub('^.', '', c)

colors <- gsub('^.', '', colors)

#pdf("./results/figures/Manuscript/Fig3B.pdf",height=4,width=4,useDingbats=FALSE)
#pdf("modconnection_dailyloading.pdf",height=4,width=4,useDingbats=FALSE)
tiff('modconnection.tiff', units="in", width=5, height=5, res=600, compression = 'lzw')
edgecolors = numbers2colors(E(g1)$weight, colors = blueWhiteRed(10, gamma=3), signed=T, centered=T, lim=c(-1,1))
plot.igraph(g1, vertex.label = c,
            vertex.label.dist=2.5, 
            vertex.size=18,
            vertex.label.color="black",
            vertex.label.family = "sans",
            vertex.color = colors,
            vertex.label.cex=0.5,
            layout=layoutFR,
            edge.color=edgecolors,
            edge.width=2,asp=1)
labs = seq(1,-1,by=-.25)
str = paste(labs,sep="\n")
#text(-1.25,-1, labels = paste(labs,collapse='\n'),pos = 4,cex = .5)
p = matrix(NA,nrow=9,ncol=4)
p[,1]=-1.35; p[,2]=-1.3
p[,3]=p[,4] = -.6-.6*seq(0,1,by=.12)
for(i in 1:9) {
  lines(x=p[i,1:2],y=p[i,3:4], lwd = 2, col=numbers2colors(as.numeric(labs[i]), colors = blueWhiteRed(10, gamma=3), signed=T, centered=T, lim=c(-1,1)))
  text(x=-1.3,y=p[i,3],labs[i],cex=.4,pos=4)
}
dev.off()



## Part 2) 
## Make Module MDS plot
## --------------------
moduleColors = moduleColors

MEs2 <- t(MEs)

datExpr2 <- t(datExpr)
modColors = data.frame(color= moduleColors,row.names = rownames(datExpr2))

cons_kme = signedKME(datExpr, MEs0)
cons_kme = cons_kme[,-which(colnames(cons_kme)=="kMEgrey")]

maxsize = 20  #plot top 10 hub genes for each module

# [1] "kMEblack"       "kMEblue"        "kMEbrown"       "kMEgreen"       "kMEgreenyellow"
# [6] "kMEmagenta"     "kMEpink"        "kMEpurple"      "kMEred"         "kMEsalmon"     
# [11] "kMEtan"         "kMEturquoise"   "kMEyellow"         



#cons_kme <- t(cons_kme)

gene_idx = order(cons_kme[which(modColors=="coral1"),1], decreasing = T)[1:maxsize]
for(i in c(5,15,16)){
  gene_idx = c(gene_idx, order(cons_kme[,i], decreasing = T)[1:maxsize])
}

hubGenes = character()
hubGenes.kme = numeric()
for(col in c("coral1","darkgreen","orangered4","salmon")) {
  #  col = labels2colors(i)
  
  modgenes = rownames(datExpr2)[which((modColors) == col)]
  kmes = cons_kme[modgenes, paste("kME", col,sep="")]
  top_hubs = modgenes[order(kmes, decreasing=T)[1:maxsize]]
  top_hubs.kme = kmes[order(kmes, decreasing=T)[1:maxsize]]
  hubGenes = c(hubGenes,top_hubs)
  hubGenes.kme = c(hubGenes.kme, top_hubs.kme)
}
gene_idx = match(hubGenes,rownames(datExpr2))

adjMat = bicor(t(datExpr2))
keepgenes = rownames(cons_kme)[gene_idx]
adjMat = adjMat[gene_idx,gene_idx]
topcors=0.65^7
adjMat[adjMat< topcors]=0

#geneSymbols = datProbes$external_gene_id[match(keepgenes, datProbes$ensembl_gene_id)]

geneSymbols = keepgenes
g1 <- graph.adjacency(as.matrix(adjMat),mode="undirected",weighted=T,diag=FALSE)
#

mds = cmdscale(dist(t(adjMat)), eig = T)
layoutFR = mds$points
edgecolors = numbers2colors(E(g1)$weight, colors = redWhiteGreen(100, gamma=4), signed=T, centered=T, lim=c(-1,1))


#pdf("./results/figures/manuscript/Fig...pdf",useDingbats = F, width=12,height=12)
tiff('network.tiff', units="in", width=8, height=8, res=600, compression = 'lzw')
plot.igraph(g1, vertex.label = geneSymbols,
            vertex.label.dist=0.75, edge.width=0.25,
            vertex.size=4, vertex.frame.color="black",
            vertex.label.color="black",
            vertex.color = moduleColors[gene_idx],
            vertex.label.cex=0.75,
            layout=layout.fruchterman.reingold(g1),
            edge.color="green")
dev.off()









##################Prepping files for visualization using Cytospace 

# Recalculate topological overlap if needed
#TOM = TOMsimilarityFromExpr(datExpr, power = 6);
# Read in the annotation file
#annot = read.csv(file = "GeneAnnotation.csv");
# Select modules

#modules = c("brown", "red");
#modules = "black"
modules = "brown4"


# Select module probes
probes = names(datExpr)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];

#modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)];

modGenes = modProbes;

# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];

dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])





