# from the Aligned ReadsPErGene files, compile a gene counts table 

rm(list=ls())

library(dplyr)

temp = list.files(pattern="*.tab")
library(gtools)
temp <- mixedsort(temp)
myfiles = lapply(temp, read.table)

a2 <- myfiles[1]
a3 <- as.data.frame(a)
b2 <- a3 %>% select(3)

length(temp)

e <- dat %>% select(1)

for (i in 1:length(temp)){
a <- myfiles[i]
b <- as.data.frame(a)
e[i+1] <- b %>% select(3)
}

write.table(e, "Gene_count_compiled.txt", sep="\t")




