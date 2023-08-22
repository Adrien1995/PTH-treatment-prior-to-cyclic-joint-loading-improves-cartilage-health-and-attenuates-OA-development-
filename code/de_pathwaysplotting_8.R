
rm(list=ls())

library(ggplot2)
library(dplyr)
library(forcats)
library(stringr)



vv1cart <- read.csv("gseresults.csv")

#######grabbing top downregulatd pathways
vv1cart_down <- vv1cart[order(vv1cart$NES),]
gene_count<- vv1cart_down %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) + 1)
vv1cart_down<- left_join(vv1cart_down, gene_count, by = "ID") %>% mutate(GeneRatio = -count/setSize)
#vv1cart_down <- vv1cart_down[order(-vv1cart_down$GeneRatio),]


########grabbing top upregulated pathways
vv1cart_up <- vv1cart[order(-vv1cart$NES),]
gene_count<- vv1cart_up %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) + 1)


## merge with the original dataframe
vv1cart_up<- left_join(vv1cart_up, gene_count, by = "ID") %>% mutate(GeneRatio = count/setSize)
#vv1cart_up <- vv1cart_up[order(-vv1cart_up$GeneRatio),]

#change number of pathways
vv1cart_up <- vv1cart_up[c(1:5),]
vv1cart_down <- vv1cart_down[c(1:5),]

vv1cart_tot <- rbind(vv1cart_up, vv1cart_down)
vv1cart_tot <- vv1cart_tot[order(vv1cart_tot$NES),]


#doptplot option
#tiff("vv1canc_new_110922.tiff", units="in", width=14, height=12, res=300)
ggplot(vv1cart_tot, aes(x = GeneRatio, y = fct_reorder(Description, GeneRatio))) + 
  geom_point(aes(size = abs(GeneRatio), color = p.adjust)) +
  #theme(axis.text = element_text(size = 30)) +
  theme_bw(base_size = 20) +
  #theme(axis.text.x=element_text(size=30))
  scale_colour_gradient(limits=c(0, 0.10), low="red") +
  ylab(NULL) +
  xlim(-1,1)+
  ggtitle("GO BP")
#dev.off()


##########barplot option

ggplot(data=vv1cart_tot, aes(x=fct_reorder(Description,NES), y=GeneRatio))+
  geom_bar(position="stack", stat="identity")+
  labs(title="GOBiologicalProcess Comparison")+
  #theme(axis.text.x=element_text(angle=90,hjust=1)) +
  coord_flip()+
  labs(y="-log10(Adjusted.P.value)")+
  scale_fill_gradient(low="grey",high="brown")
#scale_fill_gradient(low="grey",high="black")

#p + facet_grid(.~type)
