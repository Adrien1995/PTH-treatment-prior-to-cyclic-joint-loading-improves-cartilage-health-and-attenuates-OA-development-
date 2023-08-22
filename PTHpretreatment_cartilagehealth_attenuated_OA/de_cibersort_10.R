########Cibersort analysis, sample script

rm(list = ls()) 

library(ggplot2)
#install.packages("ggsignif")
library(ggsignif)
library(dplyr)
library(tidyverse)
#install.packages("phia")
library(phia)
#install.packages("ggpubr")
library(ggpubr)
#install.packages("multcompView")
library(multcompView)


#############Importing results
#m3 <- read.delim("CIBERSORT_1wkLN_relative.txt",header=T)
#m4 <- read.delim("CIBERSORTx_2wkLN_relative.txt",header=T)
m3 <- read.delim("LN1wk.txt",header=T)
m4 <- read.delim("LN2wks.txt",header=T)


m <- rbind(m3,m4)

############### Creating group factors 
x <- m %>% filter(str_detect(Mixture, "\\bVV1_RLN"))
x2 <- as.vector(x[,1])

p <- m %>% filter(str_detect(Mixture, "\\bVV1_LLN"))
p2 <- as.vector(p[,1])

y <- m %>% filter(str_detect(Mixture, "\\bPV1_RLN"))
y2 <- as.vector(y[,1])

q <- m %>% filter(str_detect(Mixture, "\\bPV1_LLN"))
q2 <- as.vector(q[,1])

z <- m %>% filter(str_detect(Mixture, "\\bVA1_RLN"))
z2 <- as.vector(z[,1])

r <- m %>% filter(str_detect(Mixture, "\\bVA1_LLN"))
r2 <- as.vector(r[,1])

v <- m %>% filter(str_detect(Mixture, "\\bPA1_RLN"))
v2 <- as.vector(v[,1])

f <- m %>% filter(str_detect(Mixture, "\\bPA1_LLN"))
f2 <- as.vector(f[,1])


x3 <- m %>% filter(str_detect(Mixture, "\\bVV2_RLN"))
x23 <- as.vector(x3[,1])

p3 <- m %>% filter(str_detect(Mixture, "\\bVV2_LLN"))
p23 <- as.vector(p3[,1])

y3 <- m %>% filter(str_detect(Mixture, "\\bPV2_RLN"))
y23 <- as.vector(y3[,1])

q3 <- m %>% filter(str_detect(Mixture, "\\bPV2_LLN"))
q23 <- as.vector(q3[,1])

z3 <- m %>% filter(str_detect(Mixture, "\\bVA2_RLN"))
z23 <- as.vector(z3[,1])

r3 <- m %>% filter(str_detect(Mixture, "\\bVA2_LLN"))
r23 <- as.vector(r3[,1])

v3 <- m %>% filter(str_detect(Mixture, "\\bPA2_RLN"))
v23 <- as.vector(v3[,1])

m2 <- m %>% mutate(Group = ifelse(Mixture %in% x2, "VV1_RLN", 
                                  ifelse(Mixture %in% p2, "VV1_LLN",
                                         ifelse(Mixture %in% y2, "PV1_RLN", 
                                                ifelse(Mixture %in% q2, "PV1_LLN",
                                                       ifelse(Mixture %in% z2, "VA1_RLN",
                                                              ifelse(Mixture %in% r2, "VA1_LLN",
                                                                     ifelse(Mixture %in% v2, "PA1_RLN",
                                                                            ifelse(Mixture %in% f2, "PA1_LLN",
                                                                                   ifelse(Mixture %in% x23, "VV2_RLN",
                                                                                          ifelse(Mixture %in% p23, "VV2_LLN",
                                                                                                 ifelse(Mixture %in% y23, "PV2_RLN", 
                                                                                                        ifelse(Mixture %in% q23, "PV2_LLN",
                                                                                                               ifelse(Mixture %in% z23, "VA2_RLN",
                                                                                                                      ifelse(Mixture %in% r23, "VA2_LLN",
                                                                                                                             ifelse(Mixture %in% v23, "PA2_RLN", "PA2_LLN"))))))))))))))))



#m2$Load <- factor(m2$Load)
m2$Group <- factor(m2$Group,levels =c("VV1_RLN","VV1_LLN","VV2_RLN","VV2_LLN","VA1_RLN","VA1_LLN","VA2_RLN","VA2_LLN","PV1_RLN","PV1_LLN","PV2_RLN","PV2_LLN","PA1_RLN","PA1_LLN","PA2_RLN","PA2_LLN"))

#m2$T.cells.regulatory..Tregs.
#m2$T.cells.follicular.helper
#m2$T.cells.CD4.memory.resting
#m2$B.cells.memory
#m2$T.cells.CD8

color4 <- c("purple") 
ct4 <- adjustcolor(color4, alpha.f = 0.4) 
ct4a <- adjustcolor(color4, alpha.f = 0.9)
color3 <- c("blue") 
ct3 <- adjustcolor(color3, alpha.f = 0.4) 
ct3a <- adjustcolor(color3, alpha.f = 0.9) 
color2 <- c("red3") 
ct2 <- adjustcolor(color2, alpha.f = 0.4) 
ct2a <- adjustcolor(color2, alpha.f = 0.85) 
color1 <- c("black") 
ct1 <- adjustcolor(color1, alpha.f = 0.4) 
ct1a <- adjustcolor(color1, alpha.f = 0.9) 

ggplot(m2, aes(x=Group, y=B.Cells, fill = Group)) + # fill=name allow to automatically dedicate a color for each group
  geom_violin( ) +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1)) +
  #scale_fill_manual(values=c(ct1,ct1,ct2,ct2,ct3,ct3,ct4,ct4)) +
  
  #geom_jitter(shape=16, position=position_jitter(0.2)) +
  #+ #+ #+
  scale_x_discrete(limits=c("VV1_RLN","VV1_LLN","VV2_RLN","VV2_LLN","VA1_RLN","VA1_LLN","VA2_RLN","VA2_LLN","PV1_RLN","PV1_LLN","PV2_RLN","PV2_LLN","PA1_RLN","PA1_LLN","PA2_RLN","PA2_LLN"))  + 
  #scale_fill_brewer(palette="Blues")
  #theme_classic()
  scale_fill_manual(values=c(ct1,ct1a,ct1,ct1a,ct2,ct2a,ct2,ct2a,ct3,ct3a,ct3,ct3a,ct4,ct4a,ct4,ct4a)) + 
  ylab("Fraction of B cells in mixture") +
  xlab("") +
  #theme(text = element_text(size=18), legend.position = "none")
  theme(text = element_text(size=20),  axis.text.x = element_text(angle = 90), legend.position = "none")




###########Running two factor ANOVA with interaction 
########p <0.05 means the factor is significant, in other words we can say with a 95% confidence interval that the factors 
######## are related 
###########

# model = lm(T.cells.CD4.mempry.resting ~ Group + Load + Group:Load, data = m2)
# plot(interactionMeans(model))
aov4 <- aov(Monocytes ~ Group + Load + Group:Load, data = m2 )
summary.lm(aov4)
summary(aov4)

########TukeyHSD gives you more specific details abt which groups 
sigG <- TukeyHSD(aov4, "Group", ordered=TRUE, conf.level=0.95)
sigG$Group

sigL <- TukeyHSD(aov4, "Load", ordered=TRUE, conf.level=0.95)
sigL$Load

sigGL <- TukeyHSD(aov4, "Group:Load", ordered=TRUE, conf.level=0.95)
sigGL$'Group:Load'

#plot(sig4 , las=1 , col="brown" )  # TukeyHSD post-hoc
#plot(sig5 , las=1 , col="brown" )

######jUst sorting the ones that are significant for viz 

aG <- as.data.frame(sigG$Group) 
dG <- aG["p adj"] 
dG2 <- dG %>% filter(dG < 0.05)
dG2

aL <- as.data.frame(sigL$Load)
dL <- aL["p adj"]
dL2 <- dL %>% filter(dL < 0.05)
dL2

aGL <- as.data.frame(sigGL$'Group:Load')
dGL <- aGL["p adj"]
dGL2 <- dGL %>% filter(dGL < 0.05)
dGL2

##################
# Tukey Test Letters

# group the treatments that are not different each other together.
generate_label_df <- function(TUKEY, variable){
  
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels_lower <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  Tukey.labels_upper <- data.frame(multcompLetters(Tukey.levels,Letters=c("p","q","r","s"))['Letters'])
  
  #Put the labels in the same order as in the boxplot :
  Tukey.labels_upper$factor=rownames(Tukey.labels_upper)
  Tukey.labels_upper=Tukey.labels_upper[order(Tukey.labels_upper$factor) , ]
  Tukey.labels_lower$factor=rownames(Tukey.labels_lower)
  Tukey.labels_lower=Tukey.labels_lower[order(Tukey.labels_lower$factor) , ]
  return(list("upper"=Tukey.labels_upper,"lower"=Tukey.labels_lower))
}
# Apply the function on dataset
LABELS_group <- generate_label_df(sigG, "Group")$upper
LABELS_group <- arrange(LABELS_group,desc(factor))
LABELS_gl <- generate_label_df(sigGL, "Group:Load")$lower
LABELS_gl$load = substr(rownames(LABELS_gl),5,8)
LABELS_gl <- arrange(LABELS_gl,desc(factor))
LABELS_gl <- arrange(LABELS_gl,load)





#############Plotting (Kate's Code)
letters_group <- rep(LABELS_group$Letters,times=2)
if(sigL$Load[4]<0.05){letters_load <- c("X","Y")} else {letters_load  <- c("X","X")}
letters_gl <- LABELS_gl$Letters


pd = position_dodge(width = .75)
labeldat_g = m2 %>% group_by(Load,Group) %>% summarize(Max=max(Monocytes)) %>% add_column(let=letters_group[1:8])
labeldat_l = m2 %>% group_by(Load) %>% summarize(Max=max(Monocytes)) %>% add_column(let=letters_load)
labeldat_gl = labeldat_g = m2 %>% group_by(Load,Group) %>% summarize(Max=max(Monocytes)) %>% add_column(let=letters_gl[1:8])
dat_ggl = labeldat_g
dat_ggl$let = paste(letters_group[1:8],letters_gl[1:8])
max_load = max(labeldat_l$Max)
line.df <- data.frame(X= seq(0.6,1.4, by=0.1),Y= rep(max_load,times=9))


plot0 <- ggplot(m2, aes(x=Load,y=Monocytes,fill=Group)) +
  stat_boxplot(geom="errorbar", position= pd, width=0.2) +
  geom_boxplot(width=0.4, position=pd, outlier.size=3, outlier.shape=18) +
  
  #geom_jitter(color="black", width=0, size=.6, alpha=0.9, position=position_jitter(width=0.05)) +
  geom_point(size=.9,alpha=0.8, position=position_jitterdodge(jitter.width=0)) +
  
  #theme_dark() +
  theme(    
    panel.grid = element_blank(), #Remove all gridlines.
    
    #panel.background = element_blank(), #Remove the plot background.
    panel.background = element_rect(fill="white"),
    plot.background = element_rect(fill="white"),
    
    panel.grid.major.y = element_line(colour="black"),
    #panel.grid.minor.y = element_line(size=(0.1), colour="white"),
    
    axis.ticks.y = element_line(color="black"), #Remove the axis ticks.
    
    #legend.position = "none", #Remove the legend.
    
    axis.title.y = element_text(color="black",face="bold",size=12), #Adjust the axis title text.
    title = element_text(color="black",face="bold",size=12),
    axis.text.x=element_text(angle=45,colour="black",face="bold",size=12),
    #axis.text.x=element_blank(),
    
    legend.background = element_rect(fill="white"),
    legend.text = element_text(color="black", size=12),
    
    #axis.text.x = element_text(color="black", size=25), #Adjust the axis tick label text.
    axis.text.y = element_text(color="black", face="bold",size=12)) +
  
  xlab('Treatment Groups') + #Add axis titles
  ylab(bquote('Monocytes')) + #bquote() can include superscript/subscript 
  
  scale_fill_manual("Legend",values=c("steelblue3","turquoise","springgreen3","palegreen","steelblue3","turquoise","springgreen3","palegreen")) + #Create a custom color scale.
  
  #geom_signif(comparisons = list(c("Ctrl","Load")),map_signif_level=TRUE, color="blue1", na.rm=T) +
  geom_text(data=dat_ggl, aes(label=let,y=Max+0.015), position=pd, show.legend=FALSE) +
  geom_text(aes(label=labeldat_l$let[1],x= 1,y=max_load*1.15), show.legend=FALSE,color="blue1") +
  geom_line(data=line.df, inherit.aes=FALSE, aes(X,Y*1.12),color="blue1") +
  geom_text(aes(label=labeldat_l$let[2],x= 2,y=max_load*1.15), show.legend=FALSE,color="blue1") +
  geom_line(data=line.df, inherit.aes=FALSE, aes(X+1,Y*1.12),color="blue1") +
  
  ggtitle("1wk marrow - Monocytes")
plot(plot0)

png(filename="Plots/1wkmarrow_Monocytes.png")
plot(plot0)
dev.off()



########################################################################

plot1 <- ggplot(m2, aes(x = factor(Mixture, levels = c("VV1_LL","VA1_LL","PV1_LL","PA1_LL","VV1_RL","VA1_RL","PV1_RL","PA1_RL")),
                        y = T.cells.follicular.helper, fill = Mixture)) #Make the base call to ggplot().
#levels = can change the order in which the groups are plotted (for aesthetic purposes).

#pos <- c("VV1_LL","VA1_LL","PV1_LL","PA1_LL","VV1_RL","VA1_RL","PV1_RL","PA1_RL")

pos <- c("VV1_LL","VA1_LL","PV1_LL","PA1_LL","VV1_RL","VA1_RL","PV1_RL","PA1_RL")

plot1 <- plot1 + geom_boxplot(outlier.shape = NA, coef = 0, width = 0.5) + scale_x_discrete(limits=pos) +
  #geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd,width=0.5),color="white",show.legend=FALSE) + 
  
  ggtitle(title) + 
  guides(fill = guide_legend(reverse = TRUE)) +
  
  scale_fill_manual("Legend",values=c("steelblue3","steelblue3","paleturquoise","paleturquoise","springgreen3","springgreen3","palegreen","palegreen"),label= c("PA1L","PA1ctrl","PV1L","PV1ctrl","VA1L","VA1ctrl","VV1L","VV1ctrl")) + #Create a custom color scale.
  
  xlab('Treatment Groups') + #Add axis titles.  
  ylab(bquote(ylabel)) + #bquote() can include superscript/subscript 
  theme_dark() + #Modify the theme.  A preset can also be used.
  theme(    
    panel.grid = element_blank(), #Remove all gridlines.
    
    #panel.background = element_blank(), #Remove the plot background.
    panel.background = element_rect(fill="white"),
    plot.background = element_rect(fill="white"),
    
    panel.grid.major.y = element_line(colour="black"),
    #panel.grid.minor.y = element_line(size=(0.1), colour="white"),
    
    axis.ticks.y = element_line(color="black"), #Remove the axis ticks.
    
    #legend.position = "none", #Remove the legend.
    
    axis.title.y = element_text(color="black",face="bold",size=12), #Adjust the axis title text.
    title = element_text(color="black",face="bold",size=12),
    #axis.text.x=element_text(angle=45,colour="white",face="bold",size=12),
    axis.text.x=element_blank(),
    
    legend.background = element_rect(fill="white"),
    legend.text = element_text(color="black", size=12),
    
    #axis.text.x = element_text(color="black", size=25), #Adjust the axis tick label text.
    axis.text.y = element_text(color="black", face="bold",size=12)) +
  geom_point(color="black")
#ylim(0, 0.5) 
#  geom_signif(comparisons=list(c("PV1_RLN", "VV1_RLN")), color="orange3",annotations="*", y_position = .425, tip_length   = 0, textsize = 12, vjust=0.4) +
# geom_signif(comparisons=list(c("VA1_LLN", "VV1_LLN")), color="orange3",annotations="*", y_position = .375, tip_length   = 0, textsize = 12, vjust=0.4) +
#geom_signif(comparisons=list(c("PA1_LLN", "VV1_LLN")), color="orange3",annotations="*", y_position = .425, tip_length   = 0, textsize = 12, vjust=0.4) 

plot(plot1)
#print(plot1)
png(filename="TcellsCD4naive.png")
plot(plot1)
dev.off()

