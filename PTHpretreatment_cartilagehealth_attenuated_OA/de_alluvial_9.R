
rm(list=ls())

#install.packages("ggalluvial")
library(ggalluvial)


df2 = read.csv("combinedDE_updated.csv", header=TRUE)

wk1 <- c("Cartilage_1wk","Cortical_1wk","Cancellous_1wk","LN_1wk")
#df3 <- df2 %>% filter(df2$?..Tissue %in% wk1)
df3 <- df2 %>% filter(df2$Tissue %in% wk1)

wk2 <- c("Cartilage_2wks","Cortical_2wks","Cancellous_2wks","LN_2wks")
df4 <- df2 %>% filter(df2$Tissue %in% wk2)

#e <- c("Cartilage","Cartilage2","Cortical Bone", "Cortical Bone2","Cancellous Bone",  "Cancellous Bone2","Lymph Node", "Lymph Node2")
#colors <- c("purple","blue","red","grey","purple","blue","red","grey","purple","blue","red","grey","purple","blue","red","grey","purple","blue","red","grey","purple","blue","red","grey","purple","blue","red","grey","purple","blue","red","grey")


e <- c("Cartilage","Cartilage2","Cortical Bone", "Cortical Bone2","Cancellous Bone",  "Cancellous Bone2","Lymph Node", "Lymph Node2")
e <- c("VEH-VEH","VEH-ALN","PTH-VEH", "PTH-ALN")

tiff("alluvial1&2wks.tiff", units="in", width=14, height=12, res=300)
ggplot(df2,
       aes(x = Treatment, stratum = Tissue, alluvium = sample,
           y = value,
           fill = Tissue, label = Tissue)) +
  #scale_x_discrete(expand = c(.1, .1)) +
  scale_x_discrete(limits=e) +
  stat_stratum(reverse=FALSE)+
  #stat_stratum(geom = "text", aes(label = Tissue), reverse = FALSE)+
  #geom_flow()+
  geom_flow(reverse = FALSE) +
  #geom_flow(geom = "text", aes(label = ?..Tissue), reverse = FALSE) +
  #geom_flow(reverse = FALSE) +
  
  #geom_stratum(alpha = .5, geom = "text", aes(label = ?..Tissue), reverse = FALSE) +
  geom_stratum(alpha = .5, reverse = FALSE) +
  #geom_text(size = 5) + 
  #geom_text(stat = "stratum", size = 3) +
  #theme(legend.position = "none", text = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  ylab("Number of DEGs") + 
  ggtitle("Differential Expression in all tissues at 1wk & 2wks")
dev.off()


tiff("alluvial1wk.tiff", units="in", width=14, height=12, res=300)
ggplot(df3,
       aes(x = Treatment, stratum = Tissue, alluvium = sample,
           y = value,
           fill = Tissue, label = Tissue)) +
  #scale_x_discrete(expand = c(.1, .1)) +
  #scale_x_discrete(limits=e) +
  #stat_stratum(reverse=FALSE)+
  geom_stratum(reverse=FALSE)+
  #stat_stratum(geom = "text", aes(label = ?..Tissue), reverse = FALSE)+
  #geom_stratum(geom = "text", aes(label = ?..Tissue), reverse = FALSE)+
  
  scale_x_discrete(limits=e) +
  #geom_flow()+
  geom_flow(reverse = FALSE) +
  #geom_flow(geom = "text", aes(label = ?..Tissue), reverse = FALSE) +
  #geom_flow(reverse = FALSE) +
  
  #geom_stratum(alpha = .5, geom = "text", aes(label = ?..Tissue), reverse = FALSE) +
  geom_stratum(alpha = .5, reverse = FALSE) +
  #geom_text(size = 5) + 
  #geom_text(stat = "stratum", size = 3) +
  #theme(legend.position = "none", text = element_text(size = 20)) +
  theme(text = element_text(size = 32)) +
  ylab("Number of DEGs") + 
  ggtitle("Differential Expression in all tissues at 1wk")
dev.off()



tiff("alluvial2wks.tiff", units="in", width=14, height=12, res=300)
ggplot(df4,
       aes(x = Treatment, stratum = Tissue, alluvium = sample,
           y = value,
           fill = Tissue, label = Tissue)) +
  #scale_x_discrete(expand = c(.1, .1)) +
  scale_x_discrete(limits=e) +
  stat_stratum(reverse=FALSE)+
  #stat_stratum(geom = "text", aes(label = Tissue), reverse = FALSE)+
  #geom_flow()+
  geom_flow(reverse = FALSE) +
  #geom_flow(geom = "text", aes(label = ?..Tissue), reverse = FALSE) +
  #geom_flow(reverse = FALSE) +
  
  #geom_stratum(alpha = .5, geom = "text", aes(label = ?..Tissue), reverse = FALSE) +
  geom_stratum(alpha = .5, reverse = FALSE) +
  #geom_text(size = 5) + 
  #geom_text(stat = "stratum", size = 3) +
  #theme(legend.position = "none", text = element_text(size = 20)) +
  theme(text = element_text(size = 20)) +
  ylab("Number of DEGs") + 
  ggtitle("Differential Expression in all tissues at 2wks")
dev.off()
