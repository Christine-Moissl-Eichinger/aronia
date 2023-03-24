#boxplots
library(ggplot2)
library(dplyr)
data <- read.table("clr_data_Holasek_genus_for_boxplot_v3.txt", header=TRUE, sep="", row.names="ID")
data$genus <- factor(data$Genus,
                     labels = c("01_Ruminococcus","02_UCG_005", "03_Monoglobus", "04_Odoribacter", "05_[Eubacterium]_coprostanoligenes_group", 
                                "06_Lachnospiraceae_UCG_004"))
#p10 <- ggplot(data, aes(x = RSV, y = value, fill=tp)) +
  #geom_boxplot() + theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))
#p10

#Tabelle anpassen
library(reshape)

#define group colors
#group.colors <- c(xxx = "red", 02_mother_pre = "blue", 03_mother_post ="palegreen3")

p10 <- ggplot(data, aes(x = Genus, y = value, fill=type)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,outlier.size=0.5, notch=FALSE) + 
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, vjust=.5, hjust=1))+
  scale_fill_brewer(palette = "BuPu")
p10
p10 +labs(x="", y="clr abundance", title="Genera")
p10


