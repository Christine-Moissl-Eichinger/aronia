#################################################################################

####################### R script ALPHA DIVERSITY ################################

####### by Christina Kumpitsch/ AG Moissl-Eichinger ########

#################################################################################
#devtools::install_github(c("david-barnett/microViz","jrnold/ggthemes","jbisanz/qiime2R","gmteunisse/Fantaxtic"))
#BiocManager::install("microbiomeMarker")

library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("plyr"); packageVersion("plyr")
library("cluster"); packageVersion("cluster")
library("igraph"); packageVersion("igraph")
library("markovchain"); packageVersion("markovchain")
library("RColorBrewer")
library("phytools")
library("gridExtra")
library("grid")
library("ggthemes")
library("ggpubr")
library("dplyr")
library("cowplot")
library(microbiomeMarker)


#################################################################################

##### Import QIIME 2 qza files #####

#################################################################################


NBA <- import_qiime2(otu_qza = "verum.qza", taxa_qza = "taxonomy.qza",
  sam_tab = "metadata.txt")


detach("package:microbiomeMarker", unload = TRUE)  #unload microbiome maker again as it interferes with other commands later on

ps <- NBA

#################################################################################

##### ALPHA DIVERSITY #####

#################################################################################

# Installing from Bioconductor
# source("http://www.bioconductor.org/biocLite.R")
#BiocManager::install("MASS")

# Installing from CRAN
#install.packages("sorvi")

# Installing from Github
library(devtools)
#install_github("antagomir/netresponse")

library(microbiome)
library(knitr)

tab <-microbiome::alpha(ps, index = "all")
kable(head(tab))

#Richness

tab_richness <- richness(ps)
kable(head(tab_richness))

plot(tab_richness)

## Testing differences in alpha diversity
# Construct the data
d <- meta(ps)
diversity_shannon <- microbiome::alpha(ps, index = "diversity_shannon")

#Split the values by group
spl <- split(diversity_shannon, d$group)




#################################################################################

###### Preparation to plot alpha diversity ######

#################################################################################

library(microbiome)
library(ggpubr)
library(knitr)
library(dplyr)

# remove all features with less than 1 read
ps1 <- prune_taxa(taxa_sums(ps) > 0, ps)

# calculation of alpha diversities
tab <- microbiome::alpha(ps1, index = "all")
kable(head(tab))

## prepare data for visualization
ps1.meta <- meta(ps1)
kable(head(ps1.meta))

#add diversity table to metadata
ps1.meta$Shannon <- tab$diversity_shannon 
ps1.meta$InverseSimpson <- tab$diversity_inverse_simpson
ps1.meta$richness <- tab$chao1
ps1.meta$evenness <- tab$evenness_pielou

#export table
write.table(ps1.meta, file = 'diversity_verum_overtime.txt', col.names = TRUE,
            row.names = FALSE, sep = "\t")

#compare alpha diversity based on specific metadata
# create a list of pairwise comparisons
polyphenol <- unique(ps1.meta$group) # for 2 variables


##############################################

###### create a boxplot using ggboxplot ######

##############################################


p2.1 <- ggboxplot(ps1.meta, x = "timepoint", y = "Shannon",
               fill = "timepoint", palette = c("#5D3A9B", "#5D3A9B", "#5D3A9B"),
               add = "jitter", 
               xlab = "timepoint", ylab = "Shannon", main = "Shannon index",
               legend = "right")
#  ggtitle("alpha diversity") 

p2.1
  
p2.2 <- ggboxplot(ps1.meta, x = "timepoint", y = "richness",
                  fill = "timepoint", palette = c("#5D3A9B", "#5D3A9B", "#5D3A9B"),
                  add = "jitter",
                  xlab = "timepoint", ylab = "Richness", main = "Richness",
                  legend = "right")
  #  ggtitle("alpha diversity") 
p2.2

p2.3 <- ggboxplot(ps1.meta, x = "timepoint", y = "evenness",
                  fill = "timepoint", palette = c("#5D3A9B", "#5D3A9B", "#5D3A9B"),
                  add = "jitter", 
                  xlab = "timepoint", ylab = "Evenness", main = "Evenness",
                legend = "right")
  #  ggtitle("alpha diversity") 

p2.3

p4 = p2.1 + p2.2 + p2.3

print(p4)


