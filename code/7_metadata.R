#### Metadata analysis

###******* LOAD UP LIBRARIES AND DATA NEEDED *************#############################################################
## libraries needed
library(reshape2) # for data manipulation
library(ggplot2) # for plotting
library(plyr) # for data manipulation
library(data.table)
library(reshape2)
library(magrittr)
library(dplyr)
library(MASS)
library(tidyverse)
library(RColorBrewer)
library(here)
library(patchwork)
library(ggpubr)
library(viridis)
library(gridExtra)
library(mlogit)
library(gmodels)

theme_set(theme_bw(base_size=12)) # theme setting for plots: black and white (bw) and font size (24)
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(14)

setwd(here::here())

## Load in grofit functions if package no longer working
source("code/grofit_functions.R") # have copied those we use into this R script

## Load in functions for this analysis
source("code/functions_for_heat_curves_additional_double_peak.R")

#### Read in data
## Read in clusterdata
c <- c()
c$parameters <- read_csv("output/clustered_parameters.csv")

## Select clusterdata for each strain. Decided to base metadata analysis of drytime 0, inoc 5 data. First check if all replicates per strain show the same cluster type
# Then deduplicate to keep that one cluster type per strain

mdd_cluster <- c$parameters %>% 
  dplyr::select(inocl, cluster, strain, drytime, rep, width_peak) %>% #or without width_peak
  filter(drytime == 0) %>%
  filter(inocl == 5) %>% 
  filter(!(strain %in% c("Newman", "RWW12", "SA3297", "SA2704", "RWW146", "SAC042W", "Mu50", "M116")))
#length(unique(mdd_cluster$strain)) #97 

mdd_cluster <- mdd_cluster %>%
  dplyr::mutate(cluster = replace_na(cluster, "unclustered"))
#unique(mdd_cluster$cluster)

mdd_cluster %>%
  group_by(strain) %>%
  summarise(u=unique(cluster))%>% 
  group_by(strain) %>%
  summarise(n = n()) %>%
  filter(n > 1)
#All strains have only one cluster for multiple replicates of drytime 0, inoc 5 -> replicates can be filtered

mdd_cluster <- mdd_cluster[!duplicated(mdd_cluster$strain),]
dim(mdd_cluster) # 97 x 5

## Read in metadata
mdd <- read_csv2("data/MACOTRA_100metadata.csv")
names(mdd)
colnames(mdd)[colnames(mdd) == 'spa...16'] <- 'spa_type'
colnames(mdd)[colnames(mdd) == 'spa...112'] <- 'spa'
mdd$success[mdd$success == "successful"] <- "Successful"

## Combine clusterdata with metadata
mdd_orig <- mdd
mdd <- as.data.frame(mdd)
mdd_cluster <- as.data.frame(mdd_cluster)

mdd$strain <- as.character(mdd$strain)
#anti_join(mdd, mdd_cluster) #data for strain 11287 present in mdd, but not in mdd_cluster
#arrange(mdd_cluster, strain) #indeed, strain 11287 not present in mdd_cluster

mdd <- mdd %>%
  filter(!(strain %in% c("11287")))

# length(unique(mdd$strain))
# length(unique(mdd_cluster$strain))
# dim(mdd) # 97 118
# dim(mdd_cluster) #97  5

mdd_all <- merge(mdd_cluster, mdd, by.x = "strain")
dim(mdd_all)
names(mdd_all)

##Subset metadata dataframes

mdd_amr <- dplyr::select(mdd_all, strain, cluster, original_ID : dfrK) #metadata + amr gene data
dim(mdd_amr) # 97 37
#names(mdd_amr)
mdd_vir <- dplyr::select(mdd_all, strain, cluster, original_ID : spa_type, adsA : vWbp) #metadata + virulence gene data
dim(mdd_vir) # 97 99
#names(mdd_vir)
#mdd1 <- dplyr::select(mdd_all, strain, cluster, original_ID : spa_type) #metadata
mdd1 <- dplyr::select(mdd_all, strain, cluster, width_peak : spa_type) #metadata
dim(mdd1) # 97 17
#names(mdd1) #cluster, country, year, success, carriage, site, lineage, mec, pvl, spa_type #width_peak

## Data description ####

mdd1%>%
  count(cluster)

# cluster  n
# 1        double  5
# 2        normal 56
# 3 post_shoulder 20
# 4         spike  6
# 5   unclustered  2
# 6          wide  8
# sum       total 97

mdd1%>%
  group_by(lineage, cluster) %>%
    count(cluster, lineage) %>% group_by(lineage) %>% dplyr::summarise(n=n())

# datatable <- mdd1%>%  group_by(lineage, cluster) %>%  count(cluster, lineage)
# datatable <- as.matrix(datatable)

table_lineage <- table(mdd1$lineage, mdd1$cluster)
# png("plots/table_cluster_lineage.png", height=4, width=7, units = "in", res = 72)
# grid.table(table_lineage)
# dev.off()

## Data exploration ####
## Main question: are any metabolic clusters correlated with infection and/or success? Correct for group sizes
## If so, are there also correlations with infection site, AMR/virulence genes and/or genetic lineage?

ggplot(mdd1) + 
  geom_bar(aes(x = cluster, fill = country), position = "fill") +
  scale_x_discrete("Cluster type", label = c("Double (n=5)","Normal (n=56)","Post-shoulder (n=20)","Spike (n=6)","Unclustered (n=2)","Wide (n=8)")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
#ggsave("plots/metadata_success.png", width = 10, height = 6)
# more unsuccessful in non-normal clusters?
ggplot(mdd1) + 
  geom_bar(aes(x = cluster, fill = carriage), position = "fill") +
  scale_x_discrete("Cluster type", label = c("Double (n=5)","Normal (n=56)","Post-shoulder (n=20)","Spike (n=6)","Unclustered (n=2)","Wide (n=8)")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
#ggsave("plots/metadata_carriage.png", width = 10, height = 6)
# more infection in non-normal clusters?
ggplot(mdd1) + 
  geom_bar(aes(x = cluster, fill = site), position = "fill") +
  scale_x_discrete("Cluster type", label = c("Double (n=5)","Normal (n=56)","Post-shoulder (n=20)","Spike (n=6)","Unclustered (n=2)","Wide (n=8)")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
#ggsave("plots/metadata_site.png", width = 10, height = 6)
# correlation with infection/carriage? focus on that first, then try to minimize categories in site
ggplot(mdd1) + 
  geom_bar(aes(x = cluster, fill = pvl), position = "fill") +
  scale_x_discrete("Cluster type", label = c("Double (n=5)","Normal (n=56)","Post-shoulder (n=20)","Spike (n=6)","Unclustered (n=2)","Wide (n=8)")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
#ggsave("plots/metadata_pvl.png", width = 10, height = 6)
# pvl positive in double and spike, related to infection/carriage?
ggplot(mdd1) + 
  geom_bar(aes(x = cluster, fill = lineage), position = "fill") +
  scale_x_discrete("Cluster type", label = c("Double (n=5)","Normal (n=56)","Post-shoulder (n=20)","Spike (n=6)","Unclustered (n=2)","Wide (n=8)")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
#ggsave("plots/metadata_lineage.png", width = 10, height = 6)
# correlation with CC30, CC22 and CC5?

##Statistical model ####
# Only keep lineages CC1, CC22, CC30, CC% and CC8
mdd2 <- mdd1 %>%
  filter(lineage %in% c("CC1", "CC22", "CC30", "CC5", "CC8"))

#Group cluster types Double, Spike and Post-shoulder together in an Odd cluster type, discard Unclustered ones
mdd3 <- mdd2 %>%
  filter(!(cluster %in% c("unclustered")))

mdd3$cluster[mdd3$cluster %in% c("double", "spike", "post_shoulder")] <- "odd" 
unique(mdd3$cluster)

## Paramatric test: multinomial logistic regression > not enough data for this
# modeldata <- mlogit.data(mdd2, choice = "cluster", shape = "wide")
# 
# # m0 <- mlogit(cluster ~ 1, data = mdd2, reflevel = "normal") #baseline model
# # summary(m0)
# m1 <- mlogit(cluster ~ 1 | lineage, data = modeldata, reflevel = "normal")
# summary(m1)
# # model with success, carriage, site or pvl does not explain variability in the data
# # model with lineage does explain variability in the data > attributed to CC30 in cluster wide

## Chi-square test

fisher.test(mdd1$cluster, mdd1$site, simulate.p.value=TRUE) #large dataset, simulated p-value for lineage and site
#fisher.test(mdd2$cluster, mdd2$lineage, simulate.p.value=TRUE)

# CrossTable(mdd1$lineage, mdd1$cluster)
# chisq.test(mdd1$lineage, mdd1$cluster)
# fisher.test(mdd1$cluster, mdd1$lineage, simulate.p.value=TRUE)
# CrossTable(mdd2$lineage, mdd2$cluster)
# chisq.test(mdd2$lineage, mdd2$cluster)
# fisher.test(mdd2$cluster, mdd2$country, simulate.p.value=TRUE)

table(mdd2$lineage, mdd2$cluster)
table(mdd3$lineage, mdd3$cluster)

#CrossTable(mdd3$cluster, mdd3$lineage)
#chisq.test(mdd3$cluster, mdd3$success)
fisher.test(mdd3$cluster, mdd3$success, simulate.p.value=FALSE)

#chisq.test(mdd3$cluster, mdd3$carriage)
fisher.test(mdd3$cluster, mdd3$carriage, simulate.p.value=FALSE)

chisq.test(mdd3$cluster, mdd3$lineage)
fisher.test(mdd3$cluster, mdd3$lineage, simulate.p.value=TRUE)

# library(rstatix)
# fisher_test(table(mdd2$lineage, mdd2$cluster), simulate.p.value=TRUE)
# pairwise_fisher_test(as.matrix(table(mdd2$lineage, mdd2$cluster)), p.adjust.method = "fdr")
# how to convert table to two-dimensional contingency table? maybe not working as 2xc comparison, too many categories?

ggplot(mdd3 , aes(cluster, width_peak)) +
  geom_boxplot(aes(colour = lineage)) +
  geom_point(aes(colour = lineage)) +
  facet_wrap(~lineage)
ggsave("plots/metadata_width_peakCC30.png", width = 10, height = 6)

#Not for further analysis
# ggplot(mdd1) + 
#   geom_bar(aes(x = cluster, fill = country), position = "fill") +
#   scale_x_discrete("Cluster type", label = c("Double (n=5)","Normal (n=56)","Post-shoulder (n=20)","Spike (n=6)","Unclustered (n=2)","Wide (n=8)"))
# # biased as mostly Dutch strains included
# mdd1$year <- as.factor(mdd1$year)
# ggplot(mdd1) + 
#   geom_bar(aes(x = cluster, fill = year), position = "fill") +
#   scale_x_discrete("Cluster type", label = c("Double (n=5)","Normal (n=56)","Post-shoulder (n=20)","Spike (n=6)","Unclustered (n=2)","Wide (n=8)"))
# #biased selection + too many variation in year to see effect
# ggplot(mdd1) + 
#   geom_bar(aes(x = cluster, fill = mec)) +
#   scale_x_discrete("Cluster type", label = c("Double (n=5)","Normal (n=56)","Post-shoulder (n=20)","Spike (n=6)","Unclustered (n=2)","Wide (n=8)"))
# # no correlation, all MRSA should have mecA (or mecC)
# ggplot(mdd1) + 
#   geom_bar(aes(x = cluster, fill = spa_type)) +
#   scale_x_discrete("Cluster type", label = c("Double (n=5)","Normal (n=56)","Post-shoulder (n=20)","Spike (n=6)","Unclustered (n=2)","Wide (n=8)"))
# # not enough info to analyse