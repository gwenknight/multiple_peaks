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
  dplyr::select(inocl, cluster, strain, drytime, rep) %>% 
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
anti_join(mdd, mdd_cluster) #data for strain 11287 present in mdd, but not in mdd_cluster
arrange(mdd_cluster, strain) #indeed, strain 11287 not present in mdd_cluster

mdd <- mdd %>%
  filter(!(strain %in% c("11287")))

length(unique(mdd$strain))
length(unique(mdd_cluster$strain))
dim(mdd) # 97 118
dim(mdd_cluster) #97  5

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
mdd1 <- dplyr::select(mdd_all, strain, cluster, original_ID : spa_type) #metadata
dim(mdd1) # 97 17
#names(mdd1)

## Data exploration ####

ggplot(mdd1) + 
  geom_bar(aes(x = cluster, fill = success))
ggplot(mdd1) + 
  geom_bar(aes(x = cluster, fill = lineage), position = "fill")
