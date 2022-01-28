#### Clustering

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

theme_set(theme_bw(base_size=12)) # theme setting for plots: black and white (bw) and font size (24)
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(14)

setwd(here::here())

## Load in grofit functions if package no longer working
source("code/grofit_functions.R") # have copied those we use into this R script

## Load in functions for this analysis
source("code/functions_for_heat_curves_additional_double_peak.R")

####### Input data
ddm <- read_csv("output/time_series_all.csv")
ddm <- ddm %>% filter(source == "Macotra")
ddm$rep_no <- as.numeric(sub(".*\\.", "", ddm$rep))

param <- read_csv("output/param_multiple_peaks.csv")
param <- param %>% filter(!is.na(strain))
param$rep_no <- as.numeric(sub(".*\\.", "", param$rep))
param$inoc <- param$inocl

# Some strain removed as not enough reps
strains_in <- as.numeric(unlist(param %>% summarise(unique(strain))))
ddm <- ddm %>% filter(strain %in% strains_in)

##### Cluster
c <- cluster(ddm, param)

write.csv(c$parameters,"output/clustered_parameters.csv")
write.csv(c$ts,"output/clustered_time_series.csv")


#### Cluster on latently corrected data
## Correct for latent period 
# removed lag time
ddm_latent <- left_join(ddm, param[,c("strain", "rep", "lag")], by = c("strain", "rep"))
ddm_latent$value_latent_correct <- ddm_latent$Time - as.numeric(ddm_latent$lag) 
ddm_latent <- ddm_latent %>% filter(value_latent_correct > 0)
# set to zero at lag time
ddm_latent <- ddm_latent %>% group_by(rep, inoc, drytime, strain) %>% mutate(initial_time = min(Time),
                                                                             initial_val_f = ifelse(Time == initial_time, value_J, 0),
                                                                             initial_val = max(initial_val_f),
                                                               value_adj_latent = value - initial_val, 
                                                               time_adj_latent = Time - initial_time)
ggplot(ddm_latent %>% filter(inoc == 5, drytime == 0) %>% ungroup() %>% arrange(desc(cluster)), aes(x=time_adj_latent, y = value_adj_latent, group = interaction(strain,rep))) + 
  geom_line(aes(col = factor(cluster))) + 
  facet_grid(cluster~.)
ggsave("plots/inoc_5_clustered_rows_correct_latent.pdf")



