##### Analysis of clusters

###******* LOAD UP LIBRARIES AND DATA NEEDED *************#############################################################
## libraries needed
#detach(package:plyr)
library(reshape2) # for data manipulation
library(ggplot2) # for plotting
#library(plyr) # for data manipulation
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

####### Input data
ts <- read_csv("output/clustered_time_series.csv")[,c(-1,-2)]
para <- read_csv("output/clustered_parameters.csv")[,c(-1,-2)]

w<-which(is.na(para$cluster))
para[w,"cluster"] <- "unclustered"
w<-which(is.na(ts$cluster))
ts[w,"cluster"] <- "unclustered"

###### FIGURE Comparison of characteristics of clusters
# if first is not max then second = max, otherwise 1st is max so take next 
para <- para %>% ungroup() %>% mutate(second_peak_h = ifelse((t_m_h_flow > (timepeak + 2)), v_m_h_flow, ifelse(mp_t2 > (timepeak + 4), mp_h2, 0))) 

cluster_data <- para %>% dplyr::select(c("cluster", "inocl","drytime","valpeak", "timepeak", "auc", "exp_gr", "shoulder_point_v","shoulder_point_t",
                                          "shoulder_point_past_v","shoulder_point_past_t", "second_peak_h")) %>% 
  pivot_longer(cols = c(valpeak:second_peak_h)) 

### STORY PIC

my_comparisons <- list( c("post_shoulder", "double"), c("double", "normal"), c("normal", "post_shoulder"),c("normal","spike"),c("spike","wide"),
                        c("double","spike"),c("double","wide"))

cluster_data$cluster <- factor(cluster_data$cluster, levels = c("normal","double","spike","post_shoulder","wide","unclustered"))
g <- ggplot(cluster_data %>% filter(drytime ==0, inocl == 5) %>% filter(name %in% c("auc", "valpeak","exp_gr","second_peak_h")) %>% 
              mutate(across(name, factor, levels=c("auc", "valpeak","exp_gr","second_peak_h"))), aes(x = cluster, y = value)) + 
  geom_boxplot(aes(fill = cluster)) +  facet_wrap(~name,ncol = 4, scales = "free") 
ggsave("plots/inoc_5_story_boxplot.pdf", width = 20, height = 8)

# var_name <- c(
#   auc = "AUC",
#   valpeak = "Max value 1st peak",
#   exp_gr = "Max exp growth rate",
#   second_peak_h = "Max value 2nd peak"
# )

# ggplot(cluster_data %>% filter(drytime ==0, inocl == 5) %>% filter(name %in% c("auc", "valpeak","exp_gr","second_peak_h")) %>% 
#          mutate(across(name, factor, levels=c("auc", "valpeak","exp_gr","second_peak_h"))), aes(x = cluster, y = value)) + 
#   geom_boxplot(aes(fill = cluster), show.legend = FALSE) +  
#   scale_x_discrete("Cluster type", label = c("Normal","Double","Spike","Post-shoulder","Wide","Unclustered")) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
#   scale_y_continuous("Value") +
#   facet_wrap(~name, ncol = 4, scales = "free", labeller = labeller(name = var_name)) +
#   labs(fill = "Cluster")
# ggsave("plots/final/figure3.png", width = 10, height = 5)

g + stat_compare_means(comparisons = my_comparisons)
ggsave("plots/inoc_5_story_boxplot_stats.pdf", width = 20, height = 8)

