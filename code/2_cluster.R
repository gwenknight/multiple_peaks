#### Clustering
# Clusters data and makes output of strain distribution in each cluster

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


########## ************************************************************************************************ #########################
##### Clustering #######
c <- cluster(ddm, param) #ddm 97 strains but param 105 strains
########## ************************************************************************************************ #########################

### Store output of clustering
write.csv(c$parameters,"output/clustered_parameters.csv")
write.csv(c$ts,"output/clustered_time_series.csv")

## Read in if doing later
c <- c()
c$parameters <- read_csv("output/clustered_parameters.csv")
c$ts <- read.csv("output/clustered_time_series.csv")


###### FIGURE: e.g.s of clusters
cluster_types <- unique(c$ts$cluster)
for(i in cluster_types){
  ggplot(c$ts %>% filter(cluster == i, inoc == 5, drytime == 0), aes(x=Time, y = value_J, group = interaction(strain, rep))) + 
    geom_line(aes(col = factor(rep_no))) + facet_wrap(~strain) + scale_color_discrete("Replicate")
  ggsave(paste0("plots/ALL_10_5_",i,".pdf"))
}

## Normal ## Strain 11001 - take first one
## Double ## Strain 1040 - clear double peaks
## Post shoulder ## Strain 11073 - rep 2
## Spike ## Strain 11210 - all 
## Wide ## 11106
# Example strains 
eg_cluster_strains = c$ts %>% filter(strain %in% c("11001", "11040", "11073", "11210","11106"), inoc == 5, drytime == 0)
ggplot(eg_cluster_strains, aes(x=Time, y = value_J, group = interaction(strain, rep))) +
  geom_line(aes(col = factor(rep_no))) + facet_wrap(~strain) + scale_color_discrete("Replicate") +
  scale_x_continuous(lim = c(0,20)) +
  geom_text(data = eg_cluster_strains, aes(x = 5, label = cluster), y = Inf, vjust = 2)
ggsave("plots/Example_each_cluster_strain_type.pdf")
# ggplot(eg_cluster_strains, aes(x=Time, y = value_J, group = interaction(strain, rep))) +
#   geom_line(aes(col = factor(rep_no))) + facet_wrap(~strain) + scale_color_discrete("Replicate") +
#   scale_x_continuous("Time (h)", minor_breaks = seq(0, 20, 5), lim = c(0,20)) +
#   scale_y_continuous(expression(paste("Heatflow (mW)")), limits = c(-0.005, 0.1), breaks = c(0, 0.02, 0.04, 0.06, 0.08, 0.10)) +
#   geom_text(data = eg_cluster_strains, aes(x = 5, label = cluster), y = Inf, vjust = 2, hjust = 0.35)
# ggsave("plots/final/figure1.png", width = 8, height = 6)
# ggsave("plots/final/figure1.tiff", width = 8, height = 6, dpi = 600)

#### TABLE: number of strains in each cluster
### CHECK: that unique cluster per strain:
#c$parameters %>% dplyr::select(inocl, cluster, strain, drytime) %>% group_by(inocl, strain, drytime) %>% summarise(u = unique(cluster)) %>% group_by(inocl, strain, drytime) %>% summarise(n = n()) %>% filter( n > 1) 
## None have more than one cluster per strain and inoculum in each drytime

#length(unique(c$parameters$strain)) ### CHECK: 97 strains. Val?rie: 105 strains


# gf1 <- c$parameters %>% filter(drytime == 0) %>% dplyr::select(inocl, cluster, strain) %>% 
#   group_by(inocl, strain) %>% summarise(u = unique(cluster)) %>% group_by(inocl, u) %>% dplyr::summarise(n=n())
gf1 <- c$parameters %>% 
  filter(drytime == 0) %>% 
  filter(!(strain %in% c("Newman", "RWW12", "SA3297", "SA2704", "RWW146", "SAC042W", "Mu50", "M116"))) %>%
  dplyr::select(inocl, cluster, strain) %>% 
  group_by(inocl, strain) %>% summarise(u = unique(cluster)) %>% group_by(inocl, u) %>% dplyr::summarise(n=n())
# length(unique(gf1$strain)) # CHECK now 97 strains
# Add in unclustered name
w<-which(gf1$u == "")
gf1[w,"u"] <- "unclustered"

table_cluster_distribution <- gf1 %>% pivot_wider(names_from = u, values_from = n, values_fill = 0) %>% arrange(desc(inocl))
rowSums(table_cluster_distribution) # CHECK: 97 strains + inocl column = 100/101/102
table_cluster_distribution <- rename(table_cluster_distribution, Inoculum = inocl, Double = double, 
                                     Normal = normal, Spike = spike, Wide = wide, Postshoulder = post_shoulder, Unclustered = "NA")
table_cluster_distribution <- table_cluster_distribution[,c("Inoculum","Normal","Double", "Spike","Postshoulder","Wide", "Unclustered")]
rownames(table_cluster_distribution) <- NULL

pdf("plots/table_cluster_distribution.pdf", height=11, width=8.5)
grid.table(table_cluster_distribution)
dev.off()
# png("plots/final/table_cluster_distribution.png", height=2, width=7, units = "in", res = 72)
# grid.table(table_cluster_distribution)
# dev.off()

### Table with how they change across inoculum 
table_changes <- c$parameters %>% ungroup() %>% dplyr::select(strain, rep, drytime, inocl, cluster) %>% 
  pivot_wider(names_from = inocl, values_from = cluster)
write.csv(table_changes, "output/table_changes_cluster_over_inoc.csv")

colnames(table_changes) <- c("strain", "rep", "drytime", "five", "four", "three")
t_changes <- table_changes %>% mutate(fivetofour = ifelse(five == four,0,1),
                         fourtothree = ifelse(four == three,0,1),
                         difference = fivetofour + fourtothree) %>% filter(difference > 0)
write.csv(t_changes, "output/table_changes_cluster_over_inoc_with_differences.csv")

tt_changes <- table_changes %>% mutate(fivetofour = ifelse(five == four,0,1),
                                      fourtothree = ifelse(four == three,0,1),
                                      fivetothree = ifelse(five == three,0,1),
                                      difference = fivetofour + fourtothree + fivetothree) #%>% filter(difference > 0)
write.csv(tt_changes, "output/tt_changes_cluster_over_inoc_with_more_differences.csv")
#write.csv2(tt_changes, "output/table_changes_cluster_over_inoc_with_more_differences.csv")

td_changes <- table_changes %>% mutate(fivetofour = ifelse(five == four,0,1),
                                       fourtothree = ifelse(four == three,0,1),
                                       fivetothree = ifelse(five == three,0,1),
                                       difference = fivetofour + fourtothree + fivetothree) %>% filter(difference > 0)
dim(td_changes) #259 x 10
td_changes <- td_changes %>%  filter(drytime == 0) # 129 x 10
td_changes <- td_changes[!td_changes$five == "",] # 124 x 10
td_changes <- td_changes[!td_changes$three == "",] # 118 x 10

write.csv(td_changes, "output/td_changes_cluster_over_inoc_with_more_differences.csv")

unique(td_changes$five)
unique(td_changes$three)
unique(td_changes[, c("five", "three")]) # 15 unique combinations of clusters between five and three

count_td_changes <- dplyr::count_(td_changes, vars = c('five','three'))
write.csv(count_td_changes, "output/count_td_changes.csv")

#### FIGURE: Patterns across inoc
gf1$u <- factor(gf1$u, levels = c("normal","double","spike","post_shoulder","wide","unclustered"))
ggplot(gf1, aes(u, inocl, fill= n)) + geom_tile() + scale_fill_viridis(discrete=FALSE, "Number of\nstrains") + 
  scale_x_discrete("Cluster type", labels = c("Normal","Double","Spike","Post shoulder","Wide","Unclustered")) + 
  scale_y_continuous("Inoculum size") 
ggsave("plots/heatmap_cluster_by_inoculum.png")

ggplot(gf1, aes(x= inocl, y = n, fill= u)) + 
  geom_bar(position="dodge", stat="identity") + 
  scale_fill_discrete("Cluster", labels = c("Normal","Double","Spike","Post shoulder","Wide","Unclustered")) +
  scale_y_continuous("Number of strains") + 
  scale_x_continuous("Inoculum") + 
  geom_text(aes(label=n),position=position_dodge(width=0.9), vjust=-0.25) 
ggsave("plots/final/figure2.png",width = 8, height = 6, dpi = 600)


# ggplot(gf1, aes(u, inocl, fill= n)) + geom_tile() + scale_fill_viridis(discrete=FALSE, "Number of\nstrains") +
#   scale_x_discrete("Cluster type", label = c("Normal","Double","Spike","Post-shoulder","Wide","Unclustered")) +
#   scale_y_continuous("Inoculum size")
# ggsave("plots/final/figure2.png", width = 8, height = 4)
# ggsave("plots/final/figure2.tiff", width = 8, height = 6, dpi = 600)
