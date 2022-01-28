##### Analysis 

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

####### Input data
ts <- read_csv("output/clustered_time_series.csv")[,-1]
para <- read_csv("output/clustered_parameters.csv")[,-1]

w<-which(is.na(para$cluster))
para[w,"cluster"] <- "no cluster"
w<-which(is.na(ts$cluster))
ts[w,"cluster"] <- "no cluster"

#### Questions
cluster_data <- para %>% dplyr::select(c("strain","rep","drytime","inocl","cluster","valpeak", "timepeak", "auc", "exp_gr", "shoulder_point_v","shoulder_point_t",
                                         "shoulder_point_past_v","shoulder_point_past_t", "mp_h2")) %>% 
  pivot_longer(cols = c(valpeak:mp_h2)) 

cluster_data_summ_all <- cluster_data %>% group_by(cluster, name) %>% summarise(mean = mean(value), sd = sd(value))

g <- ggplot(cluster_data_summ_all, aes(x = cluster, y = mean, fill = name)) + geom_bar(stat="identity", color="black", 
                                                                              position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) + 
  facet_wrap(~name, scales = "free")
ggsave("plots/ALL_clusters_mean_sd.pdf")

g <- ggplot(cluster_data, aes(x = cluster, y = value, fill = name)) + geom_boxplot() + facet_wrap(~name, scales = "free")
ggsave("plots/ALL_clusters_boxplot.pdf")

#### How many of each cluster type by inoculum and drytime? 
para2 <- para
para2$cluster <- factor(para2$cluster, levels = c("double","spike","normal","post_shoulder","wide","no cluster"))
howmany <- para2 %>% group_by(inocl, drytime, cluster,  .drop=FALSE) %>% summarise( n = n()) %>% complete(cluster, fill = list(n = 0))
#ggplot(howmany, aes(x=cluster,y = n, group = interaction(inocl,drytime))) + geom_point(aes(col = interaction(inocl,drytime))) + geom_line(aes(col = interaction(inocl,drytime))) + 
#  facet_wrap(~drytime)

ggplot(howmany, aes(x=cluster,y = n, group = interaction(inocl,drytime))) + geom_bar(stat="identity",position = position_dodge(),aes(fill = factor(inocl)))+ 
  facet_wrap(~drytime, ncol = 1) + scale_fill_discrete("Inoculum")
ggsave("plots/ana_clusters_by_inoc_drytime.pdf") # no double / spike at 3 / 5

### Does normal or double have higher total energy? (AUC)
cluster_data$cluster <- factor(cluster_data$cluster, levels = c("double","spike","normal","post_shoulder","wide","no cluster"))x
cluster_data$inocl <- factor(cluster_data$inocl)
cluster_data$drytime <- factor(cluster_data$drytime)
cluster_data$name <- factor(cluster_data$name)
cluster_data_summ <- cluster_data %>% group_by(drytime, inocl, cluster, name, .drop=FALSE) %>% summarise(mean = mean(value), sd = sd(value))

ggplot(cluster_data_summ %>% filter(name == "auc"), aes(x = cluster, y = mean, fill = factor(inocl))) + 
  geom_bar(stat="identity", color="black", position=position_dodge())+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9)) + 
  facet_wrap(~drytime , scales = "free", ncol = 1) + scale_fill_discrete("Inoculum") + 
  ggtitle("AUC")
ggsave("plots/ana_auc_summary.pdf")

ggplot(cluster_data_summ , aes(x = cluster, y = mean, fill = factor(inocl))) + 
  geom_bar(stat="identity", color="black", position=position_dodge())+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9)) + 
  facet_grid(name~drytime , scales = "free") + scale_fill_discrete("Inoculum") 
ggsave("plots/ana_all_summary.pdf")


ggplot(cluster_data_summ %>% filter(name %in% c("valpeak","auc")), aes(x = cluster, y = mean, fill = factor(inocl))) + 
  geom_bar(stat="identity", color="black", position=position_dodge())+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9)) + 
  facet_grid(name ~ drytime, scales = "free") + scale_fill_discrete("Inoculum") 
ggsave("plots/ana_auc&valpeak_summary.pdf")


