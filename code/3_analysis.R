##### Analysis 

###******* LOAD UP LIBRARIES AND DATA NEEDED *************#############################################################
## libraries needed
detach(package:plyr)
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
ts <- read_csv("output/clustered_time_series_farapart.csv")[,-1]
para <- read_csv("output/clustered_parameters_farapart.csv")[,-1]

w<-which(is.na(para$cluster))
para[w,"cluster"] <- "no cluster"
w<-which(is.na(ts$cluster))
ts[w,"cluster"] <- "no cluster"

#### Questions
### Second peaks? 
#OLD: if more than 5 from end of exp then ok - check up to third minor. mp_t1 = first peak ##ifelse(mp_h2 > 0, ifelse(mp_t2 > (timepeak + 5), mp_h2, ifelse(mp_h3>0,ifelse(mp_h3 > (timepeak + 5), mp_h3,0),0)),0)) 
para <- para %>% ungroup() %>% mutate(second_peak_h = ifelse((t_m_h_flow > (timepeak + 2)), v_m_h_flow, ifelse(mp_t2 > (timepeak + 4), mp_h2, 0))) # if first is no max then second = max, otherwise 1st is max so take next 


cluster_data <- para %>% dplyr::select(c("strain","rep","drytime","inocl","cluster","valpeak", "timepeak", "auc", "exp_gr", "shoulder_point_v","shoulder_point_t",
                                         "shoulder_point_past_v","shoulder_point_past_t", "second_peak_h")) %>% 
  pivot_longer(cols = c(valpeak:second_peak_h)) 

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
para2$cluster <- factor(para2$cluster, levels = c("normal","double","spike","post_shoulder","wide","no cluster"))
howmany <- para2 %>% group_by(inocl, drytime, cluster,  .drop=FALSE) %>% summarise( n = n()) %>% complete(cluster, fill = list(n = 0))
#ggplot(howmany, aes(x=cluster,y = n, group = interaction(inocl,drytime))) + geom_point(aes(col = interaction(inocl,drytime))) + geom_line(aes(col = interaction(inocl,drytime))) + 
#  facet_wrap(~drytime)

ggplot(howmany, aes(x=cluster,y = n, group = interaction(inocl,drytime))) + geom_bar(stat="identity",position = position_dodge(),aes(fill = factor(inocl)))+ 
  facet_wrap(~drytime, ncol = 1) + scale_fill_discrete("Inoculum")
ggsave("plots/ana_clusters_by_inoc_drytime.pdf") # no double / spike at 3 / 5

### Does normal or double have higher total energy? (AUC)
cluster_data$cluster <- factor(cluster_data$cluster, levels = c("normal","double","spike","post_shoulder","wide","no cluster"))
cluster_data$inocl <- factor(cluster_data$inocl)
cluster_data$drytime <- factor(cluster_data$drytime)
cluster_data$name <- factor(cluster_data$name)
cluster_data_summ <- cluster_data %>% group_by(drytime, inocl, cluster, name, .drop=FALSE) %>% summarise(mean = mean(value), sd = sd(value))



### STORY PIC
g <- ggplot(cluster_data %>% filter(drytime ==0, inocl == 5) %>% filter(name %in% c("auc", "valpeak","exp_gr","second_peak_h")) %>% 
              mutate(across(name, factor, levels=c("auc", "valpeak","exp_gr","second_peak_h"))), aes(x = cluster, y = value)) + 
  geom_boxplot(aes(fill = cluster)) +  facet_wrap(~name,ncol = 4, scales = "free") 
ggsave("plots/inoc_5_story_boxplot.pdf", width = 20, height = 8)
g + stat_compare_means(comparisons = my_comparisons)
ggsave("plots/inoc_5_story_boxplot_stats.pdf", width = 20, height = 8)
















###### OLD ANALYSIS

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


ggplot(cluster_data_summ %>% filter(name %in% c("valpeak","auc","timepeak")), aes(x = cluster, y = mean, fill = factor(inocl))) + 
  geom_bar(stat="identity", color="black", position=position_dodge())+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9)) + 
  facet_grid(name ~ drytime, scales = "free") + scale_fill_discrete("Inoculum") 
ggsave("plots/ana_auc&val&timepeak_summary.pdf")

ggplot(cluster_data_summ %>% filter(name %in% c("valpeak","auc","timepeak"), cluster %in% c("normal","spike","double")), aes(x = cluster, y = mean, fill = factor(inocl))) + 
  geom_bar(stat="identity", color="black", position=position_dodge())+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9)) + 
  facet_grid(name ~ drytime, scales = "free") + scale_fill_discrete("Inoculum") 
ggsave("plots/ana_auc&val&timepeak_norm_double.pdf")

my_comparisons <- list( c("post_shoulder", "double"), c("double", "normal"), c("normal", "post_shoulder"),c("normal","spike"),c("spike","wide"),
                        c("double","spike"),c("double","wide"))

ggplot(cluster_data_summ , aes(x = cluster, y = mean, fill = factor(inocl))) + 
  geom_bar(stat="identity", color="black", position=position_dodge())+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9)) + 
  facet_grid(name~drytime , scales = "free") + scale_fill_discrete("Inoculum") + 
  stat_compare_means(comparisons = my_comparisons)

#### How does AUC look vs heat of first / second peak? 
g1 <- ggplot(para %>% filter(drytime ==0, inocl == 5), aes(x=auc, y = valpeak, group = cluster)) + geom_point(aes(col = cluster)) + geom_smooth(method='lm',formula= y~x,se=FALSE, size=2, aes(col = cluster)) + 
  ggtitle("AUC vs height of first peak")
g2 <- ggplot(para %>% filter(drytime ==0, inocl == 5), aes(x=auc, y = mp_h2, group = cluster)) + geom_point(aes(col = cluster)) + geom_smooth(method='lm',formula= y~x,se=FALSE, size=2, aes(col = cluster)) + 
  ggtitle("AUC vs height of second peak")
g1 + g2
ggsave("plots/inoc_5_AUC_vs_heightfp.pdf")

ggplot(para %>% filter(drytime ==0, inocl == 5), aes(x=auc, y = mp_h2, group = cluster)) + geom_point(aes(col = cluster)) + geom_smooth(method='lm',formula= y~x,se=FALSE, size=2, aes(col = cluster)) + 
  ggtitle("AUC vs height of first peak")

ggplot(cluster_data %>% filter(drytime ==0, inocl == 5, name %in% c("mp_h2","auc","valpeak","exp_gr")) , aes(x = cluster, y = value, fill = factor(cluster))) + 
  geom_boxplot() + 
  facet_wrap(~name, scales = "free") + 
  stat_compare_means(comparisons = my_comparisons)
ggsave("plots/inoc_5_auc_vs_height_stats.pdf", width = 15, height = 6)

# which normal strains have second peaks? 
n_2nd <- para %>% filter(mp_h2 > 0, cluster == "normal", drytime ==0, inocl == 5)
ts$strain <- as.numeric(ts$strain)
ddm_n2nd<- left_join(ts, n_2nd[,c("strain", "rep", "mp_h2","mp_t2")], by = c("strain", "rep"))
ggplot(ddm_n2nd %>% filter(drytime ==0, !is.na(mp_h2), inoc == 5, mp_h2 > 0, cluster == "normal"), aes(x=Time, y = value_J, group = interaction(strain,rep))) + 
  geom_line(aes(col = factor(cluster))) + facet_wrap(~strain) + 
  geom_point(data = ddm_n2nd %>% filter(mp_h2 > 0, !is.na(mp_h2),cluster == "normal"), aes(x=mp_t2, y = mp_h2))
ggsave("plots/inoc_5_normal_with2nd_peak.pdf")

# which not double strains have second peaks? 
n_2nd <- para %>% filter(!cluster == "double", mp_h2 > 0, drytime ==0, inocl == 5)
ddm_n2nd<- left_join(ts, n_2nd[,c("strain", "rep", "mp_h2","mp_t2")], by = c("strain", "rep"))
ggplot(ddm_n2nd %>% filter(mp_h2 > 0, drytime ==0, inoc == 5, !cluster == "double"), aes(x=Time, y = value_J, group = interaction(strain,rep))) + 
  geom_line(aes(col = factor(cluster))) + facet_wrap(~strain) + 
  geom_point(data = ddm_n2nd %>% filter(mp_h2 > 0, !cluster == "double"), aes(x=mp_t2, y = mp_h2))
ggsave("plots/inoc_5_not_double_with2nd_peak.pdf")

# Just inoc 5 / drytime == 0

g <- ggplot(cluster_data%>% filter(drytime ==0, inocl == 5), aes(x = cluster, y = value)) + geom_boxplot(aes(fill = name)) +  facet_wrap(~name, scales = "free") + 
  stat_compare_means(comparisons = my_comparisons)
ggsave("plots/inoc_5_clusters_boxplot.pdf")
