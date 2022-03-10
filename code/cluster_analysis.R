##### First look 

### Analysis of growth data
# Look at just the 10^5 inoculum: what are the clusters of behaviours? 


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

## where is the data? These are the outputs from 1_data_set#.R: standardized all variable names etc in here
ddm1 <- as.data.table(read.csv("data/ddm_set1.csv")[,-1])
ddm2 <- as.data.table(read.csv("data/ddm_set2.csv")[,-1])
ddm3 <- as.data.table(read.csv("data/ddm_set3.csv")[,-1])
ddm4 <- as.data.table(read.csv("data/ddm_set4.csv")[,-1])
ddm5 <- as.data.table(read.csv("data/ddm_set5.csv")[,-1])
ddm6 <- as.data.table(read.csv("data/ddm_set6.csv")[,-1])
ddm7 <- as.data.table(read.csv("data/ddm_set7.csv")[,-1])
ddm8 <- as.data.table(read.csv("data/ddm_set8.csv")[,-1])
ddm9 <- as.data.table(read.csv("data/ddm_set9.csv")[,-1])
ddm10 <- as.data.table(read.csv("data/ddm_set10.csv")[,-1])
ddm11 <- as.data.table(read.csv("data/ddm_set11.csv")[,-1])
ddm12 <- as.data.table(read.csv("data/ddm_set12.csv")[,-1])
ddm13 <- as.data.table(read.csv("data/ddm_set13.csv")[,-1])
ddm14 <- as.data.table(read.csv("data/ddm_set14.csv")[,-1])
ddm <- as.data.frame(rbind(ddm1,ddm2,ddm3,ddm4,ddm5,ddm6,ddm7,ddm8,ddm9,ddm10,ddm11,ddm12,ddm13,ddm14) )

length(unique(ddm$strain)) 

#### Label MACOTRA vs. not
ddm$source <- "Macotra"
w<-which(ddm$strain %in% c("Newman","RWW12","SA3297","SA2704","RWW146","SAC042W", "Mu50", "M116"))
ddm[w,"source"] <- "Other"

###******** UNITS / DATA *************#######################################################################################################################################
ddm$value_J = ddm$value

###******* MODEL FITTING *************#######################################################################################################################################
## Fit growth curves to each set of data to determine the underlying parameters 
# Just 10^5 inoculum and pre-drying
ddm_5 <- ddm %>% filter(inoc == 5, drytime == 0, source == "Macotra")

#ddm_5 <- ddm_5 %>% filter(strain %in% c(11274, 11271,11001,11006))

## What are the strains?
u <- as.character(unique(ddm_5$strain,stringsasFactors = FALSE)) # strains
## How many replicates? 
r <- unique(ddm_5$rep) # replicates
# What are the inoculums? 
q <- unique(ddm_5$inoc)

# Where the parameters for each strain are stored
param <- matrix(0, length(u)*length(r)*length(q)*3, 51); # number of strains x number of replicates x number of experimental conditions
index <- 1 # for counting 
max <- c() # for calibration
drying_times <- c(0,168)

## Run thru each strain/rep/drying time/ inoculum: fit a growth curve to each
for(jj in 1:length(u)){ # for each strain
  r <- unique(ddm_5 %>% filter(strain == u[jj]) %>% dplyr::select(rep))[,1] # which replicates for this strain (all have different names)
  if(length(r) < 3){print(paste0(u[jj], " has too few replicates"))}else{ # filter out if only 2 replicates
    for(ii in 1:length(r)){ # for each replicate
      for(kk in c(1,2)){ #each of the experimental conditions
        for(ll in 1:length(q)){ #each of the inocula
          
          strain <- u[jj];
          replicate <- r[ii]
          condition <- drying_times[kk]
          inocl <- q[ll]
          
          wi <- intersect(which(ddm_5$strain == strain),which(ddm_5$rep == replicate)) # if fit to each replicate
          wj <- intersect(wi, which(ddm_5$drytime == condition))
          w <- intersect(wj, which(ddm_5$inoc == as.numeric(inocl)))
          
          if(length(w) > 0){ # if this replicate exists for this strain (i.e. there is data)
            data1 <- ddm_5[w,] # Grab data
            
            print(c(jj, strain, replicate, condition, inocl)) # output so can track how it is working
            p <- cut_extract_dp(data1, "Time", "value_J", paste(strain, replicate, condition, inocl,sep="_")) ### NEW function: runs on any timeseries: gives baseline and cut parameter analysis
            
            ## Required parameters
            
            param[index,] <- c(strain, replicate, condition, inocl, p$param)
            index <- index + 1 # counting for storing matrix - next row each iteration
          }
          
        }
      }
    }
  }
}
# The original set 
param_orig <- param

## Fitted parameters
param <- as.data.frame(param)
colnames(param) <- c("strain","rep","drytime","inocl",
                     "t_m_h_flow", "v_m_h_flow", "exp_gr","lag","auc",
                     "odd_peaks","odd_width","width_peak","odd_shoulder","odd_shoulder_past","odd_double",
                     "shoulder_point_t","shoulder_point_v", "shoulder_point_past_t","shoulder_point_past_v",
                     "cut_exp", "timepeak", "valpeak",
                     "mp_t1","mp_t2","mp_t3","mp_t4","mp_t5","mp_t6","mp_t7","mp_t8","mp_t9","mp_t10",
                     "mp_h1","mp_h2","mp_h3","mp_h4","mp_h5","mp_h6","mp_h7","mp_h8","mp_h9","mp_h10",
                     "gap1","gap2","gap3","gap4","gap5","gap6","gap7","gap8","gap9")

w<-which(param$lag!=0); param <- param[w,] # remove 0 values
param$rep <- as.numeric(param$rep)
param$odd_peaks <- as.numeric(param$odd_peaks)
param$odd_width <- as.numeric(param$odd_width)
param$odd_shoulder <- as.numeric(param$odd_shoulder)
param$odd_shoulder_past <- as.numeric(param$odd_shoulder_past)
param$auc <- as.numeric(param$auc)
dim(param)
## Store so don't have to run above
write.csv(param, "output/param_multiple_peaks_5.csv")


## Look at auc
g1 <- ggplot(param, aes(x = width_peak , y = auc)) + geom_point() + ggtitle("Width peak")
g2 <- ggplot(param, aes(x = as.numeric(timepeak) , y = auc)) + geom_point() + ggtitle("Time peak")
g3 <- ggplot(param, aes(x = as.numeric(valpeak) , y = auc)) + geom_point() + ggtitle("Value peak")
g1 + g2 + g3
ggsave("plots/auc.pdf")

mean_auc = mean(param$auc)
max_auc = max(param$auc)
big_auc = param %>% filter(param$auc > max_auc - 0.25*max_auc) %>% summarise(unique(strain))

### how many had double peaks? 
param %>% filter(odd_peaks > 0) %>% summarise(unique(strain))

strains_in <- as.numeric(unlist(param %>% summarise(unique(strain))))
### Plot
g <- ggplot(ddm_5, aes(x=Time, y = value_J, group = rep)) + geom_line() + facet_wrap(~strain)
ggsave("plots/inoc_5.pdf")


#### Systematically remove groupings
ddm_5$cluster <- 0; param$cluster <- 0
ddm_5 <- ddm_5 %>% filter(strain %in% strains_in)
# (1) clear double peaks
double_peak_curves_analysis <- as.data.frame(param %>% filter(odd_peaks == 1) %>% group_by(strain) %>% mutate(n_odd_peaks = n()) )
ggplot(double_peak_curves_analysis %>% group_by(strain) %>% slice(1), aes(n_odd_peaks)) + geom_histogram(binwidth = 1)
nd <- left_join(as.data.frame(ddm_5), double_peak_curves_analysis[,c("strain","rep","odd_peaks","n_odd_peaks")], by = c("strain", "rep"))
g <- ggplot(nd, aes(x=Time, y = value_J, group = rep)) + geom_line(aes(col = factor(n_odd_peaks))) + facet_wrap(~strain)
ggsave("plots/inoc_5_double_peaks.pdf")

peaks_double <- double_peak_curves_analysis %>% filter(n_odd_peaks == 3) %>% ungroup() %>% summarise(unique(strain)) # double peaks in all replicates
length(u) - dim(peaks_double)[1] # left to assign
ddm_5[which(ddm_5$strain %in% unlist(peaks_double)), "cluster"] = "double"
param[which(param$strain %in% unlist(peaks_double)), "cluster"] = "double"


# (2) clear "normal peaks"
normal_peak_curves_analysis <- as.data.frame(param %>% mutate(odd = odd_peaks + odd_width + odd_shoulder + odd_shoulder_past) %>% filter(odd == 0) %>% 
                                               group_by(strain) %>% mutate(n_normal = n()))
ggplot(normal_peak_curves_analysis %>% group_by(strain) %>% slice(1), aes(n_normal)) + geom_histogram(binwidth = 1)
nd <- left_join(as.data.frame(ddm_5), normal_peak_curves_analysis[,c("strain","rep","n_normal")], by = c("strain", "rep"))
g <- ggplot(nd, aes(x=Time, y = value_J, group = rep)) + geom_line(aes(col = factor(n_normal))) + facet_wrap(~strain)
ggsave("plots/inoc_5_normal_peaks.pdf")

peaks_normal <- normal_peak_curves_analysis %>% filter(n_normal > 1) %>% ungroup() %>% summarise(unique(strain)) # double peaks in two or all replicates
length(u) - dim(peaks_double)[1] - dim(peaks_normal)[1] # left to assign
ddm_5[which(ddm_5$strain %in% unlist(peaks_normal)), "cluster"] = "normal"
param[which(param$strain %in% unlist(peaks_normal)), "cluster"] = "normal"

# (3) Some have a single rep that is a "peak" and then rest are "shoulders" => basically multiple peaks
spiked_analysis <- as.data.frame(param %>% mutate(odd_spike = odd_peaks + odd_shoulder) %>% filter(odd_spike > 0) %>% group_by(strain) %>% mutate(n_spike = n())) # how many have spikes in all? 
ggplot(spiked_analysis %>% group_by(strain) %>% slice(1), aes(n_spike)) + geom_histogram(binwidth = 1)
nd <- left_join(as.data.frame(ddm_5), spiked_analysis[,c("strain","rep","n_spike")], by = c("strain", "rep"))
g <- ggplot(nd, aes(x=Time, y = value_J, group = rep)) + geom_line(aes(col = factor(n_spike))) + facet_wrap(~strain)
ggsave("plots/inoc_5_spike_peaks.pdf")

spiked_plus_analysis <- as.data.frame(param %>% mutate(odd_spike = odd_peaks + odd_shoulder + odd_shoulder_past) %>% filter(odd_spike > 0) %>% group_by(strain) %>% mutate(n_spike = n())) # how many have spikes in all? 
ggplot(spiked_plus_analysis %>% group_by(strain) %>% slice(1), aes(n_spike)) + geom_histogram(binwidth = 1)
nd <- left_join(as.data.frame(ddm_5), spiked_plus_analysis[,c("strain","rep","n_spike")], by = c("strain", "rep"))
g <- ggplot(nd, aes(x=Time, y = value_J, group = rep)) + geom_line(aes(col = factor(n_spike))) + facet_wrap(~strain)
ggsave("plots/inoc_5_spike_plus_peaks.pdf")

peaks_spike <- spiked_analysis %>% filter(n_spike > 1) %>% filter(!strain %in% unlist(peaks_double)) %>% ungroup() %>% summarise(unique(strain)) # 2 or more have spike
peaks_spike_plus <- spiked_plus_analysis %>% filter(n_spike > 1) %>% filter(!strain %in% unlist(peaks_double)) %>% ungroup() %>% summarise(unique(strain)) # 2 or more have spike
length(u) - dim(peaks_double)[1] - dim(peaks_normal)[1] - dim(peaks_spike)[1] # left to assign
ddm_5[which(ddm_5$strain %in% unlist(peaks_spike)), "cluster"] = "spike"
param[which(param$strain %in% unlist(peaks_spike)), "cluster"] = "spike"

## Those with assignation
clustered = unlist(ddm_5 %>% filter(!cluster == 0) %>% summarise(unique(strain)))

# (4) Clear post shoulders
post_should_analysis <- as.data.frame(param %>% filter(odd_shoulder_past == 1) %>% group_by(strain) %>% mutate(n_spast= n())) # how many have spikes in all? 
ggplot(post_should_analysis %>% group_by(strain) %>% slice(1), aes(n_spast)) + geom_histogram(binwidth = 1)
nd <- left_join(as.data.frame(ddm_5), post_should_analysis[,c("strain","rep","n_spast")], by = c("strain", "rep"))
#nd$n_spast <- nd$n_spast %>% replace_na(0)
g <- ggplot(nd, aes(x=Time, y = value_J, group = rep)) + geom_line(aes(col = factor(n_spast))) + facet_wrap(~strain)
ggsave("plots/inoc_5_shoulder_past.pdf")

peaks_post_shoulder <- post_should_analysis %>% filter(n_spast > 1) %>% filter(!strain %in% clustered) %>% ungroup() %>% summarise(unique(strain)) # 2 or more have post shoulders
length(u) - dim(peaks_post_shoulder)[1] - length(clustered) # left to assign
ddm_5[which(ddm_5$strain %in% unlist(peaks_post_shoulder)), "cluster"] = "post_shoulder"
param[which(param$strain %in% unlist(peaks_post_shoulder)), "cluster"] = "post_shoulder"

## Those with assignation
clustered = unlist(ddm_5 %>% filter(!cluster == 0) %>% summarise(unique(strain)))

# (5) Clear pre shoulders
pre_should_analysis <- as.data.frame(param %>% filter(odd_shoulder == 1) %>% group_by(strain) %>% mutate(n_spre= n())) # how many have spikes in all? 
ggplot(pre_should_analysis %>% group_by(strain) %>% slice(1), aes(n_spre)) + geom_histogram(binwidth = 1)
nd <- left_join(as.data.frame(ddm_5), pre_should_analysis[,c("strain","rep","n_spre")], by = c("strain", "rep"))
g <- ggplot(nd, aes(x=Time, y = value_J, group = rep)) + geom_line(aes(col = factor(n_spre))) + facet_wrap(~strain)
ggsave("plots/inoc_5_shoulder_pre.pdf")

# (6) Wide
wide_analysis <- as.data.frame(param %>% filter(odd_width == 1) %>% group_by(strain) %>% mutate(n_wide = n())) # how many have spikes in all? 
ggplot(wide_analysis %>% group_by(strain) %>% slice(1), aes(n_wide)) + geom_histogram(binwidth = 1)
nd <- left_join(as.data.frame(ddm_5), wide_analysis[,c("strain","rep","n_wide")], by = c("strain", "rep"))
g <- ggplot(nd, aes(x=Time, y = value_J, group = rep)) + geom_line(aes(col = factor(n_wide))) + facet_wrap(~strain)
ggsave("plots/inoc_5_wide.pdf")

peaks_wide <- wide_analysis %>% filter(n_wide > 1) %>% filter(!strain %in% clustered) %>% ungroup() %>% summarise(unique(strain)) # 2 or more are wide
length(u) - dim(peaks_wide)[1] - length(clustered) # left to assign
ddm_5[which(ddm_5$strain %in% unlist(peaks_wide)), "cluster"] = "wide"
param[which(param$strain %in% unlist(peaks_wide)), "cluster"] = "wide"

# Remainder... 
odd_normal_peak_curves_analysis <- as.data.frame(param %>% mutate(odd = odd_peaks + odd_width + odd_shoulder + odd_shoulder_past, 
                                                                  odd_label = paste(odd_peaks,odd_width,odd_shoulder,odd_shoulder_past)))
ggplot(odd_normal_peak_curves_analysis, aes(odd)) + geom_histogram(binwidth = 1)
nd <- left_join(as.data.frame(ddm_5 %>% filter(is.na(cluster))), odd_normal_peak_curves_analysis[,c("strain","rep","odd","odd_label")], by = c("strain", "rep"))
g <- ggplot(nd, aes(x=Time, y = value_J, group = rep)) + geom_line(aes(col = factor(odd_label))) + facet_wrap(~strain) + 
  scale_color_discrete("peak / width / shoulder / shoulder past")
ggsave("plots/inoc_5_normal_peaks_many.pdf")

### Plot
g <- ggplot(ddm_5 %>% ungroup() %>% arrange(desc(cluster)), aes(x=Time, y = value_J, group = interaction(strain,rep))) + geom_line(aes(col = factor(cluster))) + facet_wrap(~strain)
ggsave("plots/inoc_5_clustered.pdf")

g <- ggplot(ddm_5 %>% ungroup() %>% arrange(desc(cluster)), aes(x=Time, y = value_J, group = interaction(strain,rep))) + geom_line(aes(col = factor(cluster))) + facet_grid(cluster~.)
ggsave("plots/inoc_5_clustered_rows.pdf")

g <- ggplot(ddm_5, aes(x=Time, y = value_J, group = rep)) + geom_line(aes(col = factor(cluster))) + facet_wrap(~strain)
ggsave("plots/inoc_5_grouped.pdf")

ggplot(ddm_5 %>% filter(cluster == 0), aes(x=Time, y = value_J, group = rep)) + geom_line(aes(col = factor(cluster))) + facet_wrap(~strain)

table(param$cluster)

## Auc by cluster
param$auc <- as.numeric(param$auc)
param$width_peak <- as.numeric(param$width_peak)
g1 <- ggplot(param %>% filter(!cluster == 0), aes(x = width_peak , y = auc)) + geom_point(aes(col = factor(cluster))) + ggtitle("Width peak")+ 
  scale_color_manual(breaks = c("double","spike","normal","post_shoulder","wide"), values = c("red","deeppink","blue","green","yellow"))
g2 <- ggplot(param %>% filter(!cluster == 0), aes(x = as.numeric(timepeak) , y = auc)) + geom_point(aes(col = factor(cluster))) + ggtitle("Time peak")+ 
  scale_color_manual(breaks = c("double","spike","normal","post_shoulder","wide"), values = c("red","deeppink","blue","green","yellow"))
g3 <- ggplot(param %>% filter(!cluster == 0), aes(x = as.numeric(valpeak) , y = auc)) + geom_point(aes(col = factor(cluster))) + ggtitle("Value peak") + 
  scale_color_manual(breaks = c("double","spike","normal","post_shoulder","wide"), values = c("red","deeppink","blue","green","yellow"))
g1 + g2 + g3 + plot_layout(guides = "collect") 
ggsave("plots/auc_by_cluster.pdf", width = 20, height = 8)



## Correct for latent period 
ddm_5_latent <- left_join(ddm_5, param[,c("strain", "rep", "lag")], by = c("strain", "rep"))
ddm_5_latent$value_latent_correct <- ddm_5_latent$Time - as.numeric(ddm_5_latent$lag)
ddm_5_latent <- ddm_5_latent %>% filter(value_latent_correct > 0)
ggplot(ddm_5_latent %>% ungroup() %>% arrange(desc(cluster)), aes(x=value_latent_correct, y = value_J, group = interaction(strain,rep))) + geom_line(aes(col = factor(cluster))) + facet_grid(cluster~.)
ggsave("plots/inoc_5_clustered_rows_correct_latent.pdf")



###### Compare properties across clusters
param$valpeak <- as.numeric(param$valpeak)
param$timepeak <- as.numeric(param$timepeak)
param$exp_gr <- as.numeric(param$exp_gr)
param$shoulder_point_v <- as.numeric(param$shoulder_point_v)
param$shoulder_point_past_v <- as.numeric(param$shoulder_point_past_v)
param$shoulder_point_t<- as.numeric(param$shoulder_point_t)
param$shoulder_point_past_t <- as.numeric(param$shoulder_point_past_t)
param$mp_h2 <- as.numeric(param$mp_h2)

cluster_data <- param %>% dplyr::select(c("cluster","valpeak", "timepeak", "auc", "exp_gr", "shoulder_point_v","shoulder_point_t",
                                          "shoulder_point_past_v","shoulder_point_past_t", "mp_h2")) %>% 
  pivot_longer(cols = c(valpeak:mp_h2)) 

cluster_data_summ <- cluster_data %>% group_by(cluster, name) %>% summarise(mean = mean(value), sd = sd(value))

ggplot(cluster_data_summ, aes(x = cluster, y = mean, fill = name)) + geom_bar(stat="identity", color="black", 
                                                                              position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) + 
  facet_wrap(~name, scales = "free")
ggsave("plots/clusters_mean_sd.pdf")

ggplot(cluster_data, aes(x = cluster, y = value, fill = name)) + geom_boxplot() + facet_wrap(~name, scales = "free")
ggsave("plots/clusters_boxplot.pdf")


my_comparisons <- list( c("post_shoulder", "double"), c("double", "normal"), c("normal", "post_shoulder"),c("normal","spike"),c("spike","wide"),
                        c("double","spike"),c("double","wide"))

ggplot(cluster_data %>% filter(name == "valpeak"), aes(x = cluster, y = value, fill = name)) + geom_boxplot() + 
  facet_wrap(~name, scales = "free") + 
  stat_compare_means(comparisons = my_comparisons)
ggsave("plots/clusters_boxplot_comps_valpeak.pdf")

ggplot(cluster_data, aes(x = cluster, y = value)) + geom_boxplot(aes(fill = name)) + 
  facet_wrap(~name, scales = "free") + 
  stat_compare_means(comparisons = my_comparisons)
ggsave("plots/clusters_boxplot_comps_all.pdf")


#### Oddity in some normal have later peak 
ggplot(param %>% filter(cluster == "normal"), aes(x=timepeak,y=valpeak )) + geom_point(aes(col = factor(strain))) + theme(legend.position = "none")
param %>% filter(cluster == "normal",timepeak > 10)%>% summarise(unique(strain))

## Look at width
g1 <- ggplot(param, aes(x = auc , y = width_peak)) + geom_point(aes(col = factor(cluster)))  + ggtitle("AUC")+ 
  scale_color_manual(breaks = c("double","spike","normal","post_shoulder","wide","0"), values = c("red","deeppink","blue","green","yellow","black")) + geom_hline(yintercept = 88)
g2 <- ggplot(param, aes(x = as.numeric(timepeak) , y = width_peak)) + geom_point(aes(col = factor(cluster))) + ggtitle("Time peak")+ 
  scale_color_manual(breaks = c("double","spike","normal","post_shoulder","wide","0"), values = c("red","deeppink","blue","green","yellow","black"))+ geom_hline(yintercept = 88)
g3 <- ggplot(param, aes(x = as.numeric(valpeak) , y = width_peak)) + geom_point(aes(col = factor(cluster))) + ggtitle("Value peak")+ 
  scale_color_manual(breaks = c("double","spike","normal","post_shoulder","wide","0"), values = c("red","deeppink","blue","green","yellow","black"))+ geom_hline(yintercept = 88)
g1 + g2 + g3
ggsave("plots/width.pdf", width = 20, height = 8)

g1 <- ggplot(param, aes(x = auc , y = width_peak)) + geom_point(aes(col = factor(strain)))  + ggtitle("AUC")+ theme(legend.position = "none") + geom_hline(yintercept = 88)
g2 <- ggplot(param, aes(x = as.numeric(timepeak) , y = width_peak)) + geom_point(aes(col = factor(strain))) + ggtitle("Time peak")+ theme(legend.position = "none")+ geom_hline(yintercept = 88)
g3 <- ggplot(param, aes(x = as.numeric(valpeak) , y = width_peak)) + geom_point(aes(col = factor(strain))) + ggtitle("Value peak")+ theme(legend.position = "none")+ geom_hline(yintercept = 88)
g1 + g2 + g3
ggsave("plots/width_strains.pdf", width = 20, height = 8)
