##### Analysis of non-growth related heat flow
## For this using a comparison of OD and CS data

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
library(zoo)
library(imputeTS) # to give weighter moving average - exponential weighting of those further away

theme_set(theme_bw(base_size=12)) # theme setting for plots: black and white (bw) and font size (24)
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(14)

setwd(here::here())

#### Load functions
source("code/functions_for_heat_curves_additional_double_peak.R")


#######################********** Input data****************##############################################################################################################################################################
#### INPUT data read in and explore
############################################################################################################################################################################################################################
data_od_orig <- read_csv("data/growth_ODvsCS_20220224.csv")[,-1]
# remove all data before the time point that they all have which is the max of the minimum times to avoid odd completion curves prior to start and after
cutoff_time_dn = max(data_od_orig %>% dplyr::group_by(rep, exp, inoc) %>% dplyr::summarise(minn = min(Time)) %>% ungroup() %>% dplyr::select(minn))
cutoff_time_up = min(data_od_orig %>% dplyr::group_by(rep, exp, inoc) %>% dplyr::summarise(maxx = max(Time)) %>% ungroup() %>% dplyr::select(maxx))

## Look at OD & CS data
ggplot(data_od_orig, aes(x=Time, y = value, group = interaction(inoc, exp, strain, rep))) + geom_line(aes(col = factor(inoc))) + 
  facet_grid(exp~strain, scales = "free")
ggsave("plots/ODvsCS_rawdata.pdf")



#######################********** ALL ****************##############################################################################################################################################################
#### All inocula
############################################################################################################################################################################################################################

ggplot(data_od_orig, aes(x=Time, y = value, group = interaction(inoc, exp, strain, rep))) + geom_line(aes(col = factor(rep))) + 
  facet_grid(exp + inoc~strain, scales = "free")
ggsave("plots/ODvsC2_rawdata.pdf")

data_od <- data_od_orig %>% filter(Time > cutoff_time_dn, Time < cutoff_time_up) %>% 
  dplyr::select(Time, rep, exp, value, strain, inoc) %>% 
  dplyr::group_by(rep, exp, strain, inoc) %>%
  dplyr::mutate(group = interaction(rep, exp, strain,inoc))%>%
  dplyr::group_by(group) %>%
  dplyr::mutate(ma_value = rollapply(value, 10, mean,fill = NA,align = "right", partial = TRUE), #### NOT RESPECTING GROUPS??
                differ = c(0,diff(ma_value)),
                compara = ifelse(exp == "CS", value, differ)) %>%
  dplyr::ungroup() %>% filter(Time > cutoff_time_dn, Time < cutoff_time_up)

g2 <- ggplot(data_od, aes(x=Time, y = ma_value, group = interaction(rep, exp, strain))) + 
  geom_line(aes(col = factor(rep))) + 
  facet_grid(exp + inoc~strain,scales = "free") + 
  scale_x_continuous("Time (h)") 
ggsave("plots/ODvsCS_data_smoothed.pdf")

g3 <- ggplot(data_od, aes(x=Time, y = compara, group = interaction(rep, exp, strain))) + 
  geom_line(aes(col = factor(rep),lty = factor(exp))) + 
  facet_grid(inoc + exp~strain, scales = "free") + 
  scale_x_continuous("Time (h)") + scale_color_discrete("Replicate") + scale_linetype_discrete("Data source")
ggsave("plots/ODvsCS_data_compare.pdf",width = 14)

# Try to plot together? scale by max. Normalise
max_vals_norm <- data_od %>% filter(strain == "11257") %>% group_by(inoc, rep, exp )%>% summarise(max_norms = max(compara, na.rm = TRUE))

data_od <- left_join(data_od, max_vals_norm) %>% mutate(compara_norm = compara / max_norms)

ggplot(data_od, aes(x=Time, y = compara_norm, group = interaction(rep, exp, strain))) + 
  geom_line(aes(col = interaction(rep), lty = exp), lwd = 1) + 
  facet_grid(inoc ~strain, scales = "free") + 
  scale_x_continuous("Time (h)") 
ggsave("plots/ODvsCS_data_compare_norm.pdf")

# Subtract normalised data? Need to complete: measured at different time points

data_od_normd <- data_od %>% ungroup() %>% filter(Time > cutoff_time_dn, Time < cutoff_time_up) %>% # make sure cover same time for all 
  complete(rep, strain, exp, inoc, Time) %>% dplyr::mutate(compara_norm_inp = na_ma(compara_norm, k = 4, weighting = "linear", maxgap = 10)) %>% # fill in all time points and then linear imputation between (tried exponential and simple but get more odd bumps)
  dplyr::group_by(inoc, strain, rep, exp) %>% dplyr::select(Time, strain,rep, exp, inoc, compara_norm_inp) %>% # Take imputed values
  pivot_wider(id_cols = c(strain, Time, rep, inoc), names_from = exp, values_from = compara_norm_inp) %>%dplyr::mutate(nongrowth_only = CS - OD) %>% # look for difference between OD and heat output
  pivot_longer(cols = c("CS","OD"), names_to = "exp", values_to ="imput_val")

# Add in time to peak value
data_od <- data_od %>% group_by(strain, rep, exp, inoc) %>% mutate(peak = max(compara), equals_peak = ifelse(compara == peak, 1, 0), times_peak = ifelse(equals_peak == 1, Time, 0), time_peak = max(times_peak))
data_od_normd_ana <- left_join(data_od, data_od_normd)

ggplot(data_od_normd_ana, aes(x=Time, y = nongrowth_only, group = interaction(exp,rep, strain))) + 
  geom_line(aes(col = interaction(exp,rep)), lwd = 1) + 
  geom_line(aes(y = imput_val,col = interaction(exp,rep))) + 
  facet_grid(inoc~strain, scales = "free") + 
  scale_x_continuous("Time (h)") 
ggsave("plots/ODvsCS_data_nongrowth_togplot.pdf")


g1 <- ggplot(data_od_normd_ana, aes(x=Time, y = nongrowth_only, group = interaction(rep, strain))) + 
  geom_line(aes(col = interaction(rep))) + 
  facet_grid(inoc~strain) + 
  scale_x_continuous("Time (h)") + 
  scale_y_continuous("Non growth only") + 
  scale_color_discrete("Replicate") + 
  geom_hline(yintercept = 0)

g2 <- ggplot(data_od_normd_ana, aes(x=Time, group = interaction(exp,rep, strain))) + 
  geom_line(aes(y = imput_val,col = interaction(rep), linetype = exp)) + 
  facet_grid(inoc~strain, scales = "free") + 
  scale_x_continuous("Time (h)") + 
  scale_y_continuous("Normalised measure") + 
  scale_color_discrete("Experiment and\nreplicate")

g2 / g1 
ggsave("plots/ODvsCS_data_nongrowth_tog_grid.pdf")

#######################********** 10^5 ****************##############################################################################################################################################################
#### JUST DO FOR 10^5 to start! 
############################################################################################################################################################################################################################

#### Explore data and extract 
data_exploration_od_cs(data_od_orig %>% filter(inoc == 5), "inoc5")


g1 <- ggplot(data_od %>% filter(inoc == 5), aes(x=Time, y = compara, group = interaction(rep, exp, strain))) + 
  geom_line(aes(col = factor(rep),lty = factor(exp))) + 
  geom_vline(aes(xintercept = time_peak, col = factor(rep)), alpha = 0.3) + 
  facet_grid(exp~strain, scales = "free") + 
  scale_x_continuous("Time (h)") + scale_color_discrete("Replicate") + scale_linetype_discrete("Data source") + 
  scale_y_continuous("Raw data") 

g2 <- ggplot(data_od %>% filter(inoc == 5), aes(x=Time, y = compara_norm, group = interaction(rep, exp, strain))) + 
  geom_line(aes(col = interaction(rep), lty = exp)) + 
  geom_vline(aes(xintercept = time_peak, col = factor(rep)), alpha = 0.3) + 
  facet_grid(~strain, scales = "free") + 
  scale_x_continuous("Time (h)") + scale_color_discrete("Replicate") + scale_linetype_discrete("Data source") + scale_y_continuous("Normalised measure")

g3 <- ggplot(data_od_normd_ana %>% filter(inoc == 5), aes(x=Time, y = nongrowth_only, group = interaction(rep, strain))) + 
  geom_line(aes(col = interaction(rep))) + 
  facet_grid(~strain) + 
  geom_vline(aes(xintercept = time_peak, col = factor(rep)), alpha = 0.3) + 
  scale_x_continuous("Time (h)") + 
  scale_y_continuous("Non growth only") + 
  scale_color_discrete("Replicate") + 
  geom_hline(yintercept = 0)

(g1/g2 + plot_layout(guides = "collect"))/g3 + plot_layout(guides = "collect")
ggsave("plots/final/figure4.png",width = 8, height = 6, dpi = 600)

#######################********** Inoculum effect ****************##############################################################################################################################################################
#### Look at inoculum effect
############################################################################################################################################################################################################################

ggplot(data_od_normd_ana, aes(x=Time, y = nongrowth_only, group = interaction(inoc, rep, strain))) + 
  geom_line(aes(col = factor(inoc))) + 
  facet_grid(rep~strain) + 
  scale_x_continuous("Time (h)") + 
  scale_y_continuous("Non growth only") + 
  scale_color_discrete("Replicate") + 
  geom_hline(yintercept = 0)

max_nongrowth_tab <- data_od_normd_ana %>% group_by(rep, inoc, strain) %>% filter(nongrowth_only == max(nongrowth_only))

ggplot(max_nongrowth_tab, aes(x=inoc, y = nongrowth_only)) + geom_point(aes(col = factor(strain))) + 
  geom_line(aes(group = interaction(strain, rep), col = factor(strain))) + 
  scale_y_continuous("Maximum non-growth") + scale_x_continuous("Inoculum")
ggsave("plots/ODvsCS_inoc_vs_max_nongrowth.pdf")

ggplot(max_nongrowth_tab, aes(x=inoc, y = Time)) + geom_point(aes(col = factor(strain))) + 
  geom_line(aes(group = interaction(strain, rep), col = factor(strain))) + 
  scale_y_continuous("Maximum non-growth") + scale_x_continuous("Inoculum")
ggsave("plots/ODvsCS_inoc_vs_time_max_nongrowth.pdf")

ggplot(max_nongrowth_tab, aes(x=inoc, y = nongrowth_only)) + geom_point(aes(col = factor(strain))) + 
  geom_smooth(aes(group = interaction(strain), col = factor(strain), fill = factor(strain))) + 
  scale_y_continuous("Maximum non-growth") + scale_x_continuous("Inoculum")
ggsave("plots/ODvsCS_inoc_vs_max_nongrowth_smooth.pdf")

ggplot(max_nongrowth_tab, aes(x=inoc, y = Time)) + geom_point(aes(col = factor(strain))) + 
  geom_smooth(aes(group = interaction(strain), col = factor(strain), fill = factor(strain))) + 
  scale_y_continuous("Maximum non-growth") + scale_x_continuous("Inoculum")
ggsave("plots/ODvsCS_inoc_vs_time_max_nongrowth_smooth.pdf")


#######################********** Extract characteristics ****************##############################################################################################################################################################
#### Extract key parameters from the non-growth data
############################################################################################################################################################################################################################


## What are the strains?
u <- as.character(unique(data_od$strain,stringsasFactors = FALSE)) # strains
## How many replicates? 
r <- unique(data_od$rep) # replicates
## How many inoc?
inn <- unique(data_od$inoc)

# Where the parameters for each strain are stored
param <- matrix(0, length(u)*length(r)*length(inn), 9); # number of strains x number of replicates x number of experimental conditions
index <- 1 # for counting 
max <- c() # for calibration

data_od <- data_od %>% dplyr::ungroup()

# ## Remove odd one that doesn't peak.. 
# w <- intersect(which(data_od$strain == 11051), which(data_od$inoc == 1))
# data_od <- data_od[-w,]

## Run thru each strain/rep/drying time/ inoculum: fit a growth curve to each
for(jj in 1:length(u)){ # for each strain
  r <- unlist(unique(data_od %>% filter(strain == u[jj]) %>% dplyr::select(rep))[,1]) # which replicates for this strain (all have different names)
  if(length(r) < 3){print(paste0(u[jj], " has too few replicates"))}else{ # filter out if only 2 replicates
    for(ii in 1:length(r)){ # for each replicate
      for(pp in 1:length(inn)){ # for each inoc
        
        data1 <- data_od %>% filter(strain == u[jj], rep == r[ii], inoc == inn[pp], exp == "OD")
        data2 <- data_od %>% filter(strain == u[jj], rep == r[ii], inoc == inn[pp], exp == "CS")
        data_nong <- data_od_normd_ana %>% filter(strain == u[jj], rep == r[ii], inoc == inn[pp]) 
        
        if(dim(data1)[1] > 0){ # if this replicate exists for this strain (i.e. there is data)
          
          print(c(jj, u[jj], r[ii], inn[pp],ex[kk])) # output so can track how it is working
          p <- charac_extract(data_nong, "Time", "nongrowth_only", data1, data2, "compara_norm", paste(strain, replicate, condition, inocl,sep="_")) ### NEW function: simplification of more complex data recognition for clustering
          
          ## Required parameters
          param[index,] <- c(u[jj], r[ii], inn[pp], p)
          index <- index + 1 # counting for storing matrix - next row each iteration
        }
        
      }
    }
  }
}

# The original set 
param_orig <- param

## Fitted parameters
param <- param[-which(param[,1]==0),] # remove any empty rows (when replicates < 3)
param <- as.data.frame(param)
colnames(param) <- c("strain","rep","inoc",
                     "t_max_h_flow", "v_max_h_flow", 
                     "t_min_h_flow", "v_min_h_flow", 
                     "auc","lagtime")
param$rep <- as.numeric(param$rep)
param$auc <- as.numeric(param$auc)
dim(param)

## Store so don't have to run above
write.csv(param, "output/simple_extract_odcs_allinoc.csv")
param <- read_csv("output/simple_extract_odcs_allinoc.csv")[,-1]

##########################################################################################
#### *********** Analyse above characteristics*************** ####################
##########################################################################################
## Take mean over replicates
param_mean <- param %>% group_by(strain, inoc) %>% summarise(mean_t_max = mean(t_max_h_flow),mean_h_max = mean(v_max_h_flow),
                                                             mean_t_min = mean(t_min_h_flow),mean_h_min = mean(v_min_h_flow),
                                                             mean_auc = mean(auc))

# Heatmaps over inocula
g1 <- ggplot(param_mean, aes(inoc, factor(strain), fill= mean_t_max)) + geom_tile() +
  scale_fill_distiller(palette = "RdBu","Time to\nmax.\nnon-growth", ) + 
  scale_x_continuous("Inoculum") + 
  scale_y_discrete("Strain") #+ 
  #ggtitle("Time to max non-growth output")

g2 <- ggplot(param_mean, aes(inoc, factor(strain), fill= mean_h_max)) + geom_tile() +
  scale_fill_distiller(palette = "RdBu","Height of\nmax.\nnon-growth", ) + 
  scale_x_continuous("Inoculum") + 
  scale_y_discrete("Strain") #+ 
  #ggtitle("Height of max non-growth output")

g3 <- ggplot(param_mean, aes(inoc, factor(strain), fill= mean_t_min)) + geom_tile() +
  scale_fill_distiller(palette = "RdBu","Time to\nmin.\nnon-growth", ) + 
  scale_x_continuous("Inoculum") + 
  scale_y_discrete("Strain") #+ 
  #ggtitle("Time to min non-growth output")

g4 <- ggplot(param_mean, aes(inoc, factor(strain), fill= mean_h_min)) + geom_tile() +
  scale_fill_distiller(palette = "RdBu","Height of \nmin.\nnon-growth", ) + 
  scale_x_continuous("Inoculum") + 
  scale_y_discrete("Strain") #+ 
  #ggtitle("Height of min non-growth output")

g5 <- ggplot(param_mean, aes(inoc, factor(strain), fill= mean_auc)) + geom_tile() +
  scale_fill_distiller(palette = "RdBu","AUC \nnon-growth", ) + 
  scale_x_continuous("Inoculum") + 
  scale_y_discrete("Strain") #+ 
  #ggtitle("AUC non-growth output")

g1 + g3 + g5 + g2 + g4 + plot_layout(ncol = 3)
ggsave("plots/ODvsCS_summary_heatmap.pdf", width = 12, height = 8)

#### Plot over strains
param_long <- param %>% pivot_longer(cols = t_max_h_flow:auc)
param_long$name <- factor(param_long$name, levels = c("t_max_h_flow","t_min_h_flow" ,"auc", "v_max_h_flow","v_min_h_flow"))

ggplot(param_long, aes(x=inoc, y = value, group = interaction(strain,rep))) + geom_line(aes(col = factor(strain))) + 
  facet_wrap(~name, scales = "free") + 
  scale_color_discrete("Strain") + 
  scale_x_continuous("Inoculum")
ggsave("plots/ODvsCS_summary_lines.pdf")

ggplot(param_long, aes(x=inoc, y = value, group = interaction(strain))) + geom_point(aes(col = factor(strain))) + 
  geom_smooth(aes(col = factor(strain), fill = factor(strain)),method='lm', formula= y~x) + 
  facet_wrap(~name, scales = "free") + 
  scale_color_discrete("Strain") + scale_fill_discrete("Strain") + 
  scale_x_continuous("Inoculum")
ggsave("plots/ODvsCS_summary_smoothed.pdf")


#### Timing / bump analysis over inocula on one graph
data_bump <- left_join(data_od_normd_ana %>% dplyr::select(Time, nongrowth_only,inoc, rep, strain), param) %>% mutate(nolag = Time - lagtime)

ggplot(data_bump, aes(x=nolag, y = nongrowth_only, group = interaction(rep, strain))) + 
  geom_line(aes(col = interaction(rep))) + 
  facet_grid(inoc~strain) + 
  scale_x_continuous("Time (h)") + 
  scale_y_continuous("Non growth only") + 
  scale_color_discrete("Replicate") + 
  geom_hline(yintercept = 0) 
# taking off lag doesn't work... 

# take mean over replicates
ggplot(data_bump %>% filter(inoc %in% c(2,3,4,5,6)), aes(x=Time, y = nongrowth_only, group = interaction(inoc, strain))) + 
  geom_smooth(aes(fill = factor(inoc), col = factor(inoc))) + 
  facet_grid(~strain) + 
  scale_x_continuous("Time (h)") + 
  scale_y_continuous("Non growth only") + 
  scale_color_discrete("Inoculum") + scale_fill_discrete("Inoculum")
ggsave("plots/ODvsCS_non_growth_average.pdf", width = 12, height = 8)

##

# 
# NEED TO DO THE BELOW? Treat CS and OD data as time series and explore with clustering etc? 
# 
# 
# 
# #######################********** Extract characteristics ****************##############################################################################################################################################################
# #### Extract key parameters from the CS and OD data
# ############################################################################################################################################################################################################################
# 
# 
# ## What are the strains?
# u <- as.character(unique(data_od$strain,stringsasFactors = FALSE)) # strains
# ## How many replicates? 
# r <- unique(data_od$rep) # replicates
# ## How many experimental conditions? 
# ex <- unique(data_od$exp)
# ## How many inoc?
# inn <- unique(data_od$inoc)
# 
# # Where the parameters for each strain are stored
# param <- matrix(0, length(u)*length(r)*length(ex)*length(inn), 51); # number of strains x number of replicates x number of experimental conditions
# index <- 1 # for counting 
# max <- c() # for calibration
# 
# data_od <- data_od %>% dplyr::ungroup()
# 
# ## Remove odd one that doesn't peak.. 
# w <- intersect(which(data_od$strain == 11051), which(data_od$inoc == 1))
# data_od <- data_od[-w,]
# 
# ## Run thru each strain/rep/drying time/ inoculum: fit a growth curve to each
# for(jj in 1:length(u)){ # for each strain
#   r <- unlist(unique(data_od %>% filter(strain == u[jj]) %>% dplyr::select(rep))[,1]) # which replicates for this strain (all have different names)
#   if(length(r) < 3){print(paste0(u[jj], " has too few replicates"))}else{ # filter out if only 2 replicates
#     for(ii in 1:length(r)){ # for each replicate
#       for(pp in 1:length(inn)){ # for each inoc
#         for(kk in 1:length(ex)){ #each of the experimental conditions
#           
#           data1 <- data_od %>% filter(strain == u[jj], rep == r[ii], inoc == inn[pp], exp == ex[kk])
#           
#           if(dim(data1)[1] > 0){ # if this replicate exists for this strain (i.e. there is data)
#             
#             if(length(which(is.na(data1$compara)))>0){data1 <- data1[-which(is.na(data1$compara)),]} # remove any NA
#             
#             print(c(jj, u[jj], r[ii], inn[pp],ex[kk])) # output so can track how it is working
#             p <- cut_extract_dp(data1, "Time", "compara", paste(strain, replicate, condition, inocl,sep="_")) ### NEW function: runs on any timeseries: gives baseline and cut parameter analysis
#             
#             ## Required parameters
#             
#             param[index,] <- c(u[jj], r[ii], inn[pp],ex[kk], unlist(p$param))
#             index <- index + 1 # counting for storing matrix - next row each iteration
#           }
#           
#         }
#       }
#     }
#   }
# }
# 
# # The original set 
# param_orig <- param
# 
# ## Fitted parameters
# param <- as.data.frame(param)
# colnames(param) <- c("strain","rep","inoc","experiment",
#                      "t_m_h_flow", "v_m_h_flow", "exp_gr","lag","auc",
#                      "odd_peaks","odd_width","width_peak","odd_shoulder","odd_shoulder_past","odd_double",
#                      "shoulder_point_t","shoulder_point_v", "shoulder_point_past_t","shoulder_point_past_v",
#                      "cut_exp", "timepeak", "valpeak",
#                      "mp_t1","mp_t2","mp_t3","mp_t4","mp_t5","mp_t6","mp_t7","mp_t8","mp_t9","mp_t10",
#                      "mp_h1","mp_h2","mp_h3","mp_h4","mp_h5","mp_h6","mp_h7","mp_h8","mp_h9","mp_h10",
#                      "gap1","gap2","gap3","gap4","gap5","gap6","gap7","gap8","gap9")
# 
# w<-which(param$lag!=0); param <- param[w,] # remove 0 values
# param$rep <- as.numeric(param$rep)
# param$odd_peaks <- as.numeric(param$odd_peaks)
# param$odd_width <- as.numeric(param$odd_width)
# param$odd_shoulder <- as.numeric(param$odd_shoulder)
# param$odd_shoulder_past <- as.numeric(param$odd_shoulder_past)
# param$auc <- as.numeric(param$auc)
# param$valpeak <- as.numeric(param$valpeak)
# param$timepeak <- as.numeric(param$timepeak)
# dim(param)
# 
# ## Store so don't have to run above
# write.csv(param, "output/param_od_allinoc.csv")
# param <- read_csv("output/param_od_allinoc.csv")[,-1]
# 
# #######################********** Clustering ****************##############################################################################################################################################################
# #### Cluster the data
# ############################################################################################################################################################################################################################
# 
# ### CHECK PUTTING IN RIGHT DATA AND CAN the cluster function deal with all the inocula?
# 
# 
# data_od$drytime <- 0
# param$drytime <- 0
# data_od$inoc <- 5
# param$inoc <- 5
# data_od$strain <- as.character(data_od$strain)
# param$strain <- as.character(param$strain)
# data_od$value_J <- data_od$compara
# 
# 
# c_cs <- cluster(data_od %>% filter(exp == "CS"), param %>% filter(experiment == "CS"),plot_where = "plots/CS_")
# c_od <- cluster(data_od %>% filter(exp == "OD"), param %>% filter(experiment == "OD"),plot_where = "plots/OD_")
# 
# 
# 
# 
# #### Questions
# ### Second peaks? 
# #OLD: if more than 5 from end of exp then ok - check up to third minor. mp_t1 = first peak ##ifelse(mp_h2 > 0, ifelse(mp_t2 > (timepeak + 5), mp_h2, ifelse(mp_h3>0,ifelse(mp_h3 > (timepeak + 5), mp_h3,0),0)),0)) 
# para <- para %>% ungroup() %>% mutate(second_peak_h = ifelse((t_m_h_flow > (timepeak + 2)), v_m_h_flow, ifelse(mp_t2 > (timepeak + 4), mp_h2, 0))) # if first is no max then second = max, otherwise 1st is max so take next 
# 
# 
# cluster_data <- para %>% dplyr::select(c("strain","rep","drytime","inocl","cluster","valpeak", "timepeak", "auc", "exp_gr", "shoulder_point_v","shoulder_point_t",
#                                          "shoulder_point_past_v","shoulder_point_past_t", "second_peak_h")) %>% 
#   pivot_longer(cols = c(valpeak:second_peak_h)) 
# 
# cluster_data_summ_all <- cluster_data %>% group_by(cluster, name) %>% summarise(mean = mean(value), sd = sd(value))
# 
# g <- ggplot(cluster_data_summ_all, aes(x = cluster, y = mean, fill = name)) + geom_bar(stat="identity", color="black", 
#                                                                               position=position_dodge()) +
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
#                 position=position_dodge(.9)) + 
#   facet_wrap(~name, scales = "free")
# ggsave("plots/ALL_clusters_mean_sd.pdf")
# 
# g <- ggplot(cluster_data, aes(x = cluster, y = value, fill = name)) + geom_boxplot() + facet_wrap(~name, scales = "free")
# ggsave("plots/ALL_clusters_boxplot.pdf")
# 
# #### How many of each cluster type by inoculum and drytime? 
# para2 <- para
# para2$cluster <- factor(para2$cluster, levels = c("normal","double","spike","post_shoulder","wide","no cluster"))
# howmany <- para2 %>% group_by(inocl, drytime, cluster,  .drop=FALSE) %>% summarise( n = n()) %>% complete(cluster, fill = list(n = 0))
# #ggplot(howmany, aes(x=cluster,y = n, group = interaction(inocl,drytime))) + geom_point(aes(col = interaction(inocl,drytime))) + geom_line(aes(col = interaction(inocl,drytime))) + 
# #  facet_wrap(~drytime)
# 
# ggplot(howmany, aes(x=cluster,y = n, group = interaction(inocl,drytime))) + geom_bar(stat="identity",position = position_dodge(),aes(fill = factor(inocl)))+ 
#   facet_wrap(~drytime, ncol = 1) + scale_fill_discrete("Inoculum")
# ggsave("plots/ana_clusters_by_inoc_drytime.pdf") # no double / spike at 3 / 5
# 
# ### Does normal or double have higher total energy? (AUC)
# cluster_data$cluster <- factor(cluster_data$cluster, levels = c("normal","double","spike","post_shoulder","wide","no cluster"))
# cluster_data$inocl <- factor(cluster_data$inocl)
# cluster_data$drytime <- factor(cluster_data$drytime)
# cluster_data$name <- factor(cluster_data$name)
# cluster_data_summ <- cluster_data %>% group_by(drytime, inocl, cluster, name, .drop=FALSE) %>% summarise(mean = mean(value), sd = sd(value))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ###### OLD ANALYSIS
# 
# ggplot(cluster_data_summ %>% filter(name == "auc"), aes(x = cluster, y = mean, fill = factor(inocl))) + 
#   geom_bar(stat="identity", color="black", position=position_dodge())+
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9)) + 
#   facet_wrap(~drytime , scales = "free", ncol = 1) + scale_fill_discrete("Inoculum") + 
#   ggtitle("AUC")
# ggsave("plots/ana_auc_summary.pdf")
# 
# ggplot(cluster_data_summ , aes(x = cluster, y = mean, fill = factor(inocl))) + 
#   geom_bar(stat="identity", color="black", position=position_dodge())+
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9)) + 
#   facet_grid(name~drytime , scales = "free") + scale_fill_discrete("Inoculum") 
# ggsave("plots/ana_all_summary.pdf")
# 
# 
# ggplot(cluster_data_summ %>% filter(name %in% c("valpeak","auc","timepeak")), aes(x = cluster, y = mean, fill = factor(inocl))) + 
#   geom_bar(stat="identity", color="black", position=position_dodge())+
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9)) + 
#   facet_grid(name ~ drytime, scales = "free") + scale_fill_discrete("Inoculum") 
# ggsave("plots/ana_auc&val&timepeak_summary.pdf")
# 
# ggplot(cluster_data_summ %>% filter(name %in% c("valpeak","auc","timepeak"), cluster %in% c("normal","spike","double")), aes(x = cluster, y = mean, fill = factor(inocl))) + 
#   geom_bar(stat="identity", color="black", position=position_dodge())+
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9)) + 
#   facet_grid(name ~ drytime, scales = "free") + scale_fill_discrete("Inoculum") 
# ggsave("plots/ana_auc&val&timepeak_norm_double.pdf")
# 
# my_comparisons <- list( c("post_shoulder", "double"), c("double", "normal"), c("normal", "post_shoulder"),c("normal","spike"),c("spike","wide"),
#                         c("double","spike"),c("double","wide"))
# 
# ggplot(cluster_data_summ , aes(x = cluster, y = mean, fill = factor(inocl))) + 
#   geom_bar(stat="identity", color="black", position=position_dodge())+
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9)) + 
#   facet_grid(name~drytime , scales = "free") + scale_fill_discrete("Inoculum") + 
#   stat_compare_means(comparisons = my_comparisons)
# 
# #### How does AUC look vs heat of first / second peak? 
# g1 <- ggplot(para %>% filter(drytime ==0, inocl == 5), aes(x=auc, y = valpeak, group = cluster)) + geom_point(aes(col = cluster)) + geom_smooth(method='lm',formula= y~x,se=FALSE, size=2, aes(col = cluster)) + 
#   ggtitle("AUC vs height of first peak")
# g2 <- ggplot(para %>% filter(drytime ==0, inocl == 5), aes(x=auc, y = mp_h2, group = cluster)) + geom_point(aes(col = cluster)) + geom_smooth(method='lm',formula= y~x,se=FALSE, size=2, aes(col = cluster)) + 
#   ggtitle("AUC vs height of second peak")
# g1 + g2
# ggsave("plots/inoc_5_AUC_vs_heightfp.pdf")
# 
# ggplot(para %>% filter(drytime ==0, inocl == 5), aes(x=auc, y = mp_h2, group = cluster)) + geom_point(aes(col = cluster)) + geom_smooth(method='lm',formula= y~x,se=FALSE, size=2, aes(col = cluster)) + 
#   ggtitle("AUC vs height of first peak")
# 
# ggplot(cluster_data %>% filter(drytime ==0, inocl == 5, name %in% c("mp_h2","auc","valpeak","exp_gr")) , aes(x = cluster, y = value, fill = factor(cluster))) + 
#   geom_boxplot() + 
#   facet_wrap(~name, scales = "free") + 
#   stat_compare_means(comparisons = my_comparisons)
# ggsave("plots/inoc_5_auc_vs_height_stats.pdf", width = 15, height = 6)
# 
# # which normal strains have second peaks? 
# n_2nd <- para %>% filter(mp_h2 > 0, cluster == "normal", drytime ==0, inocl == 5)
# ts$strain <- as.numeric(ts$strain)
# ddm_n2nd<- left_join(ts, n_2nd[,c("strain", "rep", "mp_h2","mp_t2")], by = c("strain", "rep"))
# ggplot(ddm_n2nd %>% filter(drytime ==0, !is.na(mp_h2), inoc == 5, mp_h2 > 0, cluster == "normal"), aes(x=Time, y = value_J, group = interaction(strain,rep))) + 
#   geom_line(aes(col = factor(cluster))) + facet_wrap(~strain) + 
#   geom_point(data = ddm_n2nd %>% filter(mp_h2 > 0, !is.na(mp_h2),cluster == "normal"), aes(x=mp_t2, y = mp_h2))
# ggsave("plots/inoc_5_normal_with2nd_peak.pdf")
# 
# # which not double strains have second peaks? 
# n_2nd <- para %>% filter(!cluster == "double", mp_h2 > 0, drytime ==0, inocl == 5)
# ddm_n2nd<- left_join(ts, n_2nd[,c("strain", "rep", "mp_h2","mp_t2")], by = c("strain", "rep"))
# ggplot(ddm_n2nd %>% filter(mp_h2 > 0, drytime ==0, inoc == 5, !cluster == "double"), aes(x=Time, y = value_J, group = interaction(strain,rep))) + 
#   geom_line(aes(col = factor(cluster))) + facet_wrap(~strain) + 
#   geom_point(data = ddm_n2nd %>% filter(mp_h2 > 0, !cluster == "double"), aes(x=mp_t2, y = mp_h2))
# ggsave("plots/inoc_5_not_double_with2nd_peak.pdf")
# 
# # Just inoc 5 / drytime == 0
# 
# g <- ggplot(cluster_data%>% filter(drytime ==0, inocl == 5), aes(x = cluster, y = value)) + geom_boxplot(aes(fill = name)) +  facet_wrap(~name, scales = "free") + 
#   stat_compare_means(comparisons = my_comparisons)
# ggsave("plots/inoc_5_clusters_boxplot.pdf")
