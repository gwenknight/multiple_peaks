#### Clustering of glucose dehydration data


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
data1 <- read_csv("data/glucose_HF_noBC_t8.csv")[,-1] # without baseline correction 
data2 <- read_csv("data/glucose_HF_t8.csv")[,-1] # with baseline correction 
data01 <- read_csv("data/glucose_HF.csv")[,-1] # time 0 data 
data02 <- read_csv("data/glucose_HF_noBC.csv")[,-1] # time 0 data without baseline correction

data01$drytime <- 0
data02$drytime <- 0
data1$drytime <- 7
data2$drytime <- 7

ddm0 <- rbind(data1, data02) # only WITHOUT baseline correction
w1 <- which(ddm0$rep == "x")
w2 <- which(ddm0$rep == "xi")
w3 <- which(ddm0$rep == "xii")
ddm0$rep <- as.numeric(ddm0$rep)
ddm0[w1,"rep"] <- 1
ddm0[w2,"rep"] <- 2
ddm0[w3,"rep"] <- 3

###*** Clean those that have v high low early data (remove for this one)
ttt <- intersect(intersect(intersect(intersect(which(ddm0$glucose == 1.25), which(ddm0$rep == 2)), which(ddm0$glucose_conc == "D")),which(ddm0$inoc == 2)),which(ddm0$drytime == 7))

ppp <- intersect(intersect(intersect(intersect(which(ddm0$glucose == 5), which(ddm0$rep == 1)), which(ddm0$glucose_conc == "B")),which(ddm0$inoc == 1)),which(ddm0$drytime == 7))
ddd <- intersect(intersect(intersect(intersect(which(ddm0$glucose == 5), which(ddm0$rep == 2)), which(ddm0$glucose_conc == "B")),which(ddm0$inoc == 1)),which(ddm0$drytime == 7))
qqq <- intersect(intersect(intersect(intersect(which(ddm0$glucose == 5), which(ddm0$rep == 3)), which(ddm0$glucose_conc == "B")),which(ddm0$inoc == 2)),which(ddm0$drytime == 7))
sss <- intersect(intersect(intersect(intersect(which(ddm0$glucose == 5), which(ddm0$rep == 1)), which(ddm0$glucose_conc == "B")),which(ddm0$inoc == 2)),which(ddm0$drytime == 7))

## Which removed? 
ddm0_rem <- ddm0[c(sss,ttt,ppp,qqq,ddd),]
g1 <- ggplot(ddm0_rem, aes(x=Time, y = value, group = interaction(variable, rep, inoc, glucose, drytime))) + geom_line(aes(col = variable)) + ggtitle("Removed") + 
  scale_y_continuous(lim = c(round_any(min(ddm0$value), 10, f = floor),round_any(max(ddm0$value), 10, f = ceiling)))
g2 <- ggplot(ddm0, aes(x=Time, y = value, group = interaction(variable, rep, inoc, glucose, drytime))) + geom_line(aes(col = variable)) + ggtitle("All strains") + 
  scale_y_continuous(lim = c(round_any(min(ddm0$value), 10, f = floor),round_any(max(ddm0$value), 10, f = ceiling)))
g1 + g2
ggsave("plots/glucose_conc_outliers.jpeg")

ddm <- ddm0[-c(sss,ttt,ppp,qqq,ddd),]

# Need to have two reps at least for clustering otherwise just an outlier: 
ddm %>% group_by(variable, glucose, inoc, drytime) %>% count(rep) %>% count() %>% filter(n < 3) %>% print(n=Inf)
# Two have only one replicate => have to remove as can't call cluster with only one dataset 
# glucose 5, drytime 7, inoculum size 1 + 2 
ddm <- ddm %>% filter(!variable %in% c("B7","B8"))


###******* MODEL FITTING *************#######################################################################################################################################

## What are the glucose concentrations?
u <- as.character(unique(ddm$glucose_conc,stringsasFactors = FALSE)) # strains
## How many replicates? 
r <- unique(ddm$rep) # replicates
# What are the inoculums? 
q <- unique(ddm$inoc)


# Where the parameters for each strain are stored
param <- matrix(0, length(u)*length(r)*length(q)*3, 51); # number of strains x number of replicates x number of experimental conditions
index <- 1 # for counting 
max <- c() # for calibration
drying_times <- unique(ddm$drytime)

## Run thru each strain/rep/drying time/ inoculum: fit a growth curve to each
for(jj in 1:length(u)){ # for each strain
  r <- unlist(unique(ddm %>% filter(glucose_conc == u[jj]) %>% dplyr::select(rep))[,1]) # which replicates for this strain (all have different names)
  if(length(r) < 3){print(paste0(u[jj], " has too few replicates"))}else{ # filter out if only 2 replicates
    for(ii in 1:length(r)){ # for each replicate
      for(kk in c(1,2)){ #each of the experimental conditions
        for(ll in 1:length(q)){ #each of the inocula
          
          conc <- u[jj];
          replicate <- unlist(r)[ii]
          condition <- drying_times[kk]
          inocl <- q[ll]
          
          wi <- intersect(which(ddm$glucose_conc == conc),which(ddm$rep == replicate)) # if fit to each replicate
          wj <- intersect(wi, which(ddm$drytime == condition))
          w <- intersect(wj, which(ddm$inoc == as.numeric(inocl)))
          
          if(length(w) > 0){ # if this replicate exists for this strain (i.e. there is data)
            data1 <- ddm[w,] # Grab data
            
            print(c(jj, conc, replicate, condition, inocl)) # output so can track how it is working
            p <- cut_extract_dp(data1, "Time", "value", paste(conc, replicate, condition, inocl,sep="_")) ### NEW function: runs on any timeseries: gives baseline and cut parameter analysis
            
            ## Required parameters
            
            param[index,] <- c(conc, replicate, condition, inocl, p$param)
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
colnames(param) <- c("strain","rep","drytime","inoc",
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
write.csv(param, "output/param_multiple_peaks_glucose_dehydration.csv")

########## ************************************************************************************************ #########################
##### Clustering #######
# For now have "Strain" for "Glucose concentration" 
ddm <- rename(ddm, "strain" = "glucose_conc")
ddm <- rename(ddm, "value_J" = "value")
c <- cluster(ddm, param,plot_where = "plots/glucose_conc_") 
########## ************************************************************************************************ #########################

### Store output of clustering
write.csv(c$parameters,"output/clustered_parameters_glucose_dehydration.csv")
write.csv(c$ts,"output/clustered_time_series_glucose_dehydration.csv")

## Read in if doing later
c <- c()
c$parameters <- read_csv("output/clustered_parameters_glucose_dehydration.csv")
c$ts <- read.csv("output/clustered_time_series_glucose_dehydration.csv")

### Glucose into parameter data 
c$parameters$glucose <- 0
c$parameters[which(c$parameters$strain == "B"),"glucose"] <- 5
c$parameters[which(c$parameters$strain == "C"),"glucose"] <- 2.5
c$parameters[which(c$parameters$strain == "D"),"glucose"] <- 1.25

### NA cluster exploration
ggplot(c$ts %>% filter(cluster == ""), aes(x=Time, y = value_J, group = interaction(glucose, drytime, rep, strain, inoc))) + 
  geom_line(aes(col = factor(strain), linetype = factor(drytime))) + 
  facet_grid(glucose ~ inoc) + 
  scale_x_continuous("Time") + 
  scale_color_discrete("Strain") + 
  scale_linetype_discrete("Drytime")
ggsave("plots/glucose_conc_nonclustering.jpeg")

c$ts %>% filter(cluster == "") %>% group_by(variable, drytime, inoc, glucose, strain, rep) %>% slice(1)
# Outliers - how cluster? 
c_new <- cluster(ddm %>% filter(variable == "E7", drytime == 0), param %>% filter(strain == "E", inoc == 2, drytime == 0),plot_where = "plots/glucose_conc_") 
c_new <- cluster(ddm %>% filter(variable == "C4", drytime == 7), param %>% filter(strain == "C", inoc == 5, drytime == 7),plot_where = "plots/glucose_conc_") 

# cn <- cut_extract_dp(ddm %>% filter(variable == "C4", drytime == 7, rep == 1), "Time", "value_J", paste(conc, replicate, condition, inocl,sep="_")) ### NEW function: runs on any timeseries: gives baseline and cut parameter analysis
# dcn <- as.data.frame(t(cn$param))
# dcn <- as.numeric(dcn)
# colnames(dcn) <- c("t_m_h_flow", "v_m_h_flow", "exp_gr","lag","auc",
#                   "odd_peaks","odd_width","width_peak","odd_shoulder","odd_shoulder_past","odd_double",
#                   "shoulder_point_t","shoulder_point_v", "shoulder_point_past_t","shoulder_point_past_v",
#                   "cut_exp", "timepeak", "valpeak",
#                   "mp_t1","mp_t2","mp_t3","mp_t4","mp_t5","mp_t6","mp_t7","mp_t8","mp_t9","mp_t10",
#                   "mp_h1","mp_h2","mp_h3","mp_h4","mp_h5","mp_h6","mp_h7","mp_h8","mp_h9","mp_h10",
#                   "gap1","gap2","gap3","gap4","gap5","gap6","gap7","gap8","gap9")
# dcn$shoulder_point_past_t <- as.numeric(dcn$shoulder_point_past_t)
# dcn$shoulder_point_past_v <- as.numeric(dcn$shoulder_point_past_v)
# ggplot(ddm %>% filter(variable == "C4", drytime == 7, rep == 1), aes(x=Time, y=value_J)) + geom_line() + 
#   geom_point(data = dcn, aes(x = shoulder_point_past_t, y = shoulder_point_past_v))

### Overview by glucose
ggplot(c$ts, aes(x=Time, y = value_J, group = interaction(strain, inoc, rep, cluster))) + 
  geom_line(aes(colour = cluster)) + 
  facet_grid(inoc + drytime~glucose)

ggplot(c$parameters, aes(x=inoc, group = cluster)) + 
  geom_bar(stat = "count", position = "stack", aes(fill = cluster)) + 
  facet_grid(~drytime) + 
  scale_x_continuous("Inoculum") + 
  scale_y_continuous("Number of datasets") + 
  scale_fill_discrete("Cluster type")
ggsave("plots/glucose_cluster_drytime.jpeg")


### STORY PIC
para <- c$parameters %>% ungroup() %>% mutate(second_peak_h = ifelse((t_m_h_flow > (timepeak + 2)), v_m_h_flow, ifelse(mp_t2 > (timepeak + 4), mp_h2, 0))) 

cluster_data <- para %>% dplyr::select(c("cluster", "inoc","drytime","valpeak", "timepeak", "auc", "exp_gr", "shoulder_point_v","shoulder_point_t",
                                         "shoulder_point_past_v","shoulder_point_past_t", "second_peak_h")) %>% 
  pivot_longer(cols = c(valpeak:second_peak_h)) 

my_comparisons <- list( c("post_shoulder", "double"), c("double", "normal"), c("normal", "post_shoulder"),c("normal","spike"),c("spike","wide"),
                        c("double","spike"),c("double","wide"))

cluster_data$cluster <- factor(cluster_data$cluster, levels = c("normal","double","spike","post_shoulder","wide","unclustered"))
var_name <- c(
  auc = "AUC",
  valpeak = "Max value 1st peak",
  exp_gr = "Max exp growth rate",
  second_peak_h = "Max value 2nd peak"
)

ggplot(cluster_data %>% filter(inoc == 5) %>% filter(name %in% c("auc", "valpeak","exp_gr","second_peak_h")) %>% 
         mutate(across(name, factor, levels=c("auc", "valpeak","exp_gr","second_peak_h"))), aes(x = cluster, y = value)) + 
  geom_boxplot(aes(fill = cluster)) +  facet_wrap(drytime~name,ncol = 4, scales = "free", labeller = labeller(name = var_name))
ggsave("plots/glucose_inoc_5_story_boxplot.jpeg", width = 20, height = 8)

ggplot(cluster_data %>% filter(name %in% c("auc", "valpeak","exp_gr","second_peak_h")) %>% 
         mutate(across(name, factor, levels=c("auc", "valpeak","exp_gr","second_peak_h"))), aes(x = interaction(inoc,cluster), y = value)) + 
  geom_boxplot(aes(fill = cluster)) +  facet_wrap(drytime~name,ncol = 4, scales = "free", labeller = labeller(name = var_name))




 
