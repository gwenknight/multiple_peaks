### Analysis of growth data
# Extract all needed parameters for clustering

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
write.csv(ddm, "output/time_series_all.csv")

###******* MODEL FITTING *************#######################################################################################################################################

## What are the strains?
u <- as.character(unique(ddm$strain,stringsasFactors = FALSE)) # strains
## How many replicates? 
r <- unique(ddm$rep) # replicates
# What are the inoculums? 
q <- unique(ddm$inoc)

# Where the parameters for each strain are stored
param <- matrix(0, length(u)*length(r)*length(q)*3, 51); # number of strains x number of replicates x number of experimental conditions
index <- 1 # for counting 
max <- c() # for calibration
drying_times <- c(0,168)

## Run thru each strain/rep/drying time/ inoculum: fit a growth curve to each
for(jj in 1:length(u)){ # for each strain
  r <- unique(ddm %>% filter(strain == u[jj]) %>% dplyr::select(rep))[,1] # which replicates for this strain (all have different names)
  if(length(r) < 3){print(paste0(u[jj], " has too few replicates"))}else{ # filter out if only 2 replicates
    for(ii in 1:length(r)){ # for each replicate
      for(kk in c(1,2)){ #each of the experimental conditions
        for(ll in 1:length(q)){ #each of the inocula
          
          strain <- u[jj];
          replicate <- r[ii]
          condition <- drying_times[kk]
          inocl <- q[ll]
          
          wi <- intersect(which(ddm$strain == strain),which(ddm$rep == replicate)) # if fit to each replicate
          wj <- intersect(wi, which(ddm$drytime == condition))
          w <- intersect(wj, which(ddm$inoc == as.numeric(inocl)))
          
          if(length(w) > 0){ # if this replicate exists for this strain (i.e. there is data)
            data1 <- ddm[w,] # Grab data
            
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
write.csv(param, "output/param_multiple_peaks.csv")

