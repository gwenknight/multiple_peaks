#### Analyse data from one strain 11016 that has glucose concentration variation and 7 day dehydration 

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
library(zoo)
library(patchwork)
library(imputeTS)
library(MESS) # for auc
theme_set(theme_bw(base_size=14)) # theme setting for plots: black and white (bw) and font size (24)
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(14)

setwd(here::here())

# 11016 data
# Only HF
# 3 reps  7 inoc / 4 gluc
data1 <- read_csv("data/glucose_HF_noBC_t8.csv")[,-1] # with baseline correction 
data2 <- read_csv("data/glucose_HF_t8.csv")[,-1] # without baseline correction 
data01 <- read_csv("data/glucose_HF.csv")[,-1] # time 0 data 
data02 <- read_csv("data/glucose_HF_noBC.csv")[,-1] # time 0 data without baseline correction

data01$drytime <- 0
data02$drytime <- 0
data1$drytime <- 7
data2$drytime <- 7

data_dry <- data1 # with correction (automatic from calscreener)
data_dry$value_base <- data2$value # (without correction)

### Raw data
ggplot(data_dry, aes(x= Time, y = value, group = interaction(rep, glucose_conc, inoc_name))) + geom_line(aes(colour = factor(inoc))) + 
  facet_wrap(~glucose) + 
  scale_y_continuous("Heat flow") + 
  scale_color_discrete("Inoculum")
ggsave("plots/glucose_drying_raw_data.png")

# baseline correction looks better? 
ggplot(data_dry, aes(x= Time, y = value_base, group = interaction(rep, glucose_conc, inoc_name))) + geom_line(aes(colour = factor(inoc))) + 
  facet_wrap(~glucose) + 
  scale_y_continuous("Heat flow") + 
  scale_color_discrete("Inoculum")
ggsave("plots/glucose_drying_raw_data_base_correct.png")



#######################********** Extract characteristics ****************##############################################################################################################################################################
#### Extract key parameters from the HF data after 7 days drying (NO OD data so can't get non-growth)
############################################################################################################################################################################################################################
data_hf <- rbind(data1, data2, data01, data02) # on top of each other instead of alongside

## What are the glucose levels?
u <- as.character(unique(data_hf$glucose,stringsasFactors = FALSE)) # strains
## How many replicates? 
r <- unique(data_hf$rep) # replicates
## How many inoc?
inn <- unique(data_hf$inoc)
## Baseline correction yes / no
bc <- unique(data_hf$baseline)
## drytimes
dry <- unique(data_hf$drytime)

# Where the parameters for each strain are stored
param <- matrix(0, length(u)*length(r)*length(inn)*length(bc)*length(dry), 11); # number of strains x number of replicates x number of experimental conditions
index <- 1 # for counting 
max <- c() # for calibration

data_hf <- data_hf %>% dplyr::ungroup()

## Run thru each glucose/rep/drying time/ inoculum: fit a growth curve to each
for(jj in 1:length(u)){ # for each strain
  r <- unlist(unique(data_hf %>% filter(glucose == u[jj]) %>% dplyr::select(rep))[,1]) # which replicates for this glucose (all have different names)
  if(length(r) < 3){print(paste0(u[jj], " has too few replicates"))}else{ # filter out if only 2 replicates
    for(ii in 1:length(r)){ # for each replicate
      for(pp in 1:length(inn)){ # for each inoc
        for(qq in 1:length(bc)){ # for each baseline
          for(dp in 1:length(dry)){ # for each drytime
            
            data1 <- data_hf %>% filter(glucose == u[jj], rep == r[ii], inoc == inn[pp], baseline == bc[qq], drytime == dry[dp])
            data2 <- data1 
            data2$value <- 0 # set to zero so AUC = 0 => subtraction gives data 1 AUC in below
            
            if(dim(data1)[1] > 0){ # if this replicate exists for this glucose (i.e. there is data)
              
              print(c(jj, u[jj], r[ii], inn[pp],bc[qq],dry[dp])) # output so can track how it is working
              p <- charac_extract(data1, "Time", "value", data1, data2, "value", paste(u[jj], r[ii], inn[pp], bc[qq],dry[dp],sep="_")) ### NEW function: simplification of more complex data recognition for clustering
              
              ## Required parameters
              param[index,] <- c(u[jj], r[ii], inn[pp],bc[qq],dry[dp], p)
              index <- index + 1 # counting for storing matrix - next row each iteration
            }
          }
        }
      }
    }
  }
}

# The original set 
param_orig <- param

## Fitted parameters
#param <- param[-which(param[,3]==0),] # remove any empty rows (when replicates < 3)
param <- as.data.frame(param)
colnames(param) <- c("glucose","rep","inoc","baseline","drytime",
                     "t_max_h_flow", "v_max_h_flow", 
                     "t_min_h_flow", "v_min_h_flow", 
                     "auc","lagtime")
#param$rep <- as.numeric(param$rep)
param$auc <- as.numeric(param$auc)
dim(param)

## Store so don't have to run above
write.csv(param, "output/simple_extract_glucose_allinoc_hf.csv")
param <- read_csv("output/simple_extract_glucose_allinoc_hf.csv")[,-1]

##########################################################################################
#### *********** Analyse above characteristics*************** ####################
##########################################################################################
## Take mean over replicates
param_mean_bothbl <- param %>% group_by(glucose, inoc, baseline, drytime) %>% dplyr::summarise(mean_t_max = mean(t_max_h_flow),mean_h_max = mean(v_max_h_flow),
                                                                                        mean_t_min = mean(t_min_h_flow),mean_h_min = mean(v_min_h_flow),
                                                                                        mean_auc = mean(auc))

#### Two baselines
for(i in c("yes", "no")){
  for(j in c(0,7)){
    param_mean = param_mean_bothbl %>% filter(baseline == i, drytime == j)
    
    # Heatmaps over inocula
    g1 <- ggplot(param_mean, aes(inoc, factor(glucose), fill= mean_t_max)) + geom_tile() +
      scale_fill_distiller(palette = "RdBu","Time to\nmax", ) + 
      scale_x_continuous("Inoculum") + 
      scale_y_discrete("Glucose concentration") #+ 
    #ggtitle("Time to max non-growth output")
    
    g2 <- ggplot(param_mean, aes(inoc, factor(glucose), fill= mean_h_max)) + geom_tile() +
      scale_fill_distiller(palette = "RdBu","Height of\nmax", ) + 
      scale_x_continuous("Inoculum") + 
      scale_y_discrete("Glucose concentration") #+ 
    #ggtitle("Height of max non-growth output")
    
    g3 <- ggplot(param_mean, aes(inoc, factor(glucose), fill= mean_t_min)) + geom_tile() +
      scale_fill_distiller(palette = "RdBu","Time to\nmin", ) + 
      scale_x_continuous("Inoculum") + 
      scale_y_discrete("Glucose concentration") #+ 
    #ggtitle("Time to min non-growth output")
    
    g4 <- ggplot(param_mean, aes(inoc, factor(glucose), fill= mean_h_min)) + geom_tile() +
      scale_fill_distiller(palette = "RdBu","Height of \nmin", ) + 
      scale_x_continuous("Inoculum") + 
      scale_y_discrete("Glucose concentration") #+ 
    #ggtitle("Height of min non-growth output")
    
    g5 <- ggplot(param_mean, aes(inoc, factor(glucose), fill= mean_auc)) + geom_tile() +
      scale_fill_distiller(palette = "RdBu","AUC ", ) + 
      scale_x_continuous("Inoculum") + 
      scale_y_discrete("Glucose concentration") #+ 
    #ggtitle("AUC non-growth output")
    
    g1 + g5 + g2 + plot_layout(ncol = 1) + plot_annotation(title = paste0("Baseline correction: ", i, " Drytime: ",j))
    ggsave(paste0("plots/glucose_summary_heatmap_hf_",i,j,".png"), width = 5, height = 8)
    
    #### Plot over glucose concentrations
    param_long <- param %>% pivot_longer(cols = t_max_h_flow:auc) %>% filter(baseline == i) %>% filter(name %in% c("t_max_h_flow" ,"auc", "v_max_h_flow"))
    param_long$name <- factor(param_long$name, levels = c("t_max_h_flow" ,"auc", "v_max_h_flow"))
    
    ggplot(param_long, aes(x=inoc, y = value, group = interaction(glucose,rep))) + geom_line(aes(col = factor(glucose))) + 
      facet_wrap(~name, scales = "free") + 
      scale_color_discrete("Glucose concentration") + 
      scale_x_continuous("Inoculum") + 
      ggtitle(paste0("Baseline correction: ", i, " Drytime: ",j))
    ggsave(paste0("plots/glucose_summary_lines_hf_",i,j,".png"))
    
    ggplot(param_long, aes(x=inoc, y = value, group = interaction(glucose))) + geom_point(aes(col = factor(glucose))) + 
      geom_smooth(aes(col = factor(glucose), fill = factor(glucose)),method='lm', formula= y~x) + 
      facet_wrap(~name, scales = "free") + 
      scale_color_discrete("Glucose concentration") + scale_fill_discrete("Glucose concentration") + 
      scale_x_continuous("Inoculum") + 
      ggtitle(paste0("Baseline correction: ", i, " Drytime: ",j))
    ggsave(paste0("plots/glucose_summary_smoothed_hf_",i,j,".png"))
  }
}