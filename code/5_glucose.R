#### Analyse data from one strain 11016 that has glucose concentration variation 

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
# OD
# 3 reps  7 inoc / 4 gluc
# cleanded in double_peaks_111016_dex_cleaning.R
ds <- read.csv("data/11016_glucose.csv")[,-1]

# remove un-needed columns
ds <- ds %>% dplyr::select(-c("glucose_conc","inoc_name","variable"))

# Smooth and Differentiate OD
data_ds <- ds %>% group_by(rep, glucose, inoc) %>%
  mutate(ma_value = rollapply(value, 5, mean,fill = NA),
         differ = c(0,diff(ma_value)),
         ma_differ = rollapply(differ, 5, mean,fill = NA)) %>%
  ungroup() #%>%
# filter(Time > 1, Time < 22)


### PLOTS - can comment out after first time
ggplot(data_ds, aes(x=Time, y = ma_value, group = interaction(rep, inoc,glucose))) +
  geom_line(aes(col = factor(inoc))) +
  facet_wrap(~glucose) +
  scale_y_continuous("OD") +
  scale_x_continuous("Time (h)") +
  scale_color_discrete("Inoculum", breaks = c(1,2,3,4,5,6,7), labels = c(str2expression("10^1"),str2expression("10^2"),str2expression("10^3"),
                                                                         str2expression("10^4"),str2expression("10^5"),
                                                                         str2expression("10^6"),str2expression("10^7")))
ggsave("plots/glucose_od_basic.pdf")

ggplot(data_ds, aes(x=Time, y = differ, group = interaction(rep, inoc, glucose))) +
  geom_line(aes(col = factor(inoc))) +
  facet_wrap(glucose~rep, ncol = 3) +
  scale_y_continuous("Difference in OD per time step")+
  scale_x_continuous("Time (h)") +
  scale_color_discrete("Inoculum", breaks = c(1,2,3,4,5,6,7), labels = c(str2expression("10^1"),str2expression("10^2"),str2expression("10^3"),
                                                                         str2expression("10^4"),str2expression("10^5"),
                                                                         str2expression("10^6"),str2expression("10^7")))
ggsave("plots/glucose_od_differ.pdf")

# Smoother differ
ggplot(data_ds, aes(x=Time, y = ma_differ, group = interaction(rep, inoc, glucose))) +
  geom_line(aes(col = factor(inoc))) +
  facet_wrap(glucose~rep, ncol = 3) +
  scale_y_continuous("Smoothed Difference in OD per time step")+
  scale_x_continuous("Time (h)") +
  scale_color_discrete("Inoculum", breaks = c(1,2,3,4,5,6,7), labels = c(str2expression("10^1"),str2expression("10^2"),str2expression("10^3"),
                                                                         str2expression("10^4"),str2expression("10^5"),
                                                                         str2expression("10^6"),str2expression("10^7")))
ggsave("plots/glucose_od_differentiatedsmoothed_byrep.pdf")

ggplot(data_ds, aes(x=Time, y = differ, group = interaction(rep, inoc, glucose))) +
  geom_line(aes(col = factor(inoc))) +
  facet_wrap(~glucose, ncol = 4) +
  scale_y_continuous("Difference in OD per time step")+
  scale_x_continuous("Time (h)") +
  scale_color_discrete("Inoculum", breaks = c(1,2,3,4,5,6,7), labels = c(str2expression("10^1"),str2expression("10^2"),str2expression("10^3"),
                                                                         str2expression("10^4"),str2expression("10^5"),
                                                                         str2expression("10^6"),str2expression("10^7")))
ggsave("plots/glucose_od_differentiated_byinoc.pdf")

ggplot(data_ds, aes(x=Time, y = ma_differ, group = interaction(rep, inoc,glucose))) +
  geom_line(aes(col = factor(inoc))) +
  facet_wrap(~glucose) +
  scale_y_continuous("Difference in OD per time step") +
  scale_x_continuous("Time (h)") +
  scale_color_discrete("Inoculum", breaks = c(1,2,3,4,5,6,7), labels = c(str2expression("10^1"),str2expression("10^2"),str2expression("10^3"),
                                                                         str2expression("10^4"),str2expression("10^5"),
                                                                         str2expression("10^6"),str2expression("10^7")))
ggsave("plots/glucose_od_basic.pdf")



###### Heat flow 
hf <- read.csv("data/glucose_HF.csv")[,-1]

# remove un-needed columns
hf <- hf %>% dplyr::select(-c("glucose_conc","inoc_name","variable","baseline"))

w1<-which(hf$rep == "x")
w2<-which(hf$rep == "xi")
w3<-which(hf$rep == "xii")
hf[w1,"rep"] <- 1
hf[w2,"rep"] <- 2
hf[w3,"rep"] <- 3

# Smooth and Differentiate OD
data_hf <- hf %>% group_by(rep, glucose, inoc) %>%
  mutate(ma_value = rollapply(value, 5, mean,fill = NA),
         cumm = cumsum(ifelse(is.na(ma_value), 0, ma_value))) %>%
  ungroup() #%>%
# filter(Time > 1, Time < 22)

### PLOTS - can comment out after first time
ggplot(data_hf, aes(x=Time, y = ma_value, group = interaction(rep, inoc,glucose))) +
  geom_line(aes(col = factor(inoc))) +
  facet_wrap(~glucose) +
  scale_y_continuous("Heat flow") +
  scale_x_continuous("Time (h)") +
  scale_color_discrete("Inoculum", breaks = c(1,2,3,4,5,6,7), labels = c(str2expression("10^1"),str2expression("10^2"),str2expression("10^3"),
                                                                         str2expression("10^4"),str2expression("10^5"),
                                                                         str2expression("10^6"),str2expression("10^7")))
ggsave("plots/glucose_hf_basic.pdf")

ggplot(data_hf, aes(x=Time, y = cumm, group = interaction(rep, inoc, glucose))) +
  geom_line(aes(col = factor(inoc))) +
  facet_wrap(glucose~rep, ncol = 3) +
  scale_y_continuous("Cumulative heatflow")+
  scale_x_continuous("Time (h)") +
  scale_color_discrete("Inoculum", breaks = c(1,2,3,4,5,6,7), labels = c(str2expression("10^1"),str2expression("10^2"),str2expression("10^3"),
                                                                         str2expression("10^4"),str2expression("10^5"),
                                                                         str2expression("10^6"),str2expression("10^7")))
ggsave("plots/glucose_hf_cumm.pdf")



ggplot(data_hf, aes(x=Time, y = cumm, group = interaction(rep, inoc, glucose))) +
  geom_line(aes(col = factor(inoc))) +
  facet_wrap(inoc~glucose, ncol = 4) +
  scale_y_continuous("Cumulative heatflow")+
  scale_x_continuous("Time (h)") +
  scale_color_discrete("Inoculum", breaks = c(1,2,3,4,5,6,7), labels = c(str2expression("10^1"),str2expression("10^2"),str2expression("10^3"),
                                                                         str2expression("10^4"),str2expression("10^5"),
                                                                         str2expression("10^6"),str2expression("10^7")))
ggsave("plots/glucose_hf_cumm_byinoc.pdf")



##### Build joint data set
data_hf_j <- data_hf %>% dplyr::select(Time, rep, glucose, inoc, ma_value) %>% rename(value = ma_value) %>% mutate(exp = "hf")
data_ds_j <- data_ds %>% dplyr::select(Time, rep, glucose, inoc, ma_differ) %>% rename(value = ma_differ) %>% mutate(exp = "od")
data_hf_j$rep <- as.numeric(data_hf_j$rep)

data_all <- full_join(data_hf_j, data_ds_j)

ggplot(data_all, aes(x=Time, y = value, group = interaction(rep, inoc, exp))) + 
  geom_line(aes(linetype = factor(exp), col = factor(inoc))) + 
  facet_wrap(exp~glucose, scales = "free", ncol = 4) + 
  scale_color_discrete("Inoculum") + 
  scale_linetype("Data\ntype")
ggsave("plots/glucose_all_raw.pdf", width = 12, height = 8)

#######################********** ALL ****************##############################################################################################################################################################
#### All inocula
############################################################################################################################################################################################################################
# remove all data before the time point that they all have which is the max of the minimum times to avoid odd completion curves prior to start and after
cutoff_time_dn = max(data_all %>% filter(!is.na(value)) %>% dplyr::group_by(glucose, rep, exp, inoc) %>% dplyr::summarise(minn = min(Time)) %>% ungroup() %>% dplyr::select(minn))
cutoff_time_up = min(data_all %>% filter(!is.na(value)) %>% dplyr::group_by(glucose, rep, exp, inoc) %>% dplyr::summarise(maxx = max(Time)) %>% ungroup() %>% dplyr::select(maxx))

data_od0 <- data_all %>% filter(Time > cutoff_time_dn, Time < cutoff_time_up) 

# Normalise
max_vals_norm0 <- data_od0 %>% filter(glucose == 0) %>% group_by(inoc, rep, exp )%>% summarise(max_norms = max(value, na.rm = TRUE))
max_vals_norm2p5 <- data_od0 %>% filter(glucose == 2.5) %>% group_by(inoc, rep, exp )%>% summarise(max_norms = max(value, na.rm = TRUE))

data_od <- left_join(data_od0, max_vals_norm0) %>% mutate(compara_norm = value / max_norms)
data_od2p5 <- left_join(data_od0, max_vals_norm2p5) %>% mutate(compara_norm = value / max_norms)

ggplot(data_od, aes(x=Time, y = compara_norm, group = interaction(rep, exp, glucose))) + 
  geom_line(aes(col = interaction(rep), lty = exp), lwd = 1) + 
  facet_grid(inoc ~glucose, scales = "free") + 
  scale_x_continuous("Time (h)") 
ggsave("plots/glucose_data_compare_norm0.pdf")

ggplot(data_od2p5, aes(x=Time, y = compara_norm, group = interaction(rep, exp, glucose))) + 
  geom_line(aes(col = interaction(rep), lty = exp), lwd = 1) + 
  facet_grid(inoc ~glucose, scales = "free") + 
  scale_x_continuous("Time (h)") 
ggsave("plots/glucose_data_compare_norm2p5.pdf")

# Subtract normalised data? Need to complete: measured at different time points

data_od_normd <- data_od %>% ungroup() %>% filter(Time > cutoff_time_dn, Time < cutoff_time_up) %>% # make sure cover same time for all 
  complete(rep, glucose, exp, inoc, Time) %>% dplyr::mutate(compara_norm_inp = na_ma(compara_norm, k = 4, weighting = "linear", maxgap = 10)) %>% # fill in all time points and then linear imputation between (tried exponential and simple but get more odd bumps)
  dplyr::group_by(inoc, glucose, rep, exp) %>% dplyr::select(Time, glucose,rep, exp, inoc, compara_norm_inp) %>% # Take imputed values
  pivot_wider(id_cols = c(glucose, Time, rep, inoc), names_from = exp, values_from = compara_norm_inp) %>%dplyr::mutate(nongrowth_only = hf - od) %>% # look for difference between OD and heat output
  pivot_longer(cols = c("hf","od"), names_to = "exp", values_to ="imput_val")

data_od_normd_ana <- left_join(data_od, data_od_normd) # not sure why need all this information

data_od_normd2p5 <- data_od2p5 %>% ungroup() %>% filter(Time > cutoff_time_dn, Time < cutoff_time_up) %>% # make sure cover same time for all 
  complete(rep, glucose, exp, inoc, Time) %>% dplyr::mutate(compara_norm_inp = na_ma(compara_norm, k = 4, weighting = "linear", maxgap = 10)) %>% # fill in all time points and then linear imputation between (tried exponential and simple but get more odd bumps)
  dplyr::group_by(inoc, glucose, rep, exp) %>% dplyr::select(Time, glucose,rep, exp, inoc, compara_norm_inp) %>% # Take imputed values
  pivot_wider(id_cols = c(glucose, Time, rep, inoc), names_from = exp, values_from = compara_norm_inp) %>%dplyr::mutate(nongrowth_only = hf - od)# %>% # look for difference between OD and heat output
#pivot_longer(cols = c("hf","od"), names_to = "exp", values_to ="imput_val")

data_od_normd_ana2p5 <- left_join(data_od2p5, data_od_normd2p5) # not sure why need all this information

### Merge these two controls 
data_od_normd2p5$control = 2.5
data_od_normd$control = 0 # USE this control 

data_od_normd_ana <- data_od_normd #rbind(data_od_normd2p5, data_od_normd)

# In ana: exp affects imput_val but not nongrowth
ggplot(data_od_normd_ana, aes(x=Time, y = nongrowth_only, group = interaction(rep, glucose))) + 
  geom_line(aes(col = interaction(rep)), lwd = 1) + 
  #geom_line(aes(y = imput_val,col = interaction(exp,rep))) + 
  facet_grid(inoc~glucose+control, scales = "free") + 
  scale_x_continuous("Time (h)") 
ggsave("plots/gluc_data_nongrowth_togplot.pdf")

g1 <- ggplot(data_od_normd_ana, aes(x=Time, y = nongrowth_only, group = interaction(rep, glucose))) + 
  geom_line(aes(col = interaction(rep))) + 
  facet_grid(inoc~glucose+control) + 
  scale_x_continuous("Time (h)") + 
  scale_y_continuous("Non growth only") + 
  scale_color_discrete("Replicate") + 
  geom_hline(yintercept = 0)

g2 <- ggplot(data_od_normd_ana, aes(x=Time, group = interaction(exp,rep, glucose))) + 
  geom_line(aes(y = imput_val,col = interaction(rep), linetype = exp)) + 
  facet_grid(inoc~glucose+control, scales = "free") + 
  scale_x_continuous("Time (h)") + 
  scale_y_continuous("Normalised measure") + 
  scale_color_discrete("Replicate")

g2 / g1 
ggsave("plots/gluc_data_nongrowth_tog_grid.pdf")

### Look at links to inocula
ggplot(data_od_normd_ana, aes(x=Time, y = nongrowth_only, group = interaction(glucose, inoc, rep, control))) + 
  geom_line(aes(col = factor(inoc))) + 
  facet_grid(rep~control+glucose) + 
  scale_x_continuous("Time (h)") + 
  scale_y_continuous("Non growth only") + 
  scale_color_discrete("Inocula") + 
  geom_hline(yintercept = 0)

#######################********** Extract characteristics ****************##############################################################################################################################################################
#### Extract key parameters from the non-growth data
############################################################################################################################################################################################################################

## What are the glucose levels?
u <- as.character(unique(data_od$glucose,stringsasFactors = FALSE)) # strains
## How many replicates? 
r <- unique(data_od$rep) # replicates
## How many inoc?
inn <- unique(data_od$inoc)

# Where the parameters for each strain are stored
param <- matrix(0, length(u)*length(r)*length(inn), 9); # number of strains x number of replicates x number of experimental conditions
index <- 1 # for counting 
max <- c() # for calibration

data_od <- data_od %>% dplyr::ungroup()

## Run thru each glucose/rep/drying time/ inoculum: fit a growth curve to each
for(jj in 1:length(u)){ # for each strain
  r <- unlist(unique(data_od %>% filter(glucose == u[jj]) %>% dplyr::select(rep))[,1]) # which replicates for this glucose (all have different names)
  if(length(r) < 3){print(paste0(u[jj], " has too few replicates"))}else{ # filter out if only 2 replicates
    for(ii in 1:length(r)){ # for each replicate
      for(pp in 1:length(inn)){ # for each inoc
        
        data1 <- data_od %>% filter(glucose == u[jj], rep == r[ii], inoc == inn[pp], exp == "od")
        data2 <- data_od %>% filter(glucose == u[jj], rep == r[ii], inoc == inn[pp], exp == "hf")
        data_nong <- data_od_normd_ana %>% filter(glucose == u[jj], rep == r[ii], inoc == inn[pp]) 
        
        if(dim(data1)[1] > 0){ # if this replicate exists for this glucose (i.e. there is data)
          
          print(c(jj, u[jj], r[ii], inn[pp],ex[kk])) # output so can track how it is working
          p <- charac_extract(data_nong, "Time", "nongrowth_only", data1, data2, "compara_norm", paste(glucose, replicate, condition, inocl,sep="_")) ### NEW function: simplification of more complex data recognition for clustering
          
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
#param <- param[-which(param[,3]==0),] # remove any empty rows (when replicates < 3)
param <- as.data.frame(param)
colnames(param) <- c("glucose","rep","inoc",
                     "t_max_h_flow", "v_max_h_flow", 
                     "t_min_h_flow", "v_min_h_flow", 
                     "auc","lagtime")
param$rep <- as.numeric(param$rep)
param$auc <- as.numeric(param$auc)
dim(param)

## Store so don't have to run above
write.csv(param, "output/simple_extract_glucose_allinoc.csv")
param <- read_csv("output/simple_extract_glucose_allinoc.csv")[,-1]

##########################################################################################
#### *********** Analyse above characteristics*************** ####################
##########################################################################################
## Take mean over replicates
param_mean <- param %>% group_by(glucose, inoc) %>% summarise(mean_t_max = mean(t_max_h_flow),mean_h_max = mean(v_max_h_flow),
                                                             mean_t_min = mean(t_min_h_flow),mean_h_min = mean(v_min_h_flow),
                                                             mean_auc = mean(auc))

# Heatmaps over inocula
g1 <- ggplot(param_mean, aes(inoc, factor(glucose), fill= mean_t_max)) + geom_tile() +
  scale_fill_distiller(palette = "RdBu","Time to\nmax.\nnon-growth", ) + 
  scale_x_continuous("Inoculum") + 
  scale_y_discrete("Glucose concentration") #+ 
#ggtitle("Time to max non-growth output")

g2 <- ggplot(param_mean, aes(inoc, factor(glucose), fill= mean_h_max)) + geom_tile() +
  scale_fill_distiller(palette = "RdBu","Height of\nmax.\nnon-growth", ) + 
  scale_x_continuous("Inoculum") + 
  scale_y_discrete("Glucose concentration") #+ 
#ggtitle("Height of max non-growth output")

g3 <- ggplot(param_mean, aes(inoc, factor(glucose), fill= mean_t_min)) + geom_tile() +
  scale_fill_distiller(palette = "RdBu","Time to\nmin.\nnon-growth", ) + 
  scale_x_continuous("Inoculum") + 
  scale_y_discrete("Glucose concentration") #+ 
#ggtitle("Time to min non-growth output")

g4 <- ggplot(param_mean, aes(inoc, factor(glucose), fill= mean_h_min)) + geom_tile() +
  scale_fill_distiller(palette = "RdBu","Height of \nmin.\nnon-growth", ) + 
  scale_x_continuous("Inoculum") + 
  scale_y_discrete("Glucose concentration") #+ 
#ggtitle("Height of min non-growth output")

g5 <- ggplot(param_mean, aes(inoc, factor(glucose), fill= mean_auc)) + geom_tile() +
  scale_fill_distiller(palette = "RdBu","AUC \nnon-growth", ) + 
  scale_x_continuous("Inoculum") + 
  scale_y_discrete("Glucose concentration") #+ 
#ggtitle("AUC non-growth output")

g1 + g3 + g5 + g2 + g4 + plot_layout(ncol = 3)
ggsave("plots/glucose_summary_heatmap.pdf", width = 12, height = 8)

#### Plot over glucose concentrations
param_long <- param %>% pivot_longer(cols = t_max_h_flow:auc)
param_long$name <- factor(param_long$name, levels = c("t_max_h_flow","t_min_h_flow" ,"auc", "v_max_h_flow","v_min_h_flow"))

ggplot(param_long, aes(x=inoc, y = value, group = interaction(glucose,rep))) + geom_line(aes(col = factor(glucose))) + 
  facet_wrap(~name, scales = "free") + 
  scale_color_discrete("Glucose concentration") + 
  scale_x_continuous("Inoculum")
ggsave("plots/glucose_summary_lines.pdf")

ggplot(param_long, aes(x=inoc, y = value, group = interaction(glucose))) + geom_point(aes(col = factor(glucose))) + 
  geom_smooth(aes(col = factor(glucose), fill = factor(glucose)),method='lm', formula= y~x) + 
  facet_wrap(~name, scales = "free") + 
  scale_color_discrete("Glucose concentration") + scale_fill_discrete("Glucose concentration") + 
  scale_x_continuous("Inoculum")
ggsave("plots/glucose_summary_smoothed.pdf")

