#### First analysis / look at sugar data for 11016

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
  facet_wrap(exp~glucose, scales = "free", ncol = 4)
ggsave("plots/glucose_all_raw.pdf")

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


########## MAX NON-GROWTH

# Filter by HF as don't care about input values - same nongrowth for hf or od
max_nongrowth_tab <- data_od_normd_ana %>% filter(exp == "hf") %>% group_by(rep, inoc, glucose) %>% mutate(max_non = max(nongrowth_only)) %>%  
  filter(nongrowth_only == max_non)

g1 <- ggplot(max_nongrowth_tab, aes(x=inoc, y = nongrowth_only,group = interaction(glucose, rep))) + geom_point(aes(col = factor(glucose))) + 
  geom_line(aes(col = factor(glucose))) + 
  scale_y_continuous("Maximum non-growth") + scale_x_continuous("Inoculum") +
  scale_color_discrete("Glucose")
ggsave("plots/glucose_inoc_vs_max_nongrowth_.pdf")

g2 <- ggplot(max_nongrowth_tab, aes(x=inoc, y = Time, group = interaction(glucose, rep))) + geom_point(aes(col = factor(glucose))) +  
  geom_line(aes(col = factor(glucose))) + 
  scale_y_continuous("Time to maximum non-growth") + scale_x_continuous("Inoculum") + 
  scale_color_discrete("Glucose")
ggsave("plots/glucose_inoc_vs_time_max_nongrowth.pdf")

g1 + g2 + plot_layout(guides = "collect")
ggsave("plots/glucose_inoc_vs_time&val_max_nongrowth.pdf", width = 10, height = 5)

g1 <- ggplot(max_nongrowth_tab, aes(x=inoc, y = nongrowth_only,group = interaction(glucose))) + geom_point(aes(col = factor(glucose))) + 
  geom_smooth(aes(col = factor(glucose), fill = factor(glucose)),method='lm', formula= y~x) + 
  scale_y_continuous("Maximum non-growth") + scale_x_continuous("Inoculum") + 
  facet_wrap(~glucose, ncol =  4)+ 
  scale_color_discrete("Glucose")+ 
  scale_fill_discrete("Glucose")
ggsave("plots/glucose_inoc_vs_max_nongrowth_smooth.pdf")

g2 <- ggplot(max_nongrowth_tab, aes(x=inoc, y = Time, group = interaction(glucose))) + geom_point(aes(col = factor(glucose))) +  
  geom_smooth(aes(col = factor(glucose), fill = factor(glucose)),method='lm', formula= y~x) + 
  scale_y_continuous("Time to maximum non-growth") + scale_x_continuous("Inoculum") + 
  facet_wrap(~glucose, ncol =  4)+ 
  scale_color_discrete("Glucose")+ 
  scale_fill_discrete("Glucose")
ggsave("plots/glucose_inoc_vs_time_max_nongrowth_smooth.pdf")

g1 / g2 + plot_layout(guides = "collect")
ggsave("plots/glucose_inoc_vs_time&val_max_nongrowth_smooth.pdf", width = 10, height = 5)

#### Heat map 
# Take the mean over the replicates
max_nongrowth_tab_mean <- max_nongrowth_tab %>% dplyr::select(glucose, inoc, rep, Time, max_non) %>% group_by(glucose, inoc) %>% 
  summarise(mean_time = mean(Time), mean_max = mean(max_non))

max_nongrowth_tab_mean$glucose <- as.character(max_nongrowth_tab_mean$glucose)
g1 <- ggplot(max_nongrowth_tab_mean, aes(inoc, glucose, fill= mean_time)) + geom_tile() +
  scale_fill_distiller(palette = "RdBu","Time to\nmax. non-growth", ) + 
  scale_x_continuous("Inoculum") + 
  scale_y_discrete("Glucose concentration")

g2 <- ggplot(max_nongrowth_tab_mean, aes(inoc, glucose, fill= mean_max)) + geom_tile() +
  scale_fill_distiller(palette = "PRGn","Value of\nmax. non-growth") + 
  scale_x_continuous("Inoculum") + 
  scale_y_discrete("Glucose concentration")

g1 + g2 
ggsave("plots/glucose_inoc_vs_time&val_max_heatmap.pdf", width = 15, height = 5)

#### MIN NON_GROWTH

# Filter by HF as don't care about input values - same nongrowth for hf or od
min_nongrowth_tab <- data_od_normd_ana %>% filter(exp == "hf") %>% group_by(rep, inoc, glucose) %>% mutate(min_non = min(nongrowth_only)) %>%  
  filter(nongrowth_only == min_non)

g1 <- ggplot(min_nongrowth_tab, aes(x=inoc, y = nongrowth_only,group = interaction(glucose, rep))) + geom_point(aes(col = factor(glucose))) + 
  geom_line(aes(col = factor(glucose))) + 
  scale_y_continuous("Minimum non-growth") + scale_x_continuous("Inoculum") +
  scale_color_discrete("Glucose")
ggsave("plots/glucose_inoc_vs_min_nongrowth_.pdf")

g2 <- ggplot(min_nongrowth_tab, aes(x=inoc, y = Time, group = interaction(glucose, rep))) + geom_point(aes(col = factor(glucose))) +  
  geom_line(aes(col = factor(glucose))) + 
  scale_y_continuous("Time to minimum non-growth") + scale_x_continuous("Inoculum") + 
  scale_color_discrete("Glucose")
ggsave("plots/glucose_inoc_vs_time_min_nongrowth.pdf")

g1 + g2 + plot_layout(guides = "collect")
ggsave("plots/glucose_inoc_vs_time&val_min_nongrowth.pdf", width = 10, height = 5)

g1 <- ggplot(min_nongrowth_tab, aes(x=inoc, y = nongrowth_only,group = interaction(glucose))) + geom_point(aes(col = factor(glucose))) + 
  geom_smooth(aes(col = factor(glucose), fill = factor(glucose)),method='lm', formula= y~x) + 
  scale_y_continuous("Minimum non-growth") + scale_x_continuous("Inoculum") + 
  facet_wrap(~glucose, ncol =  4)+ 
  scale_color_discrete("Glucose")+ 
  scale_fill_discrete("Glucose")
ggsave("plots/glucose_inoc_vs_min_nongrowth_smooth.pdf")

g2 <- ggplot(min_nongrowth_tab, aes(x=inoc, y = Time, group = interaction(glucose))) + geom_point(aes(col = factor(glucose))) +  
  geom_smooth(aes(col = factor(glucose), fill = factor(glucose)),method='lm', formula= y~x) + 
  scale_y_continuous("Time to minimum non-growth") + scale_x_continuous("Inoculum") + 
  facet_wrap(~glucose, ncol =  4)+ 
  scale_color_discrete("Glucose")+ 
  scale_fill_discrete("Glucose")
ggsave("plots/glucose_inoc_vs_time_min_nongrowth_smooth.pdf")

g1 / g2 + plot_layout(guides = "collect")
ggsave("plots/glucose_inoc_vs_time&val_min_nongrowth_smooth.pdf", width = 10, height = 5)

#### Heat map 
# Take the mean over the replicates
min_nongrowth_tab_mean <- min_nongrowth_tab %>% dplyr::select(glucose, inoc, rep, Time, min_non) %>% group_by(glucose, inoc) %>% 
  summarise(mean_time = mean(Time), mean_min = mean(min_non))

min_nongrowth_tab_mean$glucose <- as.character(min_nongrowth_tab_mean$glucose)
g1 <- ggplot(min_nongrowth_tab_mean, aes(inoc, glucose, fill= mean_time)) + geom_tile() +
  scale_fill_distiller(palette = "RdBu","Time to\nmin. non-growth", ) + 
  scale_x_continuous("Inoculum") + 
  scale_y_discrete("Glucose concentration")

g2 <- ggplot(min_nongrowth_tab_mean, aes(inoc, glucose, fill= mean_min)) + geom_tile() +
  scale_fill_distiller(palette = "PRGn","Value of\nmin. non-growth") + 
  scale_x_continuous("Inoculum") + 
  scale_y_discrete("Glucose concentration")

g1 + g2 
ggsave("plots/glucose_inoc_vs_time&val_min_heatmap.pdf", width = 15, height = 5)



#### Extract AUC 
## What are the glucose concentrations?
u <- unique(data_od_normd_ana$glucose)
## How many replicates? 
r <- unique(data_od_normd_ana$rep) # replicates
# What are the inoculums? 
q <- unique(data_od_normd_ana$inoc)

# Where the parameters for each strain are stored
param <- matrix(0, length(u)*length(q)*length(r), 9); 
index <- 1 # for counting 


## Run thru each 
for(jj in 1:length(u)){ # for each glucose
  for(ii in 1:length(r)){ # for each replicate
    for(ll in 1:length(q)){ #each of the inocula
      
      data <- data_od_normd_ana %>% filter(exp == "hf") %>%
        filter(glucose == u[jj], rep == r[ii], inoc == q[ll]) %>% filter(Time > 3, Time < 18)
      data1 <- data_od_normd_ana %>% filter(exp == "hf") %>%
        filter(glucose == u[jj], rep == r[ii], inoc == q[ll]) %>% filter(Time > 3, Time < 15) # BEST ONE 
      data2 <- data_od_normd_ana %>% filter(exp == "hf") %>%
        filter(glucose == u[jj], rep == r[ii], inoc == q[ll]) %>% filter(Time > 5, Time < 22)
      
      
      if(dim(data)[1] > 0){ # if this replicate exists for this strain (i.e. there is data)
        
        print(c(u[jj], r[ii], q[ll])) # output so can track how it is working
        p <-0; p_all <- 0;p1 <-0; p_all1 <- 0;p2 <-0; p_all2 <- 0;
        try(p <- auc(data$Time, data$nongrowth_only, type = "spline", subdivisions = 100000))
        try(p_all <- auc(data$Time, data$nongrowth_only, type = "spline", subdivisions = 10000,absolutearea = TRUE))
        try(p1 <- auc(data1$Time, data1$nongrowth_only, type = "spline", subdivisions = 100000))
        try(p_all1 <- auc(data1$Time, data1$nongrowth_only, type = "spline", subdivisions = 10000,absolutearea = TRUE))
        try(p2 <- auc(data2$Time, data2$nongrowth_only, type = "spline", subdivisions = 100000))
        try(p_all2 <- auc(data2$Time, data2$nongrowth_only, type = "spline", subdivisions = 10000,absolutearea = TRUE))
        
        param[index,] <- c(u[jj], r[ii], q[ll],p, p_all,p1, p_all1,p2, p_all2)
        index <- index + 1 # counting for storing matrix - next row each iteration
        
      }
    }
  }
}


# The original set 
param_orig <- param

## Fitted parameters
param <- as.data.frame(param)
colnames(param) <- c("glucose","rep","inoc",
                     "auc1", "auc1_all","auc2", "auc2_all","auc3", "auc3_all")
which(param$auc1 == 0)
which(param$auc2 == 0)
which(param$auc3 == 0)
which(param$auc1_all == 0)
which(param$auc2_all == 0)
which(param$auc3_all == 0)

g1 <- ggplot(param, aes(x=inoc, y = auc2, group = interaction(rep,glucose),col = glucose)) + geom_point() + geom_line(aes()) + 
  scale_x_continuous("Inoculum") + 
  scale_y_continuous("AUC")
ggsave("plots/glucose_auc.pdf", width = 10, height = 5)

g2 <- ggplot(param, aes(x=inoc, y = auc2_all, group = interaction(rep,glucose),col = glucose)) + geom_point() + geom_line(aes()) + 
  scale_x_continuous("Inoculum") + 
  scale_y_continuous("AUC absoulte")
ggsave("plots/glucose_abs_auc.pdf", width = 10, height = 5)

g1 + g2 + plot_layout(guides = "collect")
ggsave("plots/glucose_auc_both_raw.pdf", width = 15, height = 5)

g1 <- ggplot(param, aes(x=inoc, y = auc2,group = interaction(glucose))) + geom_point(aes(col = factor(glucose))) + 
  geom_smooth(aes(col = factor(glucose), fill = factor(glucose)),method='lm', formula= y~x) + 
  scale_y_continuous("AUC") + scale_x_continuous("Inoculum") + 
  facet_wrap(~glucose, ncol =  4)+ 
  scale_color_discrete("Glucose")+ 
  scale_fill_discrete("Glucose")

g2 <- ggplot(param, aes(x=inoc, y = auc2_all,group = interaction(glucose))) + geom_point(aes(col = factor(glucose))) + 
  geom_smooth(aes(col = factor(glucose), fill = factor(glucose)),method='lm', formula= y~x) + 
  scale_y_continuous("AUC absolute") + scale_x_continuous("Inoculum") + 
  facet_wrap(~glucose, ncol =  4)+ 
  scale_color_discrete("Glucose")+ 
  scale_fill_discrete("Glucose")

g1 + g2 + plot_layout(guides = "collect")
ggsave("plots/glucose_auc_both_fit.pdf", width = 15, height = 5)


## Heat map 
param_mean <- param %>% group_by(glucose, inoc) %>% 
  summarise(mean_auc2 = mean(auc2), mean_auc2_all = mean(auc2_all))

param_mean$glucose <- as.character(param_mean$glucose)
g1 <- ggplot(param_mean, aes(inoc, glucose, fill= mean_auc2)) + geom_tile() +
  scale_fill_distiller(palette = "PRGn","AUC") + 
  scale_x_continuous("Inoculum") + 
  scale_y_discrete("Glucose concentration")

g2 <- ggplot(param_mean, aes(inoc, glucose, fill= mean_auc2_all)) + geom_tile() +
  scale_fill_distiller(palette = "PRGn","AUC\nabsolute") + 
  scale_x_continuous("Inoculum") + 
  scale_y_discrete("Glucose concentration")

g1 + g2 
ggsave("plots/glucose_auc_abs_heatmap.pdf", width = 15, height = 5)

