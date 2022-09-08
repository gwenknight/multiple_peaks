### Exploring OD data

library(zoo)
library(imputeTS) # to give weighter moving average - exponential weighting of those further away
library(patchwork)
library(tidyverse)
setwd(here::here())
theme_set(theme_bw(base_size = 11))

# Data to explore
data_od_orig <- read_csv("data/growth_ODvsCS_20220224.csv")[,-1]

# remove all data before the time point that they all have which is the max of the minimum times to avoid odd completion curves prior to start and after
cutoff_time_dn = max(data_od_orig %>% dplyr::group_by(rep, exp, inoc) %>% dplyr::summarise(minn = min(Time)) %>% ungroup() %>% dplyr::select(minn))
cutoff_time_up = min(data_od_orig %>% dplyr::group_by(rep, exp, inoc) %>% dplyr::summarise(maxx = max(Time)) %>% ungroup() %>% dplyr::select(maxx))


## Look at OD & CS data
ggplot(data_od_orig, aes(x=Time, y = value, group = interaction(inoc, exp, strain, rep))) + geom_line(aes(col = factor(inoc))) + 
  facet_grid(exp~strain, scales = "free")
ggsave("plots/ODvsC2_rawdata.pdf")

#######################********** 10^5 ****************##############################################################################################################################################################
#### JUST DO FOR 10^5 to start! 
############################################################################################################################################################################################################################

ggplot(data_od_orig %>% filter(inoc == 5), aes(x=Time, y = value, group = interaction(inoc, exp, strain, rep))) + geom_line(aes(col = factor(rep))) + 
  facet_grid(exp~strain, scales = "free")
ggsave("plots/ODvsC2_inoc5_rawdata.pdf")

data_od <- data_od_orig %>% filter(Time > cutoff_time_dn, Time < cutoff_time_up) %>% 
  dplyr::select(Time, rep, exp, value, strain, inoc) %>% 
  filter(inoc == 5) %>% # only 10^5 for this analysis
  group_by(rep, exp, strain) %>%
  mutate(ma_value = zoo::rollapply(value, 5, mean,fill = NA),
         differ = c(0,diff(ma_value)),
         compara = ifelse(exp == "CS", value, differ)) %>%
  ungroup() 

g2 <- ggplot(data_od, aes(x=Time, y = ma_value, group = interaction(rep, exp, strain))) + 
  geom_line(aes(col = factor(rep))) + 
  facet_grid(exp~strain,scales = "free") + 
  scale_x_continuous("Time (h)") 
ggsave("plots/ODvsCS_inoc5_data_smoothed.pdf")

g3a <- ggplot(data_od, aes(x=Time, y = compara, group = interaction(rep, exp, strain))) + 
  geom_line(aes(col = factor(rep))) + 
  facet_grid(exp~strain, scales = "free") + 
  scale_x_continuous("Time (h)") 
ggsave("plots/ODvsCS_inoc5_data_compare.pdf")

# Try to plot together? Need to normalise... 
## First tried scale by max. Normalise - but then don't get comparison between strains
# data_od <- data_od %>% group_by(strain, rep, exp) %>% mutate(max_v = max(compara, na.rm = TRUE),
#                                                              compara_norm = compara / max_v)
# Second find normal - divide by max of that then subtract
strains <- c("11016","11051", "11210", "11257")
clusters <- c("spike", "double", "spike", "normal")
max_vals_norm <- data_od %>% filter(strain == "11257") %>% group_by(rep, exp )%>% summarise(max_norms = max(compara, na.rm = TRUE))

data_od <- left_join(data_od, max_vals_norm) %>% mutate(compara_norm = compara / max_norms)


ggplot(data_od, aes(x=Time, y = compara_norm, group = interaction(rep, exp, strain))) + 
  geom_line(aes(col = factor(rep), lty = exp), lwd = 1) + 
  facet_wrap(~strain, scales = "free", ncol = 2) + 
  scale_x_continuous("Time (h)") + 
  scale_color_discrete("Replicate") + 
  scale_linetype("Data")
ggsave("plots/ODvsCS_inoc5_data_compare_norm.pdf")

# Subtract normalised data? Need to complete: measured at different time points

#data_play <- data_od[c(1:100,1149:1259),] %>% filter(Time > cutoff_time)
#data_play <- data_play %>% ungroup() %>% complete(exp, Time) %>% mutate(compara_norm_inp = na_ma(compara_norm, k = 2, weighting = "exponential", maxgap = Inf)) %>% print(n=Inf)
#ggplot(data_play, aes(x=Time, y = compara_norm, group = exp)) + geom_line(aes(group = exp, col = exp)) + 
#  geom_line(aes(y = compara_norm_inp, group = exp, col = exp), linetype = 2) 

data_od_normd <- data_od %>% ungroup() %>% filter(Time > cutoff_time_dn, Time < cutoff_time_up) %>% # make sure cover same time for all 
  complete(rep, strain, exp, Time) %>% mutate(compara_norm_inp = na_ma(compara_norm, k = 4, weighting = "linear", maxgap = 10)) %>% # fill in all time points and then linear imputation between (tried exponential and simple but get more odd bumps)
  group_by(strain, rep, exp) %>% dplyr::select(Time, strain,rep, exp, compara_norm_inp) %>% # Take imputed values
  pivot_wider(id_cols = c(strain, Time, rep), names_from = exp, values_from = compara_norm_inp) %>% mutate(nongrowth_only = CS - OD) %>% # look for difference between OD and heat output
  pivot_longer(cols = c("CS","OD"), names_to = "exp", values_to ="imput_val")

data_od_normd_ana <- left_join(data_od, data_od_normd)

ggplot(data_od_normd_ana, aes(x=Time, y = nongrowth_only, group = interaction(exp,rep, strain))) + 
  geom_line(aes(col = factor(rep)), lwd = 1) + 
  geom_line(aes(y = imput_val,col = interaction(rep))) + 
  facet_wrap(~strain, scales = "free", ncol = 2) + 
  scale_x_continuous("Time (h)") 
ggsave("plots/ODvsCS_inoc5_data_nongrowth_togplot.pdf")


g1 <- ggplot(data_od_normd_ana, aes(x=Time, y = nongrowth_only, group = interaction(rep, strain))) + 
  geom_line(aes(col = interaction(rep))) + 
  facet_wrap(~strain, ncol = 4) + 
  scale_x_continuous("Time (h)") + 
  scale_y_continuous("Non growth only") + 
  scale_color_discrete("Replicate") + 
  geom_hline(yintercept = c(0,0.5)) + 
  geom_vline(xintercept = c(25000,30000))

g2 <- ggplot(data_od_normd_ana, aes(x=Time, group = interaction(exp,rep, strain))) + 
  geom_line(aes(y = imput_val,col = interaction(rep), linetype = exp)) + 
  facet_wrap(~strain, ncol = 4) + 
  scale_x_continuous("Time (h)") + 
  scale_y_continuous("Normalised measure") + 
  scale_color_discrete("Experiment and\nreplicate")

g3a / g2 / g1 
ggsave("plots/ODvsCS_inoc5_data_nongrowth_tog_grid.pdf", width = 20, height = 20)

###### Norm by subtracting the "normal curve" = 11257
# # Normalise
# data_od_11257 <- data_od %>% filter(strain == 11257) %>% ungroup()
# data_11257 <- data_od_11257 %>% select(Time, rep, exp, compara) %>% rename(comp_11257 = compara)
# 
# data_od_sub <- data_od %>% group_by(strain, rep, exp) %>% left_join(data_11257) %>% 
#   mutate(compara_sub = compara - comp_11257,# subtract "normal" curve
#          max_v_s = max(compara_sub, na.rm = TRUE), # then normalise
#         compara_sub_norm = compara / max_v_s) 
# 
# ggplot(data_od_sub, aes(x=Time, y = compara_sub, group = interaction(rep, exp, strain))) + 
#   geom_line(aes(col = factor(rep), lty = exp), lwd = 1) + 
#   facet_wrap(~strain, scales = "free", ncol = 2) + 
#   scale_x_continuous("Time (h)") + 
#   scale_color_discrete("Replicate") + 
#   scale_linetype("Data")
# ggsave("plots/ODvsCS_inoc5_data_compare_sub.pdf")
# 
# ggplot(data_od_sub, aes(x=Time, y = compara_sub_norm, group = interaction(rep, exp, strain))) + 
#   geom_line(aes(col = factor(rep), lty = exp), lwd = 1) + 
#   facet_wrap(~strain, scales = "free", ncol = 2) + 
#   scale_x_continuous("Time (h)") + 
#   scale_color_discrete("Replicate") + 
#   scale_linetype("Data")
# ggsave("plots/ODvsCS_inoc5_data_compare_subn.pdf")

# Subtract normalised data? Need to complete: measured at different time points

#data_play <- data_od[c(1:100,1149:1259),] %>% filter(Time > cutoff_time)
#data_play <- data_play %>% ungroup() %>% complete(exp, Time) %>% mutate(compara_norm_inp = na_ma(compara_norm, k = 2, weighting = "exponential", maxgap = Inf)) %>% print(n=Inf)
#ggplot(data_play, aes(x=Time, y = compara_norm, group = exp)) + geom_line(aes(group = exp, col = exp)) + 
#  geom_line(aes(y = compara_norm_inp, group = exp, col = exp), linetype = 2) 

# data_od_sub_normd <- data_od_sub %>% ungroup() %>% 
#   complete(rep, strain, exp, Time) %>% mutate(compara_norm_inp = na_ma(compara_sub_norm, k = 4, weighting = "linear", maxgap = 10)) %>% # fill in all time points and then linear imputation between (tried exponential and simple but get more odd bumps)
#   group_by(strain, rep, exp) %>% dplyr::select(Time, strain,rep, exp, compara_norm_inp) %>% # Take imputed values
#   pivot_wider(id_cols = c(strain, Time, rep), names_from = exp, values_from = compara_norm_inp) %>% mutate(nongrowth_only = CS - OD) %>% # look for difference between OD and heat output
#   pivot_longer(cols = c("CS","OD"), names_to = "exp", values_to ="imput_val")
# 
# data_od_sub_normd_ana <- left_join(data_od, data_od_sub_normd)
# 
# ggplot(data_od_sub_normd_ana, aes(x=Time, y = nongrowth_only, group = interaction(exp,rep, strain))) + 
#   geom_line(aes(col = factor(rep)), lwd = 1) + 
#   geom_line(aes(y = imput_val,col = interaction(rep))) + 
#   facet_wrap(~strain, scales = "free", ncol = 2) + 
#   scale_x_continuous("Time (h)") 
# ggsave("plots/ODvsCS_inoc5_data_nongrowth_togplot_sub.pdf")
# 
# 
# g1 <- ggplot(data_od_sub_normd_ana, aes(x=Time, y = nongrowth_only, group = interaction(rep, strain))) + 
#   geom_line(aes(col = interaction(rep))) + 
#   facet_wrap(~strain, ncol = 4) + 
#   scale_x_continuous("Time (h)") + 
#   scale_y_continuous("Non growth only") + 
#   scale_color_discrete("Replicate") + 
#   geom_hline(yintercept = 0)
# 
# g2 <- ggplot(data_od_sub_normd_ana, aes(x=Time, group = interaction(exp,rep, strain))) + 
#   geom_line(aes(y = imput_val,col = interaction(rep), linetype = exp)) + 
#   facet_wrap(~strain, scales = "free", ncol = 4) + 
#   scale_x_continuous("Time (h)") + 
#   scale_y_continuous("Normalised measure") + 
#   scale_color_discrete("Experiment and\nreplicate")
# 
# g2 / g1 
# ggsave("plots/ODvsCS_inoc5_data_nongrowth_tog_grid_sub.pdf", width = 20)

### Extract characteristics - not sure needed now / priority. 

## What are the strains?
u <- as.character(unique(data_od$strain,stringsasFactors = FALSE)) # strains
## How many replicates? 
r <- unique(data_od$rep) # replicates
## How many experimental conditions? 
ex <- unique(data_od$exp)

# Where the parameters for each strain are stored
param <- matrix(0, length(u)*length(r)*length(ex), 50); # number of strains x number of replicates x number of experimental conditions
index <- 1 # for counting 
max <- c() # for calibration

data_od <- data_od %>% ungroup()

## Run thru each strain/rep/drying time/ inoculum: fit a growth curve to each
for(jj in 1:length(u)){ # for each strain
  r <- unlist(unique(data_od %>% filter(strain == u[jj]) %>% dplyr::select(rep))[,1]) # which replicates for this strain (all have different names)
  if(length(r) < 3){print(paste0(u[jj], " has too few replicates"))}else{ # filter out if only 2 replicates
    for(ii in 1:length(r)){ # for each replicate
      for(kk in 1:length(ex)){ #each of the experimental conditions
        
        strain <- u[jj];
        replicate <- r[ii]
        exper <- ex[kk]
        
        wi <- intersect(which(data_od$strain == strain),which(data_od$rep == replicate)) # if fit to each replicate
        w <- intersect(wi, which(data_od$exp == exper))
        
        if(length(w) > 0){ # if this replicate exists for this strain (i.e. there is data)
          data1 <- data_od[w,] # Grab data
          if(length(which(is.na(data1$compara)))>0){data1 <- data1[-which(is.na(data1$compara)),]} # remove any NA
          
          print(c(jj, strain, replicate, exper)) # output so can track how it is working
          p <- cut_extract_dp(data1, "Time", "compara", paste(strain, replicate, condition, inocl,sep="_")) ### NEW function: runs on any timeseries: gives baseline and cut parameter analysis
          
          ## Required parameters
          
          param[index,] <- c(strain, replicate, exper, unlist(p$param))
          index <- index + 1 # counting for storing matrix - next row each iteration
        }
        
      }
    }
  }
}
# The original set 
param_orig <- param

## Fitted parameters
param <- as.data.frame(param)
colnames(param) <- c("strain","rep","experiment",
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
param$valpeak <- as.numeric(param$valpeak)
param$timepeak <- as.numeric(param$timepeak)
dim(param)

## Store so don't have to run above
write.csv(param, "output/param_od.csv")
#param <- read_csv("output/param_od.csv")[,-1]

### look at OD vs CS
param_wide_timepeak <- param %>% pivot_wider(id_cols = c(strain, rep),names_from = experiment, values_from = c(timepeak))
ggplot(param, aes(x=interaction(strain,rep), y = timepeak)) + geom_point(aes(col = experiment))

ggplot(param_wide_timepeak, aes(x=CS, y = OD)) + geom_point(aes(col = strain)) + 
  scale_y_continuous(lim = c(0,31000)) + 
  geom_smooth(method='lm', formula= y~x)

ggplot(param, aes(x=interaction(strain,rep), y = exp_gr)) + geom_point(aes(col = experiment))

### clustering? 
data_od$drytime <- 0
param$drytime <- 0
data_od$inoc <- 5
param$inoc <- 5
data_od$strain <- as.character(data_od$strain)
param$strain <- as.character(param$strain)
data_od$value_J <- data_od$compara


c_cs <- cluster(data_od %>% filter(exp == "CS"), param %>% filter(experiment == "CS"),plot_where = "plots/CS_")
c_od <- cluster(data_od %>% filter(exp == "OD"), param %>% filter(experiment == "OD"),plot_where = "plots/OD_")


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

### Look at links to inocula
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


### Extract characteristics

## What are the strains?
u <- as.character(unique(data_od$strain,stringsasFactors = FALSE)) # strains
## How many replicates? 
r <- unique(data_od$rep) # replicates
## How many experimental conditions? 
ex <- unique(data_od$exp)
## How many inoc?
inn <- unique(data_od$inoc)

# Where the parameters for each strain are stored
param <- matrix(0, length(u)*length(r)*length(ex)*length(inn), 51); # number of strains x number of replicates x number of experimental conditions
index <- 1 # for counting 
max <- c() # for calibration

data_od <- data_od %>% dplyr::ungroup()

## Remove odd one that doesn't peak.. 
w <- intersect(which(data_od$strain == 11051), which(data_od$inoc == 1))
data_od <- data_od[-w,]

## Run thru each strain/rep/drying time/ inoculum: fit a growth curve to each
for(jj in 1:length(u)){ # for each strain
  r <- unlist(unique(data_od %>% filter(strain == u[jj]) %>% dplyr::select(rep))[,1]) # which replicates for this strain (all have different names)
  if(length(r) < 3){print(paste0(u[jj], " has too few replicates"))}else{ # filter out if only 2 replicates
    for(ii in 1:length(r)){ # for each replicate
      for(pp in 1:length(inn)){ # for each inoc
        for(kk in 1:length(ex)){ #each of the experimental conditions
          
          data1 <- data_od %>% filter(strain == u[jj], rep == r[ii], inoc == inn[pp], exp == ex[kk])
          
          if(dim(data1)[1] > 0){ # if this replicate exists for this strain (i.e. there is data)
            
            if(length(which(is.na(data1$compara)))>0){data1 <- data1[-which(is.na(data1$compara)),]} # remove any NA
            
            print(c(jj, u[jj], r[ii], inn[pp],ex[kk])) # output so can track how it is working
            p <- cut_extract_dp(data1, "Time", "compara", paste(strain, replicate, condition, inocl,sep="_")) ### NEW function: runs on any timeseries: gives baseline and cut parameter analysis
            
            ## Required parameters
            
            param[index,] <- c(u[jj], r[ii], inn[pp],ex[kk], unlist(p$param))
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
colnames(param) <- c("strain","rep","inoc","experiment",
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
param$valpeak <- as.numeric(param$valpeak)
param$timepeak <- as.numeric(param$timepeak)
dim(param)

## Store so don't have to run above
write.csv(param, "output/param_od_allinoc.csv")
param <- read_csv("output/param_od.csv")[,-1]

### look at OD vs CS
param_wide_timepeak <- param %>% pivot_wider(id_cols = c(strain, rep),names_from = experiment, values_from = c(timepeak))
ggplot(param, aes(x=interaction(strain,rep), y = timepeak)) + geom_point(aes(col = experiment))

ggplot(param_wide_timepeak, aes(x=CS, y = OD)) + geom_point(aes(col = strain)) + 
  scale_y_continuous(lim = c(0,31000)) + 
  geom_smooth(method='lm', formula= y~x)

ggplot(param, aes(x=interaction(strain,rep), y = exp_gr)) + geom_point(aes(col = experiment))

### clustering? 
data_od$drytime <- 0
param$drytime <- 0
data_od$inoc <- 5
param$inoc <- 5
data_od$strain <- as.character(data_od$strain)
param$strain <- as.character(param$strain)
data_od$value_J <- data_od$compara


c_cs <- cluster(data_od %>% filter(exp == "CS"), param %>% filter(experiment == "CS"),plot_where = "plots/CS_")
c_od <- cluster(data_od %>% filter(exp == "OD"), param %>% filter(experiment == "OD"),plot_where = "plots/OD_")


