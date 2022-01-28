###### Analysis of delay / hidden second peaks
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
theme_set(theme_bw(base_size=14)) # theme setting for plots: black and white (bw) and font size (24)
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(14)

setwd(here::here())

## data? 
ddm <- read.csv(paste0("output/","cut_","all_time_series_fit_params.csv"))[,-1]

ddm %>% filter(odd_type_db != "0") %>% summarise(unique(strain)) # ave some double peak
ddm %>% group_by(strain) %>% mutate(no_double = ifelse(all(odd_type_db == 0),1,0)) %>% ungroup() %>% 
  filter(no_double == 1) %>% summarise(unique(strain)) # have no double peak

ddm %>% filter(strain == "11068") %>% ggplot(aes(x=Time, y = value, group = rep)) + geom_line(aes(linetype = factor(odd_type_db),col = factor(rep))) + facet_wrap(inoc ~ drytime, ncol = 2)

ddm %>% filter(strain == "11273") %>% ggplot(aes(x=Time, y = value, group = rep)) + geom_line(aes(linetype = factor(odd_type_db),col = factor(rep))) + facet_wrap(inoc ~ drytime, ncol = 2)

### param = data? 
param <- read.csv("output/cut_all_model_fit_params.csv")[,-1]



### odd_double <- 1 when can fit two normal or some other combo to the curve - is this really double peak? 
# For those without a double peak... width increases with inoculum... 
ggplot(param, aes(x=inocl, y = width_peak)) + geom_boxplot(aes(group = inocl)) + 
  facet_wrap(~odd_double)
ggsave("data_paper2/plots/hidden_widthvsinoc.pdf")

ggplot(param, aes(x=inocl, y = width_peak)) + geom_point(aes(group = inocl)) + 
  geom_smooth(method = lm) +
  facet_wrap(~odd_double)
ggsave("data_paper2/plots/hidden_widthvsinoc_line.pdf")

ggplot(param, aes(x=inocl, y = v_m_h_flow)) + geom_boxplot(aes(group = inocl)) + 
  facet_wrap(~odd_double)
ggsave("data_paper2/plots/hidden_heightvsinoc.pdf")

ggplot(param, aes(x=inocl, y = v_m_h_flow)) + geom_point(aes(group = inocl)) + 
  geom_smooth(method = lm) +
  facet_wrap(~odd_double)
ggsave("data_paper2/plots/hidden_heightvsinoc_line.pdf")

param2 <- read.csv("output/param_all_odd_labelled.csv") # just macotra
ggplot(param2, aes(x=inocl, y = v_m_h_flow)) + geom_point(aes(group = inocl)) + 
  geom_smooth(method = lm) +
  facet_wrap(~odd_double)
ggsave("data_paper2/plots/hidden_heightvsinoc_line_macotra.pdf")

### odd_peak <- 1 when have two close peaks 
# For those without a double peak... width increases with inoculum... 
ggplot(param, aes(x=inocl, y = width_peak)) + geom_boxplot(aes(group = inocl)) + 
  facet_wrap(~odd_peaks)
ggsave("data_paper2/plots/hidden_op_widthvsinoc.pdf")

ggplot(param, aes(x=inocl, y = width_peak)) + geom_point(aes(group = inocl)) + 
  geom_smooth(method = lm) +
  facet_wrap(~odd_peaks)
ggsave("data_paper2/plots/hidden_op_widthvsinoc_line.pdf")

ggplot(param, aes(x=inocl, y = v_m_h_flow)) + geom_boxplot(aes(group = inocl)) + 
  facet_wrap(~odd_peaks)
ggsave("data_paper2/plots/hidden_op_heightvsinoc.pdf")

ggplot(param, aes(x=inocl, y = v_m_h_flow)) + geom_point(aes(group = inocl)) + 
  geom_smooth(method = lm) +
  facet_wrap(~odd_peaks)
ggsave("data_paper2/plots/hidden_op_heightvsinoc_line.pdf")



##### Height of MAX peak vs inoculum for all strains
ggplot(param, aes(x=inocl, y = v_m_h_flow, group = interaction(rep, drytime))) + geom_point(aes(col = factor(drytime))) + 
  geom_line(aes(col = factor(drytime))) + facet_wrap(~strain_name) + 
  ggtitle("Height of max peak vs inoculum") + 
  scale_color_discrete("Drytime")
ggsave("data_paper2/plots/09_height_peak_vs_inoculum.pdf")

ggplot(param %>% filter(!strain_name %in% c("Newman","RWW12","SA3297","SA2704","RWW146","SAC042W", "Mu50", "M116")), 
       aes(x=inocl, y = v_m_h_flow, group = interaction(rep, drytime))) + geom_point(aes(col = factor(drytime))) + 
  geom_line(aes(col = factor(drytime))) + facet_wrap(~strain_name) + 
  ggtitle("Height of max peak vs inoculum just MACOTRA strains") + 
  scale_color_discrete("Drytime")
ggsave("data_paper2/plots/09_height_peak_vs_inoculum_macotra.pdf")

##### Height of FIRST peak/shoulder vs inoculum for all strains
param <- param %>% mutate(first_is_max = ifelse(v_m_h_flow == valpeak,1,0))

param %>% filter(first_is_max == 0) %>% summarise(n())
param %>% filter(first_is_max == 1) %>% summarise(n())


ggplot(param, aes(x=inocl, y = valpeak, group = interaction(rep, drytime))) + geom_point(aes(col = factor(drytime))) + 
  geom_line(aes(col = factor(drytime))) + facet_wrap(~strain_name) + 
  ggtitle("Height of first peak vs inoculum") + 
  scale_color_discrete("Drytime")
ggsave("data_paper2/plots/09_height_first_peak_vs_inoculum.pdf")

ggplot(param %>% filter(!strain_name %in% c("Newman","RWW12","SA3297","SA2704","RWW146","SAC042W", "Mu50", "M116")), 
       aes(x=inocl, y = valpeak, group = interaction(rep, drytime))) + geom_point(aes(col = factor(drytime))) + 
  geom_line(aes(col = factor(drytime))) + facet_wrap(~strain_name) + 
  ggtitle("Height of first peak vs inoculum just MACOTRA strains") + 
  scale_color_discrete("Drytime")
ggsave("data_paper2/plots/09_height_first_peak_vs_inoculum_macotra.pdf")

ggplot(param %>% filter(!strain_name %in% c("Newman","RWW12","SA3297","SA2704","RWW146","SAC042W", "Mu50", "M116")), 
       aes(x=inocl, y = valpeak, group = interaction(inocl, drytime))) + geom_boxplot(aes(col = factor(drytime))) + 
  ggtitle("Height of first peak vs inoculum just MACOTRA strains") + 
  scale_color_discrete("Drytime")
ggsave("data_paper2/plots/09_height_first_peak_vs_inoculum_macotra_boxplot.pdf")


# ggplot(param%>% filter(!strain_name %in% c("Newman","RWW12","SA3297","SA2704","RWW146","SAC042W", "Mu50", "M116")), 
#        aes(x=inocl, y = valpeak, group = interaction(strain_name, rep, drytime))) + geom_point(aes(col = factor(drytime))) + 
#   geom_line(aes(col = factor(drytime))) + facet_wrap(~first_is_max) + 
#   ggtitle("Height of first peak vs inoculum split by if first is max") + 
#   scale_color_discrete("Drytime")

ggplot(param%>% filter(!strain_name %in% c("Newman","RWW12","SA3297","SA2704","RWW146","SAC042W", "Mu50", "M116")), 
       aes(x=inocl, y = valpeak, group = interaction(inocl, drytime))) + geom_boxplot(aes(col = factor(drytime))) + 
  facet_wrap(~first_is_max) + 
  ggtitle("Height of first peak vs inoculum split by if first is max") + 
  scale_color_discrete("Drytime")
ggsave("data_paper2/plots/09_height_first_peak_vs_inoculum_macotra_boxplot_byfirstmax.pdf")

##### ADD filter to say time has to be large too 
param <- param %>% mutate(gap = t_m_h_flow - timepeak)

ggplot(param, aes(x=inocl, y = gap, group = inocl)) + geom_boxplot()
h <- hist(param$gap, breaks = seq(0,20,0.1))
h$counts

param <- param %>% mutate(first_diff_to_max8 = ifelse(gap > 0.8,1,0),first_diff_to_max5 = ifelse(gap > 0.5,1,0))


ggplot(param %>% filter(!strain_name %in% c("Newman","RWW12","SA3297","SA2704","RWW146","SAC042W", "Mu50", "M116")), 
       aes(x=inocl, y = valpeak, group = interaction(inocl, drytime))) + geom_boxplot(aes(col = factor(drytime))) + 
  facet_wrap(~first_diff_to_max8) + 
  ggtitle("Height of first peak vs inoculum just MACOTRA strains") + 
  scale_color_discrete("Drytime")
ggsave("data_paper2/plots/09_height_first_peak_vs_inoculum_macotra_boxplot_0.8_difference.pdf", width = 10, height = 5)

ggplot(param %>% filter(!strain_name %in% c("Newman","RWW12","SA3297","SA2704","RWW146","SAC042W", "Mu50", "M116")), 
       aes(x=inocl, y = valpeak, group = interaction(inocl, drytime))) + geom_boxplot(aes(col = factor(drytime))) + 
  facet_wrap(~first_diff_to_max5) + 
  ggtitle("Height of first peak vs inoculum just MACOTRA strains") + 
  scale_color_discrete("Drytime")
ggsave("data_paper2/plots/09_height_first_peak_vs_inoculum_macotra_boxplot_0.5_difference.pdf", width = 10, height = 5)


# ggplot(param%>% filter(!strain_name %in% c("Newman","RWW12","SA3297","SA2704","RWW146","SAC042W", "Mu50", "M116")), 
#        aes(x=inocl, y = valpeak, group = interaction(strain_name, rep, drytime))) + geom_point(aes(col = factor(drytime))) + 
#   geom_line(aes(col = factor(drytime))) + facet_wrap(~first_is_max) + 
#   ggtitle("Height of first peak vs inoculum split by if first is max") + 
#   scale_color_discrete("Drytime")

ggplot(param%>% filter(!strain_name %in% c("Newman","RWW12","SA3297","SA2704","RWW146","SAC042W", "Mu50", "M116")), 
       aes(x=inocl, y = valpeak, group = interaction(inocl, drytime))) + geom_boxplot(aes(col = factor(drytime))) + 
  facet_wrap(~first_is_max) + 
  ggtitle("Height of first peak vs inoculum split by if first is max") + 
  scale_color_discrete("Drytime")
ggsave("data_paper2/plots/09_height_first_peak_vs_inoculum_macotra_boxplot_byfirstmax.pdf")


##### Look at peaks for Valerie's odd ones 
val_odd = c(11006, 11016, 11051, 11142, 11210, 11214, 11280, 11040, 11274, 11014, 11285, 11046, 11161, 11048, 11057, 11050, 11165, 11283)
param_odd <- param2 %>% mutate(val_odd = ifelse(strain_name %in% val_odd, "Val odd", "Not odd" ))

ggplot(param_odd, 
       aes(x=inocl, y = valpeak, group = interaction(inocl, drytime))) + geom_boxplot(aes(col = factor(drytime))) + 
  facet_wrap(~val_odd) + 
  ggtitle("Height of first peak vs inoculum split by valerie's by eye odd classification") + 
  scale_color_discrete("Drytime")
ggsave("data_paper2/plots/09_height_first_peak_vs_inoculum_macotra_boxplot_byvalodd.pdf")

ggplot(param_odd %>% filter(!strain_name %in% c("Newman","RWW12","SA3297","SA2704","RWW146","SAC042W", "Mu50", "M116")), 
       aes(x=inocl, y = valpeak, group = interaction(inocl, val_odd, drytime))) + geom_boxplot(aes(col = factor(val_odd))) + 
  facet_wrap(~drytime) + 
  ggtitle("Height of first peak vs inoculum split by valerie's by eye odd classification (facet now drytime)") + 
  scale_color_manual("Odd", values = c("Val odd" = "darkgreen", "Not odd" = "orange"))
ggsave("data_paper2/plots/09_height_first_peak_vs_inoculum_macotra_boxplot_byvalodd_flip.pdf")

###### Check how Valerie and my odd peaks match up
# param2 = param with just macotra and some extra columns
param2 %>% filter(any_odd > 0) %>% summarise(unique(strain_name)) # Lots! 
#param_cnott <- 
param_odd <- param2 %>% 
  dplyr::select(strain_name, rep, drytime, inocl, any_odd, odd_type_db, remove_dataset_exp_iter, total_rep_rem) %>% group_by(strain_name, drytime,rep) %>% 
  mutate(n_datasets = n(), non_zero = sum(any_odd > 0), perc_odd = 100 * non_zero / n_datasets) %>% 
  #mutate(n_datasets = n(), non_zero = sum(odd_type_db > 0), perc_odd = 100 * non_zero / n_datasets) %>% # get way too many
  mutate(odd_rep_dt = ifelse(sum(perc_odd > 66)>1,1,ifelse(total_rep_rem>1,1,0))) %>% # if more than 2/3 odd or rep removed due to exponential growth issues then
  ungroup() %>%
  group_by(strain_name) %>%
  mutate(odd_strain = ifelse(sum(odd_rep_dt)>5,1,0)) %>% # this is three odd rep/drytimes
  mutate(number_odd_data_sets = sum(any_odd | remove_dataset_exp_iter > 0))
  
odd_strains <- param_odd %>% 
    filter(odd_strain == 1) %>% 
  ungroup() %>% 
  summarise(unique(strain_name))

odd_strains <- as.numeric(unlist(odd_strains[,1]))
length(odd_strains)
length(val_odd)
sort(val_odd)
sort(odd_strains)

setdiff(odd_strains, val_odd)
setdiff(val_odd, odd_strains)



###******* MODEL FITTING *************#######################################################################################################################################
## Fit growth curves to each set of data to determine the underlying parameters 
source("code/functions_for_heat_curves_additional_double_peak.R")
#### ADD INTO TO ORIGINAL
## (1) Peak 1 / peak 2 etc: time and value of peak
## (2) Delay between peaks

## What are the strains?
u <- as.character(unique(ddm$strain,stringsasFactors = FALSE)) # strains
## How many replicates? 
r <- unique(ddm$rep) # replicates
# What are the inoculums? 
q <- unique(ddm$inoc)

# Where the parameters for each strain are stored
param_dp <- matrix(0, length(u)*length(r)*length(q)*3, 48); # number of strains x number of replicates x number of experimental conditions
index <- 1 # for counting 
max <- c() # for calibration
drying_times <- c(0,168)

## Run thru each strain/rep/drying time/ inoculum: fit a growth curve to each
for(jj in 1:length(u)){ # for each strain
  r <- unique(ddm %>% filter(strain == u[jj]) %>% dplyr::select(rep))[,1] # which replicates for this strain (all have different names)
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
          
          print(c(strain, replicate, condition, inocl)) # output so can track how it is working
          p <- cut_extract_dp(data1, "Time", "value_J", paste(strain, replicate, condition, inocl,sep="_")) ### NEW function: runs on any timeseries: gives baseline and cut parameter analysis
          
          ## Required parameters
          
          param_dp[index,] <- c(strain, replicate, condition, inocl, p$param)
          index <- index + 1 # counting for storing matrix - next row each iteration
        }
        
      }
    }
  }
}

# The original set 
param_orig_dp <- param_dp

## Fitted parameters
param_dp <- as.data.frame(param_dp)
colnames(param_dp) <- c("strain_name","rep","drytime","inocl",
                     "t_m_h_flow", "v_m_h_flow", "exp_gr","lag","auc",
                     "odd_peaks","odd_width","width_peak","odd_shoulder","odd_double","shoulder_point_t","shoulder_point_v",
                     "cut_exp", "timepeak", "valpeak",
                     "mp_t1","mp_t2","mp_t3","mp_t4","mp_t5","mp_t6","mp_t7","mp_t8","mp_t9","mp_t10",
                     "mp_h1","mp_h2","mp_h3","mp_h4","mp_h5","mp_h6","mp_h7","mp_h8","mp_h9","mp_h10",
                     "gap1","gap2","gap3","gap4","gap5","gap6","gap7","gap8","gap9"
                     )

w<-which(param_dp$lag!=0); param_dp <- param_dp[w,] # remove 0 values
dim(param_dp)

write.csv(param_dp, "data_paper2/param_dp.csv")


#### Look at minor peaks
param_dp_mp <- param_dp %>% filter(!strain_name %in% c("Newman","RWW12","SA3297","SA2704","RWW146","SAC042W", "Mu50", "M116")) %>% 
  dplyr::select(-c(mp_t10, mp_h10)) %>% # makes the substr below complex
  dplyr::select(c("strain_name", "rep", "drytime", "inocl","cut_exp","timepeak","valpeak","mp_t1":"gap9")) %>% pivot_longer(cols = mp_t1:gap9) %>%
  mutate(peak = substr(name, nchar(name)-0, nchar(name))) %>%
  mutate(type = substr(name, nchar(name)-1, nchar(name)-1)) %>%
  filter(value > 0) %>%
  dplyr::select(-c(name)) %>% 
  pivot_wider(names_from = type, values_from = value) %>% 
  group_by(strain_name, rep, drytime, inocl) %>%
  mutate(n_peaks = max(peak)) %>%
  ungroup()

max(param_dp_mp$peak)
param_dp_mp %>% filter(peak>2)
param_dp_mp$n_peaks <- as.numeric(param_dp_mp$n_peaks)

### number of peaks at each inoculum? 
ggplot(param_dp_mp, aes(x=n_peaks, fill = inocl)) + geom_histogram(position = "identity", bins = 4) + facet_wrap(drytime~inocl, ncol = 3)
ggsave("data_paper2/plots/09_npeaks_by_inocl.pdf")

ggplot(param_dp_mp %>% filter(strain_name %in% val_odd), aes(x=n_peaks, fill = inocl)) + geom_histogram(position = "identity", bins = 4) + 
  facet_wrap(drytime~inocl, ncol = 3)
ggsave("data_paper2/plots/09_npeaks_by_inocl_val_odd.pdf")

ggplot(param_dp_mp %>% filter(!strain_name %in% val_odd), aes(x=n_peaks, fill = inocl)) + geom_histogram(position = "identity", bins = 4) + 
  facet_wrap(drytime~inocl, ncol = 3)
ggsave("data_paper2/plots/09_npeaks_by_inocl_val_notodd.pdf")

## gaps
param_dp_mp$time_peak <- as.numeric(param_dp_mp$t)
param_dp_mp$height_peak <- as.numeric(param_dp_mp$h)
param_dp_mp$gap <- as.numeric(param_dp_mp$p)

ggplot(param_dp_mp, aes(x=inocl, y = gap, group = interaction(peak, drytime))) + geom_point(aes(col = drytime)) + 
  facet_wrap(~drytime)

ggplot(param_dp_mp, aes(x=inocl, y = gap, group = interaction(peak, drytime))) + geom_boxplot() + 
  facet_wrap(~drytime)

ggplot(param_dp_mp %>% filter(n_peaks == 2), aes(x=inocl, y = gap, group = interaction(drytime))) + geom_point() + 
  geom_smooth(method = lm) +
  facet_wrap(~drytime)

ggplot(param_dp_mp, aes(x=inocl, y = gap, group = interaction(n_peaks, drytime))) + geom_point() + 
  geom_smooth(method = lm) +
  facet_wrap(drytime~n_peaks, ncol = 4)
ggsave("data_paper2/plots/09_gaps.pdf")

ggplot(param_dp_mp %>% filter(strain_name %in% val_odd), aes(x=inocl, y = gap, group = interaction(n_peaks, drytime))) + geom_point() + 
  geom_smooth(method = lm) +
  facet_wrap(drytime~n_peaks, ncol = 4)
ggsave("data_paper2/plots/09_gaps_val_odd.pdf")

ggplot(param_dp_mp  %>% filter(!strain_name %in% val_odd), aes(x=inocl, y = gap, group = interaction(n_peaks, drytime))) + geom_point() + 
  geom_smooth(method = lm) +
  facet_wrap(drytime~n_peaks, ncol = 4)
ggsave("data_paper2/plots/09_gaps_val_not_odd.pdf")


ggplot(param_dp_mp, aes(x=inocl, y = gap, group = interaction(n_peaks, drytime))) + geom_point(aes(col = factor(n_peaks))) + 
  geom_smooth(method = lm) +
  facet_wrap(~drytime, nrow = 2)

ggplot(param_dp_mp, aes(x=inocl, y = gap, group = interaction(n_peaks, inocl, drytime))) + geom_boxplot() + 
  facet_wrap(drytime~n_peaks, ncol = 4)
ggsave("data_paper2/plots/09_gaps_boxplot.pdf")

### Height? 
ggplot(param_dp_mp, aes(x=inocl, y = height_peak, group = interaction(n_peaks, drytime))) + geom_point() + 
  geom_smooth(method = lm) +
  facet_wrap(drytime~n_peaks, ncol = 4)
ggsave("data_paper2/plots/09_heights.pdf")

ggplot(param_dp_mp %>% filter(strain_name %in% val_odd), aes(x=inocl, y = height_peak, group = interaction(n_peaks, drytime))) + geom_point() + 
  geom_smooth(method = lm) +
  facet_wrap(drytime~n_peaks, ncol = 4)
ggsave("data_paper2/plots/09_heights_val_odd.pdf")

ggplot(param_dp_mp  %>% filter(!strain_name %in% val_odd), aes(x=inocl, y = height_peak, group = interaction(n_peaks, drytime))) + geom_point() + 
  geom_smooth(method = lm) +
  facet_wrap(drytime~n_peaks, ncol = 4)
ggsave("data_paper2/plots/09_heights_val_not_odd.pdf")


ggplot(param_dp_mp, aes(x=inocl, y = height_peak, group = interaction(n_peaks, inocl, drytime))) + geom_boxplot() + 
  facet_wrap(drytime~n_peaks, ncol = 4)
ggsave("data_paper2/plots/09_height_boxplot.pdf")

## Timing of peaks
ggplot(param_dp_mp, aes(x=peak, y = height_peak, group = interaction(n_peaks, drytime))) + geom_point(aes(col = factor(inocl))) + 
  geom_smooth(method = lm) +
  facet_wrap(drytime~n_peaks, ncol = 4)
ggsave("data_paper2/plots/09_basic_height.pdf")

ggplot(param_dp_mp, aes(x=peak, y = time_peak, group = interaction(n_peaks, drytime))) + geom_point(aes(col = factor(inocl))) + 
  geom_smooth(method = lm) +
  facet_wrap(drytime~n_peaks, ncol = 4) + 
  scale_y_continuous(lim = c(0,30))
ggsave("data_paper2/plots/09_basic_time.pdf")

write.csv(param_dp_mp, "data_paper2/param_dp_mp.csv")


