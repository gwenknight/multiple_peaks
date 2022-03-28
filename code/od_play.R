### Exploring OD data

library(zoo)
library(imputeTS) # to give weighter moving average - exponential weighting of those further away
library(patchwork)
setwd(here::here())

# Data to explore
data_od_orig <- read_csv("data/growth_ODvsCS_20220224.csv")[,-1]

## Look at OD & CS data
ggplot(data_od_orig, aes(x=Time, y = value, group = interaction(inoc, exp, strain, rep))) + geom_line(aes(col = factor(inoc))) + 
  facet_grid(exp~strain, scales = "free")
ggsave("plots/ODvsC2_rawdata.pdf")

ggplot(data_od_orig %>% filter(inoc == 5), aes(x=Time, y = value, group = interaction(inoc, exp, strain, rep))) + geom_line(aes(col = factor(rep))) + 
  facet_grid(exp~strain, scales = "free")
ggsave("plots/ODvsC2_inoc5_rawdata.pdf")

data_od <- data_od_orig %>% dplyr::select(Time, rep, exp, value, strain, inoc) %>% 
  filter(inoc == 5) %>% # only 10^5 for this analysis
  group_by(rep, exp, strain) %>%
  mutate(ma_value = rollapply(value, 10, mean,fill = NA),
         differ = c(0,diff(ma_value)),
         compara = ifelse(exp == "CS", value, differ)) %>%
  ungroup() 



g2 <- ggplot(data_od, aes(x=Time, y = ma_value, group = interaction(rep, exp, strain))) + 
  geom_line(aes(col = factor(rep))) + 
  facet_grid(exp~strain,scales = "free") + 
  scale_x_continuous("Time (h)") 
ggsave("plots/ODvsCS_inoc5_data_smoothed.pdf")

g3 <- ggplot(data_od, aes(x=Time, y = compara, group = interaction(rep, exp, strain))) + 
  geom_line(aes(col = factor(rep))) + 
  facet_grid(exp~strain, scales = "free") + 
  scale_x_continuous("Time (h)") 
ggsave("plots/ODvsCS_inoc5_data_compare.pdf")

# Try to plot together? scale by max. Normalise
data_od <- data_od %>% group_by(strain, rep, exp) %>% mutate(max_v = max(compara, na.rm = TRUE),
                                                             compara_norm = compara / max_v)

ggplot(data_od, aes(x=Time, y = compara_norm, group = interaction(rep, exp, strain))) + 
  geom_line(aes(col = interaction(exp, rep), lty = exp), lwd = 1) + 
  facet_wrap(~strain, scales = "free", ncol = 2) + 
  scale_x_continuous("Time (h)") 
ggsave("plots/ODvsCS_inoc5_data_compare_norm.pdf")

# Subtract normalised data? Need to complete: measured at different time points
# remove all data before the time point that they all have which is the max of the minimum times to avoid odd completion curves prior to start and after
cutoff_time_dn = max(data_od %>% summarise(minn = min(Time)) %>% ungroup() %>% dplyr::select(minn))
cutoff_time_up = min(data_od %>% summarise(maxx = max(Time)) %>% ungroup() %>% dplyr::select(maxx))

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
  geom_line(aes(col = interaction(exp,rep)), lwd = 1) + 
  geom_line(aes(y = imput_val,col = interaction(exp,rep))) + 
  facet_wrap(~strain, scales = "free", ncol = 2) + 
  scale_x_continuous("Time (h)") 
ggsave("plots/ODvsCS_inoc5_data_nongrowth_togplot.pdf")


g1 <- ggplot(data_od_normd_ana, aes(x=Time, y = nongrowth_only, group = interaction(rep, strain))) + 
  geom_line(aes(col = interaction(rep))) + 
  facet_wrap(~strain, ncol = 4) + 
  scale_x_continuous("Time (h)") + 
  scale_y_continuous("Non growth only") + 
  scale_color_discrete("Replicate")

g2 <- ggplot(data_od_normd_ana, aes(x=Time, group = interaction(exp,rep, strain))) + 
  geom_line(aes(y = imput_val,col = interaction(rep), linetype = exp)) + 
  facet_wrap(~strain, scales = "free", ncol = 4) + 
  scale_x_continuous("Time (h)") + 
  scale_y_continuous("Normalised measure") + 
  scale_color_discrete("Experiment and\nreplicate")

g2 / g1 
ggsave("plots/ODvsCS_inoc5_data_nongrowth_tog_grid.pdf")


### Extract characteristics

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


