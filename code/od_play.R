### Exploring OD data

library(zoo)
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


### look at OD vs CS
param_wide_timepeak <- param %>% pivot_wider(id_cols = c(strain, rep),names_from = experiment, values_from = c(timepeak))
ggplot(param, aes(x=interaction(strain,rep), y = timepeak)) + geom_point(aes(col = experiment))

ggplot(param_wide_timepeak, aes(x=CS, y = OD)) + geom_point(aes(col = strain)) + 
  scale_y_continuous(lim = c(0,31000)) + 
  geom_smooth(method='lm', formula= y~x)

### clustering? 
data_od$drytime <- 0
param$drytime <- 0
data_od$inoc <- 5
param$inoc <- 5
data_od$strain <- as.character(data_od$strain)
data_od$value_J <- data_od$compara


c_cs <- cluster(data_od %>% filter(exp == "CS"), param %>% filter(experiment == "CS"),plot_where = "plots/CS_")
c_od <- cluster(data_od %>% filter(exp == "OD"), param %>% filter(experiment == "OD"),plot_where = "plots/OD_")


