#### Pathway 1

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
theme_set(theme_bw(base_size=12)) # theme setting for plots: black and white (bw) and font size (24)
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(14)

setwd(here::here())

# From OD to CS
ddm_cs_orig <- read.csv("data_paper2/ddm_CS.csv")[,-1]
ddm_od_orig <- read.csv("data_paper2/ddm_OD.csv")[,-1]

# Smooth and Differentiate OD
data_od <- ddm_od_orig %>% group_by(variable, strain, inoc) %>%
  mutate(ma_value = rollapply(value, 10, mean,fill = NA),
         differ = c(0,diff(ma_value))) %>%
  ungroup() %>%
  filter(Time > 1, Time < 22)

ggplot(data_od, aes(x=Time, y = ma_value, group = interaction(rep, inoc,strain))) + 
  geom_line(aes(col = factor(inoc))) + 
  facet_wrap(~strain) + 
  scale_y_continuous("Optical density (600 nm)") +
  scale_x_continuous("Time (h)") + 
  scale_color_discrete("Inoculum", breaks = c(1,2,3,4,5,6,7), labels = c(str2expression("10^1"),str2expression("10^2"),str2expression("10^3"),
                                                                         str2expression("10^4"),str2expression("10^5"),
                                                                         str2expression("10^6"),str2expression("10^7"))) 

ggplot(data_od, aes(x=Time, y = differ, group = interaction(rep, inoc))) + 
  geom_line(aes(col = factor(inoc))) + 
  facet_wrap(~strain) + 
  scale_y_continuous("Difference in optical density (600 nm) per time step")+ 
  scale_x_continuous("Time (h)") + 
  scale_color_discrete("Inoculum", breaks = c(2,3,4,5,6), labels = c(str2expression("10^2"),str2expression("10^3"),
                                                                     str2expression("10^4"),str2expression("10^5"),
                                                                     str2expression("10^6"))) 
ggsave("data_paper2/plots/path_1_od_differentiated.pdf")


###################### CS
ggplot(ddm_cs_orig, aes(x=Time,y = value_J,group = interaction(rep, inoc,strain))) + geom_line(aes(col = factor(inoc))) + facet_wrap(~strain) + 
  scale_y_continuous("Heat flow") +
  scale_x_continuous("Time (h)") + 
  scale_color_discrete("Inoculum", breaks = c(1,2,3,4,5,6,7), labels = c(str2expression("10^1"),str2expression("10^2"),str2expression("10^3"),
                                                                         str2expression("10^4"),str2expression("10^5"),
                                                                         str2expression("10^6"),str2expression("10^7"))) 

ggsave("data_paper2/plots/path_1_cs_orig.pdf")


## Parameters
param <- read.csv("data_paper2/param_orig.csv")[,-1]
param$strain <- param$strain_name
param$inoc <- param$inocl
data_cs <- left_join(ddm_cs_orig, param, by = c("strain","rep","inoc"))

ggplot(data_cs, aes(x=Time,y = value_J,group = interaction(rep, inoc,strain))) + geom_line(aes(col = factor(inoc))) + 
  facet_grid(inoc~strain) + 
  scale_y_continuous("Heat flow") +
  scale_x_continuous("Time (h)") + 
  scale_color_discrete("Inoculum", breaks = c(1,2,3,4,5,6,7), labels = c(str2expression("10^1"),str2expression("10^2"),str2expression("10^3"),
                                                                         str2expression("10^4"),str2expression("10^5"),
                                                                         str2expression("10^6"),str2expression("10^7"))) +
  geom_vline(aes(xintercept = timepeak, col = factor(inoc))) + 
  geom_hline(aes(yintercept = valpeak, col = factor(inoc))) 
ggsave("data_paper2/plots/path_1_cs_cross.pdf")  


###******* MODEL FITTING *************#######################################################################################################################################
## Fit to differentiated OD data

# strain names
data_od$rep <- "0";
data_od[,c("strain_label","inoc_name")]<-colsplit(data_od$variable, "(?<=\\p{L})(?=[\\d+$])", c("char", "digit"))

w<- which(data_od$strain_label == "A")
data_od[w,"rep"] = "A"
w<- which(data_od$strain_label == "B")
data_od[w,"rep"] = "B"
w<- which(data_od$strain_label == "D")
data_od[w,"rep"] = "A"
w<- which(data_od$strain_label == "E")
data_od[w,"rep"] = "B"

## What are the strains?
u <- as.character(unique(data_od$strain,stringsasFactors = FALSE)) # strains
## How many replicates? 
r <- unique(data_od$rep) # replicates
# What are the inoculums? 
q <- unique(data_od$inoc)

# Where the parameters for each strain are stored
param <- matrix(0, length(u)*length(q)*length(r), 18); # number of strains x number of replicates x number of experimental conditions
index <- 1 # for counting 
max <- c() # for calibration

## Run thru each strain/rep/drying time/ inoculum: fit a growth curve to each
for(jj in 1:length(u)){ # for each strain
  r <- unique(data_od %>% filter(strain == u[jj]) %>% dplyr::select(rep))[,1] # which replicates for this strain (all have different names)
  for(ii in 1:length(r)){ # for each replicate
    #for(kk in c(1,2)){ #each of the experimental conditions
    for(ll in 1:length(q)){ #each of the inocula
      
      strain <- u[jj];
      replicate <- as.character(r[ii])
      #   condition <- drying_times[kk]
      inocl <- q[ll]
      
      wi <- intersect(which(data_od$strain == strain),which(data_od$rep == replicate)) # if fit to each replicate
      #  wj <- intersect(wi, which(data_od$drytime == condition))
      w <- intersect(wi, which(data_od$inoc == as.numeric(inocl)))
      
      if(length(w) > 0){ # if this replicate exists for this strain (i.e. there is data)
        data1 <- data_od[w,] # Grab data
        wna<-which(is.na(data1$differ)) # remove NA at beginning
        data1 <- data1[-wna,]
        
        print(c(strain, replicate, inocl)) # output so can track how it is working
        
        p <- cut_extract(data1, "Time", "differ", paste(strain, replicate, inocl,sep="_"),plot = 0, plot_where = "data_paper2/plots/path_1_", early_cut = 0) ### NEW function: runs on any timeseries: gives baseline and cut parameter analysis
        
        # print(p$param) # was getting errors generating 
        
        ## Required parameters
        param[index,] <- c(strain, replicate, inocl, p$param)
        index <- index + 1 # counting for storing matrix - next row each iteration
      }
    }
  }
}

# The original set 
param_orig <- param

## Fitted parameters
param <- as.data.frame(param)
colnames(param) <- c("strain_name","rep","inocl",
                     "t_m_h_flow", "v_m_h_flow", "exp_gr","lag","auc",
                     "odd_peaks","odd_width","width_peak","odd_shoulder","odd_double","shoulder_point_t","shoulder_point_v",
                     "cut_exp", "timepeak", "valpeak")

w<-which(param$lag!=0); param <- param[w,] # remove 0 values
dim(param)
write.csv(param,"data_paper2/param_od.csv")

## Parameters
param$strain <- param$strain_name
param$inoc <- as.numeric(param$inocl)
data_od <- left_join(data_od, param, by = c("strain","rep","inoc"))

ggplot(data_od, aes(x=Time,y = differ,group = interaction(rep, inoc,strain))) + geom_line(aes(col = factor(inoc))) + 
  facet_grid(inoc~strain) + 
  scale_y_continuous("Difference in optical density (600 nm) per time step") +
  scale_x_continuous("Time (h)") + 
  scale_color_discrete("Inoculum", breaks = c(1,2,3,4,5,6,7), labels = c(str2expression("10^1"),str2expression("10^2"),str2expression("10^3"),
                                                                         str2expression("10^4"),str2expression("10^5"),
                                                                         str2expression("10^6"),str2expression("10^7"))) +
  geom_vline(aes(xintercept = as.numeric(timepeak), col = factor(inoc))) + 
  geom_hline(aes(yintercept = as.numeric(valpeak), col = factor(inoc)))
ggsave("data_paper2/plots/path_1_od_cross.pdf")

######## To scale up differentiated OD compare valpeaks 
param_od <- read.csv("data_paper2/param_od.csv")[,-1]
param_cs <- read.csv("data_paper2/param_orig.csv")[,-1]

param_od$od_valpeak <- param_od$valpeak
param_cs$cs_valpeak <- param_cs$valpeak

param_both <- left_join(param_od[,c("strain_name","rep","inocl","od_valpeak")], param_cs[,c("strain_name","rep","inocl","cs_valpeak")]) %>%
  mutate(ratio = od_valpeak/cs_valpeak)

param_both$inoc <- param_both$inocl

########## Add to data
data_od_b <- data_od
data_od_b$meas <- "OD"
data_od_b$value <- data_od_b$differ

data_cs_b <- data_cs %>% dplyr::select(-c("value"))
data_cs_b$value <- data_cs_b$value_J
data_cs_b <- data_cs_b %>% dplyr::select(-c("value_J","csum"))
data_cs_b$meas <- "CS"

keep <- c("strain_name","Time","value","variable","inoc","rep","timepeak","valpeak","meas","lag")
data_both <- rbind(data_od_b[,keep],data_cs_b[,keep])
data_both <- left_join(data_both, param_both)
data_both <- data_both %>% mutate(conv_value = ifelse(meas == "CS",value * ratio,value)) 

ggplot(data_both, aes(x=Time, y = conv_value, group = interaction(meas, rep, inoc,strain_name))) + 
  geom_line(aes(col = factor(inoc),linetype = factor(meas))) + 
  facet_wrap(~strain_name) + 
  scale_y_continuous("CS + converted OD") +
  scale_x_continuous("Time (h)") + 
  scale_color_discrete("Inoculum", breaks = c(1,2,3,4,5,6,7), labels = c(str2expression("10^1"),str2expression("10^2"),str2expression("10^3"),
                                                                         str2expression("10^4"),str2expression("10^5"),
                                                                         str2expression("10^6"),str2expression("10^7")))
ggsave("data_paper2/plots/path_1_od_cs_peak_val_scaled.pdf")

ggplot(data_both, aes(x=Time, y = conv_value, group = interaction(meas, rep, inoc,strain_name))) + 
  geom_line(aes(col = factor(inoc),linetype = factor(meas))) + 
  facet_grid(inoc~strain_name) + 
  scale_y_continuous("CS + converted OD") +
  scale_x_continuous("Time (h)") + 
  scale_color_discrete("Inoculum", breaks = c(1,2,3,4,5,6,7), labels = c(str2expression("10^1"),str2expression("10^2"),str2expression("10^3"),
                                                                         str2expression("10^4"),str2expression("10^5"),
                                                                         str2expression("10^6"),str2expression("10^7")))
ggsave("data_paper2/plots/path_1_od_cs_peak_val_scaled_panels.pdf")

#plot(data_od$Time, data_od$differ)

## subtract lag time
data_both$lag <- as.numeric(data_both$lag)
data_both$adj_time <- data_both$Time - data_both$lag
ggplot(data_both, aes(x=adj_time, y = conv_value, group = interaction(meas, rep, inoc,strain_name))) + 
  geom_line(aes(col = factor(inoc),linetype = factor(meas))) + 
  facet_grid(inoc~strain_name) + 
  scale_y_continuous("CS + converted OD") +
  scale_x_continuous("Time (h)") + 
  scale_color_discrete("Inoculum", breaks = c(1,2,3,4,5,6,7), labels = c(str2expression("10^1"),str2expression("10^2"),str2expression("10^3"),
                                                                         str2expression("10^4"),str2expression("10^5"),
                                                                         str2expression("10^6"),str2expression("10^7")))
ggsave("data_paper2/plots/path_1_od_cs_peak_val_scaled_panels_adj_lag.pdf")

## Subtract OD from CS
# issue is that time points not the same
data_both_sub <- data_both %>% filter(!is.na(conv_value)) %>% 
  dplyr::select(c(strain_name, inoc, rep, meas, adj_time, conv_value)) %>%
  group_by(strain_name, inoc, rep) 

meass <- c("OD","CS")
strains <- unique(data_both_sub$strain_name)
inocs <- unique(data_both_sub$inoc)
reps <- unique(data_both_sub$rep)

interp_data <- c()

minn <- min(data_both_sub$adj_time)
maxx <- max(data_both_sub$adj_time)

for(i in 1:length(strains)){
  for(j in 1:length(meass)){
    for(k in 1:length(inocs)){
      for(l in 1:length(reps)){
        
        this_data <- data_both_sub %>% filter(strain_name == strains[i], meas == meass[j], inoc == inocs[k], rep == reps[l])
        if(dim(this_data)[1]>0){
          new_data <- approx(this_data$adj_time, this_data$conv_value, xout = seq(round(minn,0),round(maxx,0),0.1))
          w <- which(!is.na(new_data$y))
          interp_data <- rbind(interp_data, 
                               cbind(strains[i], meass[j], inocs[k], reps[l], new_data$x[w],new_data$y[w]))
        }
        
      }
    }
  }
}

interp_data <- as.data.frame(interp_data)
colnames(interp_data) <- c("strain_name","meas","inoc","rep","intrp_time","intrp_value")
interp_data$intrp_value <- as.numeric(interp_data$intrp_value)
interp_data$intrp_time <- as.numeric(interp_data$intrp_time)

## using complete etc didn't work
# data_both_sub <- data_both %>% filter(!is.na(conv_value)) %>% 
#   dplyr::select(c(strain_name, adj_time, inoc, rep, meas, conv_value)) %>%
#   group_by(strain_name, inoc, rep) %>% 
#   complete(adj_time, nesting(meas)) %>% 
#   group_by(strain_name, inoc, rep, meas) %>% 
#   mutate(conv_value = na.approx(conv_value))
# pivot_wider(id_cols = c(strain_name, adj_time, variable, inoc, rep), names_from = meas, values_from = conv_value)

ggplot(interp_data, aes(x=intrp_time, y = intrp_value, group = interaction(meas, rep, inoc,strain_name))) + 
  geom_point(aes(col = factor(inoc),linetype = factor(meas))) + 
  facet_grid(inoc~strain_name) + 
  scale_y_continuous("CS + converted OD") +
  scale_x_continuous("Time (h)") + 
  scale_color_discrete("Inoculum", breaks = c(1,2,3,4,5,6,7), labels = c(str2expression("10^1"),str2expression("10^2"),str2expression("10^3"),
                                                                         str2expression("10^4"),str2expression("10^5"),
                                                                         str2expression("10^6"),str2expression("10^7")))

### subtract!
subtracted_data <- interp_data %>% pivot_wider(id_cols = c(strain_name, inoc, rep, intrp_time), names_from = meas, values_from = intrp_value) %>%
  mutate(growth_heat = OD, non_growth_heat = CS - OD) %>% 
  pivot_longer(cols = c(OD,CS,growth_heat,non_growth_heat))

ggplot(subtracted_data %>% filter(name %in% c("growth_heat","non_growth_heat")), aes(x=intrp_time, y = value)) + 
  geom_line(aes(col = factor(name))) + 
  facet_grid(inoc~strain_name) + 
  scale_y_continuous("CS + converted OD") +
  scale_x_continuous("Time (h)") + 
  scale_color_discrete("Type of heat output") + 
  geom_line(data = subtracted_data %>% filter(name == "CS"), col = "grey")

ggsave("data_paper2/plots/path_1_subtracted.pdf")

subtracted_data <- subtracted_data %>% mutate(low_inoc = ifelse(inoc < 3,1,0))
ggplot(subtracted_data %>% filter(name %in% c("growth_heat","non_growth_heat")), aes(x=intrp_time, y = value, group = interaction(strain_name, name, inoc))) + 
  geom_line(aes(col = factor(inoc))) + 
  facet_grid(name~low_inoc + strain_name, scales = "free") + 
  scale_y_continuous("CS + converted OD") +
  scale_x_continuous("Time (h)") + 
  scale_color_discrete("Type of heat output") + 
  geom_line(data = subtracted_data %>% filter(name == "CS"), col = "grey")

