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
theme_set(theme_bw(base_size=8)) # theme setting for plots: black and white (bw) and font size (24)
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(14)

setwd(here::here())

# 11016 data
# OD
# 3 reps  7 inoc / 4 gluc
# cleanded in double_peaks_111016_dex_cleaning.R
ds <- read.csv("data_paper2/11016_glucose.csv")[,-1]

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
# ggplot(data_ds, aes(x=Time, y = ma_value, group = interaction(rep, inoc,glucose))) +
#   geom_line(aes(col = factor(glucose))) +
#   facet_wrap(~rep) +
#   scale_y_continuous("Heat flow") +
#   scale_x_continuous("Time (h)") #+
# #scale_color_discrete("Inoculum", breaks = c(1,2,3,4,5,6,7), labels = c(str2expression("10^1"),str2expression("10^2"),str2expression("10^3"),
# #                                     str2expression("10^4"),str2expression("10^5"),
# #                                     str2expression("10^6"),str2expression("10^7"))) 
# 
# ggplot(data_ds, aes(x=Time, y = differ, group = interaction(rep, inoc, glucose))) +
#   geom_line(aes(col = factor(inoc))) +
#   facet_wrap(glucose~rep, ncol = 3) +
#   scale_y_continuous("Difference in heatflow per time step")+
#   scale_x_continuous("Time (h)") +
#   scale_color_discrete("Inoculum", breaks = c(1,2,3,4,5,6,7), labels = c(str2expression("10^1"),str2expression("10^2"),str2expression("10^3"),
#                                                                          str2expression("10^4"),str2expression("10^5"),
#                                                                          str2expression("10^6"),str2expression("10^7")))
# ggsave("data_paper2/plots/11016sugar_od_differentiated_byrep.pdf")
## Smoother differ
# ggplot(data_ds, aes(x=Time, y = ma_differ, group = interaction(rep, inoc, glucose))) +
#   geom_line(aes(col = factor(inoc))) +
#   facet_wrap(glucose~rep, ncol = 3) +
#   scale_y_continuous("Difference in heatflow per time step")+
#   scale_x_continuous("Time (h)") +
#   scale_color_discrete("Inoculum", breaks = c(1,2,3,4,5,6,7), labels = c(str2expression("10^1"),str2expression("10^2"),str2expression("10^3"),
#                                                                          str2expression("10^4"),str2expression("10^5"),
#                                                                          str2expression("10^6"),str2expression("10^7")))
# ggsave("data_paper2/plots/11016sugar_od_differentiatedsmoothed_byrep.pdf")
# 
# ggplot(data_ds, aes(x=Time, y = differ, group = interaction(rep, inoc, glucose))) + 
#   geom_line(aes(col = factor(inoc))) + 
#   facet_wrap(inoc~glucose, ncol = 4) + 
#   scale_y_continuous("Difference in heatflow per time step")+ 
#   scale_x_continuous("Time (h)") + 
#   scale_color_discrete("Inoculum", breaks = c(1,2,3,4,5,6,7), labels = c(str2expression("10^1"),str2expression("10^2"),str2expression("10^3"),
#                                                                          str2expression("10^4"),str2expression("10^5"),
#                                                                          str2expression("10^6"),str2expression("10^7"))) 
# ggsave("data_paper2/plots/11016sugar_od_differentiated_byinoc.pdf")


###### Find peaks
### Find peaks function: from https://github.com/stas-g/findPeaks
# if bigger than "m" number of points either side of it
find_peaks <- function (x, m = 3){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    # It must be bigger than or equal to all points m to the left and to the right
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}

## What are the glucose concentrations?
u <- unique(data_ds$glucose)
## How many replicates? 
r <- unique(data_ds$rep) # replicates
# What are the inoculums? 
q <- unique(data_ds$inoc)

# Where the parameters for each strain are stored
param <- matrix(0, length(u)*length(q)*length(r), 18); 
index <- 1 # for counting 

# Where store? 
data_ds$peaks <- 0
data_ds$max_peak <- 0
data_ds$scnd_peak <- 0


## Run thru each 
for(jj in 1:length(u)){ # for each glucose
  for(ii in 1:length(r)){ # for each replicate
    for(ll in 1:length(q)){ #each of the inocula
      
      gluc <- u[jj];
      replicate <- as.character(r[ii])
      inocl <- q[ll]
      
      wi <- intersect(which(data_ds$glucose == gluc),
                      which(data_ds$rep == replicate)) 
      w <- intersect(wi, which(data_ds$inoc == as.numeric(inocl)))
      
      if(length(w) > 0){ # if this replicate exists for this strain (i.e. there is data)
        data1 <- data_ds[w,] # Grab data
        wna<-which(is.na(data1$ma_differ)) # remove NA at beginning
        data1 <- data1[-wna,]
        wnew_ds <- w[-wna]
        
        print(c(gluc, replicate, inocl)) # output so can track how it is working
        
        peaks <- find_peaks(data1$ma_differ, m = 3)
        
        # Are the peaks really peaks? 
        # must be above a low level: 
        heights = data1[peaks,"ma_differ"]
        # correct for non-zero beginning
        max_height = max(heights)
        w_max = which(heights == max_height)
        # Pick out second even if low
        if(length(peaks) > 1){
          scnd_height = max(unlist(heights)[-w_max])
          w_scnd <- which(heights == scnd_height)
          data_ds[wnew_ds[peaks[w_scnd]],"scnd_peak"] <- 1
        }
        thresh_height = 0.07 * max_height + as.numeric(data1[1,"ma_differ"])
        w_min <- which(abs(heights) < thresh_height)
        if(length(w_min)>0){
          peaks <- peaks[-w_min]
        }
        # new heights with just the subset of "high" peaks
        heights = data1[peaks,"ma_differ"]
        w_max <- which(heights == max(heights))
        
        if(length(peaks)>1){
          data_ds[wnew_ds[peaks],"peaks"] <- 1
          data_ds[wnew_ds[peaks[w_max]],"max_peak"] <- 1
        }
       ggplot(data1, aes(x=Time, y = ma_differ)) + geom_line() +
         geom_point(data = data1[peaks,], col = "red")
      }
    }
  }
}

### Add points
ggplot(data_ds, aes(x=Time, y = ma_differ, group = interaction(rep, inoc, glucose))) + 
  geom_line(aes(col = factor(inoc))) + 
  geom_point(data = data_ds %>% filter(peaks>0), aes(shape = factor(rep))) + 
  facet_wrap(inoc~glucose, ncol = 4) + 
  scale_y_continuous("Difference in heatflow per time step")+ 
  scale_x_continuous("Time (h)") + 
  scale_color_discrete("Inoculum", breaks = c(1,2,3,4,5,6,7), labels = c(str2expression("10^1"),str2expression("10^2"),str2expression("10^3"),
                                                                         str2expression("10^4"),str2expression("10^5"),
                                                                         str2expression("10^6"),str2expression("10^7"))) 
ggsave("data_paper2/plots/11016sugar_od_differentiated_byinoc_withpeaks.pdf")
### Look at one
data1 <- data_ds %>% filter(rep == 1, glucose == 5, inoc == 2)
data1 <- data_ds %>% filter(rep == 1, glucose == 1.25, inoc == 7)
ggplot(data1, aes(x=Time, y = ma_differ)) + geom_line() +
  geom_point(data = data1 %>% filter(peaks>0), 
             aes(shape = factor(rep)), col = "red")

ggplot(data1, aes(x=Time, y = ma_differ)) + geom_line() +
  geom_point(data = data1 %>% filter(max_peak>0), 
             aes(shape = factor(rep)), col = "red") + 
  geom_point(data = data1 %>% filter(scnd_peak>0),
             aes(shape = factor(rep)), col = "orange")


### Plot max and second
ggplot(data_ds, aes(x=Time, y = ma_differ, group = interaction(rep, inoc, glucose))) + 
  geom_line(aes(col = factor(inoc))) + 
  geom_point(data = data_ds %>% filter(max_peak>0), 
             aes(shape = factor(rep)), col = "red") + 
  geom_point(data = data_ds %>% filter(scnd_peak>0),
             aes(shape = factor(rep)), col = "orange") + 
  facet_wrap(inoc~glucose, ncol = 4) + 
  scale_y_continuous("Difference in heatflow per time step")+ 
  scale_x_continuous("Time (h)") + 
  scale_color_discrete("Inoculum", breaks = c(1,2,3,4,5,6,7), labels = c(str2expression("10^1"),str2expression("10^2"),str2expression("10^3"),
                                                                         str2expression("10^4"),str2expression("10^5"),
                                                                         str2expression("10^6"),str2expression("10^7"))) 
ggsave("data_paper2/plots/11016sugar_od_differentiated_byinoc_withmax_2nd.pdf")

### Plot height 
ggplot(data_ds %>% filter(scnd_peak > 0), 
       aes(x=inoc, y = ma_differ, group = interaction(rep, glucose))) + 
  geom_point(aes(col = factor(glucose))) +
  geom_line(aes(col = factor(glucose))) + 
  scale_y_continuous("Value of difference in heatflow at second peak") + 
  scale_color_discrete("Glucose conc.\n(mg/L)")
  scale_x_continuous("Inoculum", breaks = c(1,2,3,4,5,6,7), labels = c(str2expression("10^1"),str2expression("10^2"),str2expression("10^3"),
                                                                         str2expression("10^4"),str2expression("10^5"),
                                                                         str2expression("10^6"),str2expression("10^7"))) 
ggsave("data_paper2/plots/11016sugar_od_differentiated_2ndpeak.pdf")

# Average over the reps
mean_2nd <- data_ds %>% filter(scnd_peak > 0) %>% group_by(inoc, glucose) %>%
  summarise(mean = mean(ma_differ), sd = sd(ma_differ))

ggplot(mean_2nd, 
       aes(x=inoc, y = mean, group = interaction(rep, glucose))) + 
  geom_errorbar(aes(ymin = mean - sd/sqrt(3), ymax = mean + sd/sqrt(3), col = factor(glucose))) + 
  geom_point(aes(col = factor(glucose))) +
  geom_line(aes(col = factor(glucose))) + 
  scale_y_continuous("Value of difference in heatflow at second peak", lim = c(0,0.025)) + 
  scale_color_discrete("Glucose conc.\n(mg/L)") + 
scale_x_continuous("Inoculum", breaks = c(1,2,3,4,5,6,7), labels = c(str2expression("10^1"),str2expression("10^2"),str2expression("10^3"),
                                                                     str2expression("10^4"),str2expression("10^5"),
                                                                     str2expression("10^6"),str2expression("10^7"))) 
ggsave("data_paper2/plots/11016sugar_od_differentiated_inocvs2ndpeak.pdf")
