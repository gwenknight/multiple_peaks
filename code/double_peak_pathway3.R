#### Pathway 3

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

## CS data
ddm_cs_orig <- read.csv("data_paper2/ddm_CS.csv")[,-1]
ddm_cs_orig$strain_name <- substr(ddm_cs_orig$strain, 1, 5)

## Times same? seems so
ddm_wide <- ddm_cs_orig %>% dplyr::select(strain_name, inoc, rep, Time, value_J) #%>% pivot_wider(id_cols = c("inoc", "rep", "Time"),names_from = strain_name, values_from = value_J)

ddm_wide_min <- ddm_wide %>% filter(inoc == 1) %>% plyr::rename(replace = c("value_J" = "value_inoc1"))


ddm_all <- left_join(ddm_wide, ddm_wide_min, by = c("rep","Time","strain_name")) %>% 
  dplyr::select(-inoc.y) %>% 
  mutate(non_growth = value_J - value_inoc1, growth = value_inoc1) %>% 
  pivot_longer(cols = c(non_growth,growth)) %>% 
  plyr::rename(replace = c("inoc.x" = "inocl"))

ggplot(ddm_all, aes(x=Time, y = value)) + 
  geom_line(aes(col = factor(name))) + 
  facet_grid(inocl~rep + strain_name) + 
  scale_y_continuous("CS") +
  scale_x_continuous("Time (h)") + 
  scale_color_discrete("Type of heat output") 

ggsave("data_paper2/plots/path_3_subtracted_simple.pdf")  ### NEED TO scale time too...


### Put inoc 1 at end of lag time... 
param_cs <- read.csv("data_paper2/param_orig.csv")[,-1]
param_cs$strain_name <- substr(param_cs$strain_name, 1, 5)
param_cs$inoc <- param_cs$inoc

ddm_all_min <- left_join(ddm_wide, param_cs[,c("strain_name","rep","inoc","lag")]) %>%
  #mutate(lag1 = ifelse(inocl == 1, lag, 0)) %>% 
  #group_by(strain_name, inoc, rep) %>% 
  #mutate(lag_inoc1 = max(lag1) ) %>% 
  #ungroup() %>%
  #mutate(time_diff = Time  - (lag - lag_inoc1)) ### is this correct? 
  mutate(time_diff = Time - lag) # just put all to start at end of lag

### Need to interpolate time...
strains <- unique(ddm_all_min$strain_name)
inocs <- unique(ddm_all_min$inoc)
reps <- unique(ddm_all_min$rep)

interp_data <- c()

minn <- min(ddm_all_min$time_diff)
maxx <- max(ddm_all_min$time_diff)

for(i in 1:length(strains)){
  for(k in 1:length(inocs)){
    for(l in 1:length(reps)){
      
      this_data <- ddm_all_min %>% filter(strain_name == strains[i], inoc == inocs[k], rep == reps[l])
      if(dim(this_data)[1]>0){
        new_data <- approx(this_data$time_diff, this_data$value_J, xout = seq(round(minn,0),round(maxx,0),0.1))
        w <- which(!is.na(new_data$y))
        interp_data <- rbind(interp_data, 
                             cbind(strains[i], inocs[k], reps[l], new_data$x[w],new_data$y[w]))
      }
      
    }
  }
}

interp_data <- as.data.frame(interp_data)
colnames(interp_data) <- c("strain_name","inoc","rep","intrp_time","intrp_value")
interp_data$intrp_value <- as.numeric(interp_data$intrp_value)
interp_data$intrp_time <- as.numeric(interp_data$intrp_time)

interp1 <- interp_data %>% filter(inoc == 1)
interp1$inoc1val = interp1$intrp_value

interp_data2 <- left_join(interp_data, interp1 %>% dplyr::select(-intrp_value,-inoc), by = c("strain_name","rep","intrp_time")) %>%
  mutate(non_growth = intrp_value - inoc1val, growth = inoc1val, total = intrp_value) %>% 
  dplyr::select(strain_name, inoc, rep, intrp_time, non_growth, growth, total) %>% 
  pivot_longer(cols = c(non_growth,growth, total)) 

ggplot(interp_data2, aes(x=intrp_time, y = value, group = interaction(inoc, rep, strain_name, name))) + 
  geom_line(aes(col = factor(name))) + 
  facet_grid(inoc~rep + strain_name) + 
  scale_y_continuous("CS") +
  scale_x_continuous("Time (h)") + 
  scale_color_discrete("Type of heat output") 

ggsave("data_paper2/plots/path_3_subtracted_timeadj.pdf")  ### NEED TO scale time too...

