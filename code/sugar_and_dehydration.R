#### Sugar and dehydration
library(tidyverse)
library(here)
library(patchwork)
library(MESS) # for auc
theme_set(theme_bw(base_size = 14))

setwd(here())
data1 <- read_csv("data/glucose_HF_noBC_t8.csv")[,-1]
data2 <- read_csv("data/glucose_HF_t8.csv")[,-1]
data0 <- read_csv("data/glucose_HF.csv")

data_dry <- data1
data_dry$value_base <- data2$value

### Raw data
ggplot(data_dry, aes(x= Time, y = value, group = interaction(rep, glucose_conc, inoc_name))) + geom_line(aes(colour = factor(inoc))) + 
  facet_wrap(~glucose) + 
  scale_y_continuous("Heat flow") + 
  scale_color_discrete("Inoculum")
ggsave("plots/glucose_drying_raw_data.pdf")

# baseline correction looks better? 
ggplot(data_dry, aes(x= Time, y = value_base, group = interaction(rep, glucose_conc, inoc_name))) + geom_line(aes(colour = factor(inoc))) + 
  facet_wrap(~glucose) + 
  scale_y_continuous("Heat flow") + 
  scale_color_discrete("Inoculum")
ggsave("plots/glucose_drying_raw_data_base_correct.pdf")


########## MAX HF

# Filter by HF as don't care about input values - same nongrowth for hf or od
max_hf_tab <- data_dry  %>% group_by(rep, inoc, glucose) %>% mutate(max_hf = max(value_base)) %>%  
  filter(value_base == max_hf)

g1 <- ggplot(max_hf_tab, aes(x=inoc, y = value_base,group = interaction(glucose, rep))) + geom_point(aes(col = factor(glucose))) + 
  geom_line(aes(col = factor(glucose))) + 
  scale_y_continuous("Maximum HF") + scale_x_continuous("Inoculum") +
  scale_color_discrete("Glucose")
ggsave("plots/glucose_dry_inoc_vs_max_hf.pdf")

g2 <- ggplot(max_hf_tab, aes(x=inoc, y = Time, group = interaction(glucose, rep))) + geom_point(aes(col = factor(glucose))) +  
  geom_line(aes(col = factor(glucose))) + 
  scale_y_continuous("Time to maximum HF") + scale_x_continuous("Inoculum") + 
  scale_color_discrete("Glucose")
ggsave("plots/glucose_dry_time_vs_max_hf.pdf")

g1 + g2 + plot_layout(guides = "collect")
ggsave("plots/glucose_dry_both_time&inoc_vs_max_hf.pdf", width = 10, height = 5)

g1 <- ggplot(max_hf_tab, aes(x=inoc, y = value_base,group = interaction(glucose))) + geom_point(aes(col = factor(glucose))) + 
  geom_smooth(aes(col = factor(glucose), fill = factor(glucose)),method='lm', formula= y~x) + 
  scale_y_continuous("Maximum HF") + scale_x_continuous("Inoculum") + 
  facet_wrap(~glucose, ncol =  4)+ 
  scale_color_discrete("Glucose")+ 
  scale_fill_discrete("Glucose")
ggsave("plots/glucose_dryinoc_vs_max_HF_smooth.pdf")

g2 <- ggplot(max_hf_tab, aes(x=inoc, y = Time, group = interaction(glucose))) + geom_point(aes(col = factor(glucose))) +  
  geom_smooth(aes(col = factor(glucose), fill = factor(glucose)),method='lm', formula= y~x) + 
  scale_y_continuous("Time to maximum HF") + scale_x_continuous("Inoculum") + 
  facet_wrap(~glucose, ncol =  4)+ 
  scale_color_discrete("Glucose")+ 
  scale_fill_discrete("Glucose")
ggsave("plots/glucose_dryinoc_vs_time_max_HF_smooth.pdf")

g1 / g2 + plot_layout(guides = "collect")
ggsave("plots/glucose_dryinoc_vs_time&val_max_HF_smooth.pdf", width = 10, height = 5)

#### Heat map 
# Take the mean over the replicates
max_hf_tab_mean <- max_hf_tab %>% dplyr::select(glucose, inoc, rep, Time, max_hf) %>% group_by(glucose, inoc) %>% 
  summarise(mean_time = mean(Time), mean_max = mean(max_hf))

max_hf_tab_mean$glucose <- as.character(max_hf_tab_mean$glucose)
g1 <- ggplot(max_nongrowth_tab_mean, aes(inoc, glucose, fill= mean_time)) + geom_tile() +
  scale_fill_distiller(palette = "RdBu","Time to\nmax. HF", ) + 
  scale_x_continuous("Inoculum") + 
  scale_y_discrete("Glucose concentration")

g2 <- ggplot(max_nongrowth_tab_mean, aes(inoc, glucose, fill= mean_max)) + geom_tile() +
  scale_fill_distiller(palette = "PRGn","Value of\nmax. HF") + 
  scale_x_continuous("Inoculum") + 
  scale_y_discrete("Glucose concentration")

g1 + g2 
ggsave("plots/glucose_dry_inoc_vs_time&val_max_heatmap.pdf", width = 15, height = 5)

#### MIN NON_GROWTH

# Filter by HF as don't care about input values - same nongrowth for hf or od
min_hf_tab <- data_dry %>% group_by(rep, inoc, glucose) %>% mutate(min_non = min(value_base)) %>%  
  filter(value_base == min_non)

g1 <- ggplot(min_hf_tab, aes(x=inoc, y = value_base,group = interaction(glucose, rep))) + geom_point(aes(col = factor(glucose))) + 
  geom_line(aes(col = factor(glucose))) + 
  scale_y_continuous("Minimum non-growth") + scale_x_continuous("Inoculum") +
  scale_color_discrete("Glucose")
ggsave("plots/glucose_dry_inoc_vs_min_HF_.pdf")

g2 <- ggplot(min_hf_tab, aes(x=inoc, y = Time, group = interaction(glucose, rep))) + geom_point(aes(col = factor(glucose))) +  
  geom_line(aes(col = factor(glucose))) + 
  scale_y_continuous("Time to minimum non-growth") + scale_x_continuous("Inoculum") + 
  scale_color_discrete("Glucose")
ggsave("plots/glucose_dry_inoc_vs_time_min_HF.pdf")

g1 + g2 + plot_layout(guides = "collect")
ggsave("plots/glucose_dry_inoc_vs_time&val_min_HF.pdf", width = 10, height = 5)

g1 <- ggplot(min_hf_tab, aes(x=inoc, y = value_base,group = interaction(glucose))) + geom_point(aes(col = factor(glucose))) + 
  geom_smooth(aes(col = factor(glucose), fill = factor(glucose)),method='lm', formula= y~x) + 
  scale_y_continuous("Minimum non-growth") + scale_x_continuous("Inoculum") + 
  facet_wrap(~glucose, ncol =  4)+ 
  scale_color_discrete("Glucose")+ 
  scale_fill_discrete("Glucose")
ggsave("plots/glucose_dry_inoc_vs_min_HF_smooth.pdf")

g2 <- ggplot(min_hf_tab, aes(x=inoc, y = Time, group = interaction(glucose))) + geom_point(aes(col = factor(glucose))) +  
  geom_smooth(aes(col = factor(glucose), fill = factor(glucose)),method='lm', formula= y~x) + 
  scale_y_continuous("Time to minimum non-growth") + scale_x_continuous("Inoculum") + 
  facet_wrap(~glucose, ncol =  4)+ 
  scale_color_discrete("Glucose")+ 
  scale_fill_discrete("Glucose")
ggsave("plots/glucose_dry_inoc_vs_time_min_HF_smooth.pdf")

g1 / g2 + plot_layout(guides = "collect")
ggsave("plots/glucose_dry_inoc_vs_time&val_min_HF_smooth.pdf", width = 10, height = 5)

#### Heat map 
# Take the mean over the replicates
min_hf_tab_mean <- min_hf_tab %>% dplyr::select(glucose, inoc, rep, Time, min_non) %>% group_by(glucose, inoc) %>% 
  summarise(mean_time = mean(Time), mean_min = mean(min_non))

min_hf_tab_mean$glucose <- as.character(min_hf_tab_mean$glucose)
g1 <- ggplot(min_hf_tab_mean, aes(inoc, glucose, fill= mean_time)) + geom_tile() +
  scale_fill_distiller(palette = "RdBu","Time to\nmin. non-growth", ) + 
  scale_x_continuous("Inoculum") + 
  scale_y_discrete("Glucose concentration")

g2 <- ggplot(min_hf_tab_mean, aes(inoc, glucose, fill= mean_min)) + geom_tile() +
  scale_fill_distiller(palette = "PRGn","Value of\nmin. non-growth") + 
  scale_x_continuous("Inoculum") + 
  scale_y_discrete("Glucose concentration")

g1 + g2 
ggsave("plots/glucose_dry_inoc_vs_time&val_min_heatmap.pdf", width = 15, height = 5)



#### Extract AUC 
## What are the glucose concentrations?
u <- unique(data_dry$glucose)
## How many replicates? 
r <- unique(data_dry$rep) # replicates
# What are the inoculums? 
q <- unique(data_dry$inoc)

# Where the parameters for each strain are stored
param <- matrix(0, length(u)*length(q)*length(r), 9); 
index <- 1 # for counting 


## Run thru each 
for(jj in 1:length(u)){ # for each glucose
  for(ii in 1:length(r)){ # for each replicate
    for(ll in 1:length(q)){ #each of the inocula
      
      data <- data_dry %>% 
        filter(glucose == u[jj], rep == r[ii], inoc == q[ll]) %>% filter(Time > 3, Time < 18)
      data1 <- data_dry %>% 
        filter(glucose == u[jj], rep == r[ii], inoc == q[ll]) %>% filter(Time > 3, Time < 15) # BEST ONE 
      data2 <- data_dry %>% 
        filter(glucose == u[jj], rep == r[ii], inoc == q[ll]) %>% filter(Time > 5, Time < 22)
      
      
      if(dim(data)[1] > 0){ # if this replicate exists for this strain (i.e. there is data)
        
        print(c(u[jj], r[ii], q[ll])) # output so can track how it is working
        p <-0; p_all <- 0;p1 <-0; p_all1 <- 0;p2 <-0; p_all2 <- 0;
        try(p <- auc(data$Time, data$value_base, type = "spline", subdivisions = 100000))
        try(p_all <- auc(data$Time, data$value_base, type = "spline", subdivisions = 10000,absolutearea = TRUE))
        try(p1 <- auc(data1$Time, data1$value_base, type = "spline", subdivisions = 100000))
        try(p_all1 <- auc(data1$Time, data1$value_base, type = "spline", subdivisions = 10000,absolutearea = TRUE))
        try(p2 <- auc(data2$Time, data2$value_base, type = "spline", subdivisions = 100000))
        try(p_all2 <- auc(data2$Time, data2$value_base, type = "spline", subdivisions = 10000,absolutearea = TRUE))
        
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
param$auc2 <- as.numeric(param$auc2)
param$auc2_all <- as.numeric(param$auc2_all)
param$inoc <- as.numeric(param$inoc)

g1 <- ggplot(param, aes(x=inoc, y = auc2, group = interaction(rep,glucose),col = glucose)) + geom_point() + geom_line(aes()) + 
  scale_x_continuous("Inoculum") + 
  scale_y_continuous("AUC")
ggsave("plots/glucose_dry_auc.pdf", width = 10, height = 5)

g2 <- ggplot(param, aes(x=inoc, y = auc2_all, group = interaction(rep,glucose),col = glucose)) + geom_point() + geom_line(aes()) + 
  scale_x_continuous("Inoculum") + 
  scale_y_continuous("AUC absoulte")
ggsave("plots/glucose_dry_abs_auc.pdf", width = 10, height = 5)

g1 + g2 + plot_layout(guides = "collect")
ggsave("plots/glucose_dry_auc_both_raw.pdf", width = 15, height = 5)

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
ggsave("plots/glucose_dry_auc_both_fit.pdf", width = 15, height = 5)


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
ggsave("plots/glucose_dry_auc_abs_heatmap.pdf", width = 15, height = 5)







