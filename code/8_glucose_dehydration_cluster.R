#### Clustering of glucose dehydration data


###******* LOAD UP LIBRARIES AND DATA NEEDED *************#############################################################
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
library(patchwork)
library(ggpubr)
library(viridis)
library(gridExtra)

theme_set(theme_bw(base_size=12)) # theme setting for plots: black and white (bw) and font size (24)
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(14)

setwd(here::here())

## Load in grofit functions if package no longer working
source("code/grofit_functions.R") # have copied those we use into this R script

## Load in functions for this analysis
source("code/functions_for_heat_curves_additional_double_peak.R")

####### Input data
data1 <- read_csv("data/glucose_HF_noBC_t8.csv")[,-1] # without baseline correction 
data2 <- read_csv("data/glucose_HF_t8.csv")[,-1] # with baseline correction 
data01 <- read_csv("data/glucose_HF.csv")[,-1] # time 0 data 
data02 <- read_csv("data/glucose_HF_noBC.csv")[,-1] # time 0 data without baseline correction

data01$drytime <- 0
data02$drytime <- 0
data1$drytime <- 7
data2$drytime <- 7

ddm <- rbind(data1, data02) # only WITHOUT baseline correction
w1 <- which(ddm$rep == "x")
w2 <- which(ddm$rep == "xi")
w3 <- which(ddm$rep == "xii")
ddm$rep <- as.numeric(ddm$rep)
ddm[w1,"rep"] <- 1
ddm[w2,"rep"] <- 2
ddm[w3,"rep"] <- 3

###*** Clean those that have v high low early data (remove for this one)
ttt <- intersect(intersect(intersect(intersect(which(ddm$glucose == 1.25), which(ddm$rep == 2)), which(ddm$glucose_conc == "D")),which(ddm$inoc == 2)),which(ddm$drytime == 7))
sss <- intersect(intersect(intersect(intersect(which(ddm$glucose == 5), which(ddm$rep == 1)), which(ddm$glucose_conc == "B")),which(ddm$inoc == 2)),which(ddm$drytime == 7))
ppp <- intersect(intersect(intersect(intersect(which(ddm$glucose == 5), which(ddm$rep == 1)), which(ddm$glucose_conc == "B")),which(ddm$inoc == 1)),which(ddm$drytime == 7))
qqq <- intersect(intersect(intersect(intersect(which(ddm$glucose == 5), which(ddm$rep == 3)), which(ddm$glucose_conc == "B")),which(ddm$inoc == 2)),which(ddm$drytime == 7))
ddd <- intersect(intersect(intersect(intersect(which(ddm$glucose == 5), which(ddm$rep == 2)), which(ddm$glucose_conc == "B")),which(ddm$inoc == 1)),which(ddm$drytime == 7))
kkk <- intersect(intersect(intersect(intersect(which(ddm$glucose == 0), which(ddm$rep == 1)), which(ddm$glucose_conc == "E")),which(ddm$inoc == 6)),which(ddm$drytime == 0))
ggg <- intersect(intersect(intersect(intersect(which(ddm$glucose == 0), which(ddm$rep == 3)), which(ddm$glucose_conc == "E")),which(ddm$inoc == 7)),which(ddm$drytime == 0))
hhh <- intersect(intersect(intersect(intersect(which(ddm$glucose == 0), which(ddm$rep == 2)), which(ddm$glucose_conc == "E")),which(ddm$inoc == 6)),which(ddm$drytime == 0))
ddm <- ddm[-c(sss,ttt,ppp,qqq,ddd,kkk,ggg,hhh),]


###******* MODEL FITTING *************#######################################################################################################################################

## What are the glucose concentrations?
u <- as.character(unique(ddm$glucose_conc,stringsasFactors = FALSE)) # strains
## How many replicates? 
r <- unique(ddm$rep) # replicates
# What are the inoculums? 
q <- unique(ddm$inoc)


# Where the parameters for each strain are stored
param <- matrix(0, length(u)*length(r)*length(q)*3, 51); # number of strains x number of replicates x number of experimental conditions
index <- 1 # for counting 
max <- c() # for calibration
drying_times <- unique(ddm$drytime)

## Run thru each strain/rep/drying time/ inoculum: fit a growth curve to each
for(jj in 1:length(u)){ # for each strain
  r <- unlist(unique(ddm %>% filter(glucose_conc == u[jj]) %>% dplyr::select(rep))[,1]) # which replicates for this strain (all have different names)
  if(length(r) < 3){print(paste0(u[jj], " has too few replicates"))}else{ # filter out if only 2 replicates
    for(ii in 1:length(r)){ # for each replicate
      for(kk in c(1,2)){ #each of the experimental conditions
        for(ll in 1:length(q)){ #each of the inocula
          
          conc <- u[jj];
          replicate <- unlist(r)[ii]
          condition <- drying_times[kk]
          inocl <- q[ll]
          
          wi <- intersect(which(ddm$glucose_conc == conc),which(ddm$rep == replicate)) # if fit to each replicate
          wj <- intersect(wi, which(ddm$drytime == condition))
          w <- intersect(wj, which(ddm$inoc == as.numeric(inocl)))
          
          if(length(w) > 0){ # if this replicate exists for this strain (i.e. there is data)
            data1 <- ddm[w,] # Grab data
            
            print(c(jj, conc, replicate, condition, inocl)) # output so can track how it is working
            p <- cut_extract_dp(data1, "Time", "value", paste(conc, replicate, condition, inocl,sep="_")) ### NEW function: runs on any timeseries: gives baseline and cut parameter analysis
            
            ## Required parameters
            
            param[index,] <- c(conc, replicate, condition, inocl, p$param)
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
colnames(param) <- c("strain","rep","drytime","inoc",
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
dim(param)

## Store so don't have to run above
write.csv(param, "output/param_multiple_peaks_glucose_dehydration.csv")

########## ************************************************************************************************ #########################
##### Clustering #######
# For now have "Strain" for "Glucose concentration" 
ddm <- rename(ddm, "strain" = "glucose_conc")
ddm <- rename(ddm, "value_J" = "value")
c <- cluster(ddm, param,plot_where = "plots/glucose_conc_") 
########## ************************************************************************************************ #########################

### Store output of clustering
write.csv(c$parameters,"output/clustered_parameters_glucose_dehydration.csv")
write.csv(c$ts,"output/clustered_time_series_glucose_dehydration.csv")

## Read in if doing later
c <- c()
c$parameters <- read_csv("output/clustered_parameters_glucose_dehydration.csv")
c$ts <- read.csv("output/clustered_time_series_glucose_dehydration.csv")

### Overview by glucose
ggplot(c$ts, aes(x=Time, y = value_J, group = interaction(strain, inoc, rep, cluster))) + 
  geom_line(aes(colour = cluster)) + 
  facet_grid(inoc + drytime~glucose)

ggplot(c$parameters, aes(x=inoc, group = cluster)) + 
  geom_bar(stat = "count", position = "dodge", aes(fill = cluster)) + 
  facet_grid(~drytime) + 
  scale_x_continuous("Inoculum") + 
  scale_y_continuous("Number of datasets") + 
  scale_fill_discrete("Cluster type")
ggsave("plots/glucose_cluster_drytime.jpeg")


### STORY PIC
para <- c$parameters %>% ungroup() %>% mutate(second_peak_h = ifelse((t_m_h_flow > (timepeak + 2)), v_m_h_flow, ifelse(mp_t2 > (timepeak + 4), mp_h2, 0))) 

cluster_data <- para %>% dplyr::select(c("cluster", "inoc","drytime","valpeak", "timepeak", "auc", "exp_gr", "shoulder_point_v","shoulder_point_t",
                                         "shoulder_point_past_v","shoulder_point_past_t", "second_peak_h")) %>% 
  pivot_longer(cols = c(valpeak:second_peak_h)) 

my_comparisons <- list( c("post_shoulder", "double"), c("double", "normal"), c("normal", "post_shoulder"),c("normal","spike"),c("spike","wide"),
                        c("double","spike"),c("double","wide"))

cluster_data$cluster <- factor(cluster_data$cluster, levels = c("normal","double","spike","post_shoulder","wide","unclustered"))
var_name <- c(
  auc = "AUC",
  valpeak = "Max value 1st peak",
  exp_gr = "Max exp growth rate",
  second_peak_h = "Max value 2nd peak"
)

ggplot(cluster_data %>% filter(inoc == 5) %>% filter(name %in% c("auc", "valpeak","exp_gr","second_peak_h")) %>% 
         mutate(across(name, factor, levels=c("auc", "valpeak","exp_gr","second_peak_h"))), aes(x = cluster, y = value)) + 
  geom_boxplot(aes(fill = cluster)) +  facet_wrap(drytime~name,ncol = 4, scales = "free", labeller = labeller(name = var_name))
ggsave("plots/glucose_inoc_5_story_boxplot.jpeg", width = 20, height = 8)

ggplot(cluster_data %>% filter(name %in% c("auc", "valpeak","exp_gr","second_peak_h")) %>% 
         mutate(across(name, factor, levels=c("auc", "valpeak","exp_gr","second_peak_h"))), aes(x = interaction(inoc,cluster), y = value)) + 
  geom_boxplot(aes(fill = cluster)) +  facet_wrap(drytime~name,ncol = 4, scales = "free", labeller = labeller(name = var_name))

ggplot(c$ts %>% filter(cluster == ""), aes(x=Time,y =value_J,group = interaction(rep, inoc, strain))) + 
  geom_line(aes(col = factor(rep))) + 
  facet_wrap(~glucose, ncol = 4)



# 
# 
# gf1 <- c$parameters %>% 
#   filter(drytime == 0) %>% 
#   filter(!(strain %in% c("Newman", "RWW12", "SA3297", "SA2704", "RWW146", "SAC042W", "Mu50", "M116"))) %>%
#   dplyr::select(inoc, cluster, strain) %>% 
#   group_by(inoc, strain) %>% summarise(u = unique(cluster)) %>% group_by(inoc, u) %>% dplyr::summarise(n=n())
# length(unique(gf1$strain)) # CHECK now 97 strains
# # Add in unclustered name
# w<-which(gf1$u == "")
# gf1[w,"u"] <- "unclustered"
# 
# table_cluster_distribution <- gf1 %>% pivot_wider(names_from = u, values_from = n, values_fill = 0) %>% arrange(desc(inocl))
# rowSums(table_cluster_distribution) # CHECK: 97 strains + inocl column = 100/101/102
# table_cluster_distribution <- rename(table_cluster_distribution, Inoculum = inocl, Double = double, 
#                                      Normal = normal, Spike = spike, Wide = wide, Postshoulder = post_shoulder, Unclustered = "NA")
# table_cluster_distribution <- table_cluster_distribution[,c("Inoculum","Normal","Double", "Spike","Postshoulder","Wide", "Unclustered")]
# rownames(table_cluster_distribution) <- NULL
# 
# pdf("plots/table_cluster_distribution.pdf", height=11, width=8.5)
# grid.table(table_cluster_distribution)
# dev.off()
# # png("plots/final/table_cluster_distribution.png", height=2, width=7, units = "in", res = 72)
# # grid.table(table_cluster_distribution)
# # dev.off()
# 
# ### Table with how they change across inoculum 
# table_changes <- c$parameters %>% ungroup() %>% dplyr::select(strain, rep, drytime, inocl, cluster) %>% 
#   pivot_wider(names_from = inocl, values_from = cluster)
# write.csv(table_changes, "output/table_changes_cluster_over_inoc.csv")
# 
# colnames(table_changes) <- c("strain", "rep", "drytime", "five", "four", "three")
# t_changes <- table_changes %>% mutate(fivetofour = ifelse(five == four,0,1),
#                                       fourtothree = ifelse(four == three,0,1),
#                                       difference = fivetofour + fourtothree) %>% filter(difference > 0)
# write.csv(t_changes, "output/table_changes_cluster_over_inoc_with_differences.csv")
# 
# tt_changes <- table_changes %>% mutate(fivetofour = ifelse(five == four,0,1),
#                                        fourtothree = ifelse(four == three,0,1),
#                                        fivetothree = ifelse(five == three,0,1),
#                                        difference = fivetofour + fourtothree + fivetothree) #%>% filter(difference > 0)
# write.csv(tt_changes, "output/tt_changes_cluster_over_inoc_with_more_differences.csv")
# #write.csv2(tt_changes, "output/table_changes_cluster_over_inoc_with_more_differences.csv")
# 
# td_changes <- table_changes %>% mutate(fivetofour = ifelse(five == four,0,1),
#                                        fourtothree = ifelse(four == three,0,1),
#                                        fivetothree = ifelse(five == three,0,1),
#                                        difference = fivetofour + fourtothree + fivetothree) %>% filter(difference > 0)
# dim(td_changes) #259 x 10
# td_changes <- td_changes %>%  filter(drytime == 0) # 129 x 10
# td_changes <- td_changes[!td_changes$five == "",] # 124 x 10
# td_changes <- td_changes[!td_changes$three == "",] # 118 x 10
# 
# write.csv(td_changes, "output/td_changes_cluster_over_inoc_with_more_differences.csv")
# 
# unique(td_changes$five)
# unique(td_changes$three)
# unique(td_changes[, c("five", "three")]) # 15 unique combinations of clusters between five and three
# 
# count_td_changes <- dplyr::count_(td_changes, vars = c('five','three'))
# write.csv(count_td_changes, "output/count_td_changes.csv")
# 
# #### FIGURE: Patterns across inoc
# gf1$u <- factor(gf1$u, levels = c("normal","double","spike","post_shoulder","wide","unclustered"))
# ggplot(gf1, aes(u, inocl, fill= n)) + geom_tile() + scale_fill_viridis(discrete=FALSE, "Number of\nstrains") + 
#   scale_x_discrete("Cluster type", labels = c("Normal","Double","Spike","Post shoulder","Wide","Unclustered")) + 
#   scale_y_continuous("Inoculum size") 
# ggsave("plots/heatmap_cluster_by_inoculum.pdf")
# 
# # ggplot(gf1, aes(u, inocl, fill= n)) + geom_tile() + scale_fill_viridis(discrete=FALSE, "Number of\nstrains") +
# #   scale_x_discrete("Cluster type", label = c("Normal","Double","Spike","Post-shoulder","Wide","Unclustered")) +
# #   scale_y_continuous("Inoculum size")
# # ggsave("plots/final/figure2.png", width = 8, height = 4)
# # ggsave("plots/final/figure2.tiff", width = 8, height = 6, dpi = 600)
