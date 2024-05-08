#### Metadata analysis

###******* LOAD UP LIBRARIES AND DATA NEEDED *************#############################################################
## libraries needed
library(reshape2) # for data manipulation
#library(ggplot2) # for plotting
#library(plyr) # for data manipulation
library(data.table)
library(magrittr)
#library(dplyr)
library(MASS)
library(tidyverse)
library(RColorBrewer)
library(here)
library(patchwork)
library(ggpubr)
library(viridis)
library(gridExtra)
library(mlogit)
library(gmodels)
library(glmmTMB)

library(lme4)
library(nlme)
library(corrplot)
library(glmmTMB)
library(emmeans)
library(performance)
library(MASS)
library(corrplot)
library(psych)
library(DHARMa)
library(VGAM)

theme_set(theme_bw(base_size=12)) # theme setting for plots: black and white (bw) and font size (24)
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(14)

setwd(here::here())

## Load in grofit functions if package no longer working
source("code/grofit_functions.R") # have copied those we use into this R script

## Load in functions for this analysis
source("code/functions_for_heat_curves_additional_double_peak.R")

################################################################################
## Data cleaning ##
################################################################################

#### Read in data ##############################################################
## Read in data on clusters ####
c <- c()
c$parameters <- read_csv("output/clustered_parameters.csv")

## Select clusterdata for each strain. Decided to base metadata analysis of drytime 0, inoc 5 data. First check if all replicates per strain show the same cluster type
# Then deduplicate to keep that one cluster type per strain

mdd_cluster <- c$parameters %>% 
  dplyr::select(inocl, cluster, strain, drytime, rep, width_peak) %>% #or without width_peak
  filter(drytime == 0) %>%
  filter(inocl == 5) %>% 
  filter(!(strain %in% c("Newman", "RWW12", "SA3297", "SA2704", "RWW146", "SAC042W", "Mu50", "M116")))
#length(unique(mdd_cluster$strain)) #97 

mdd_cluster <- mdd_cluster %>%
  dplyr::mutate(cluster = replace_na(cluster, "unclustered"))
#unique(mdd_cluster$cluster)

mdd_cluster %>%
  group_by(strain) %>%
  summarise(u=unique(cluster))%>% 
  group_by(strain) %>%
  summarise(n = n()) %>%
  filter(n > 1)
#All strains have only one cluster for multiple replicates of drytime 0, inoc 5 -> replicates can be filtered

mdd_cluster <- mdd_cluster[!duplicated(mdd_cluster$strain),]
dim(mdd_cluster) # 97 x 5

## Read in metadata ####
mdd <- read_csv2("data/MACOTRA_100metadata.csv")
names(mdd)

## Combine cluster data with metadata ####
mdd_orig <- mdd
mdd <- as.data.frame(mdd)
mdd_cluster <- as.data.frame(mdd_cluster)

mdd$strain <- as.character(mdd$strain)
#anti_join(mdd, mdd_cluster) #data for strain 11287 present in mdd, but not in mdd_cluster
#arrange(mdd_cluster, strain) #indeed, strain 11287 not present in mdd_cluster

mdd <- mdd %>%
  filter(!(strain %in% c("11287"))) # remove strain 11287

# length(unique(mdd$strain))
# length(unique(mdd_cluster$strain))
# dim(mdd) # 97 118
# dim(mdd_cluster) #97  5

mdd_all <- merge(mdd_cluster, mdd, by.x = "strain")
dim(mdd_all)
names(mdd_all)
str(mdd_all)

unique(mdd_all$cluster) # which cluster types in dataset

mdd_all$odd <- 0
mdd_all <- mdd_all %>% 
  mutate(odd = ifelse(cluster != "normal",1,0)) #convert all non-normal cluster types in odd variable as odd = 1
# unique(mdd_all$odd)
# dim(mdd_all)

# Make new variable which combines CC398, CC45, CC80 and Other (ST59) in Other
mdd_all$lineage_short <- mdd_all$lineage
mdd_all$lineage_short[!(mdd_all$lineage_short %in% c("CC1", "CC22", "CC30", "CC5", "CC8"))] <- "Other" 
unique(mdd_all$lineage_short)

# Convert grouping and binary variables to factor, convert numerical variables to numeric #######
mdd_all <- mdd_all %>% mutate(across(2:5, as.factor))
mdd_all <- mdd_all %>% mutate(across(8:125, as.factor)) # select columns 22 to 124 (by index)
str(mdd_all)

## Which data contains NA?
mdd_all %>%
  summarise_each(list(~sum(is.na(.)))) %>%
  gather()  
#NA data in carriage (n=11), MLVA (n=39)/MLST(n=58/78), PVL (n=19), spa_type(n=78)

# ##Subset metadata dataframes
# mdd1 <- dplyr::select(mdd_all, strain, cluster, width_peak : spa_type) #metadata

################################################################################
## Data exploration ############################################################
################################################################################

## Data description ####
mdd_all%>%
  count(cluster)

# cluster  n
# 1        double  5
# 2        normal 56
# 3 post_shoulder 20
# 4         spike  6
# 5   unclustered  2
# 6          wide  8
# sum       total 97

mdd_all%>%
  group_by(lineage, cluster) %>%
    count(cluster, lineage) %>% group_by(lineage) %>% dplyr::summarise(n=n())

table(mdd_all$lineage, mdd_all$cluster) #overview strains for each genetic lineage and cluster
# table_lineage <- table(mdd_all$lineage, mdd_all$cluster)
# png("plots/table_cluster_lineage.png", height=4, width=7, units = "in", res = 72)
# grid.table(table_lineage)
# dev.off()

## Data exploration plots ####
## Main question: are any metabolic clusters correlated with infection and/or success?
## If so, are there also correlations with infection site, AMR/virulence genes and/or genetic lineage?

ggplot(mdd_all) + 
  geom_bar(aes(x = cluster, fill = country), position = "fill") +
  scale_x_discrete("Cluster type", label = c("Double (n=5)","Normal (n=56)","Post-shoulder (n=20)","Spike (n=6)","Unclustered (n=2)","Wide (n=8)")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# more unsuccessful in non-normal clusters?
ggplot(mdd_all) + 
  geom_bar(aes(x = cluster, fill = carriage), position = "fill") +
  scale_x_discrete("Cluster type", label = c("Double (n=5)","Normal (n=56)","Post-shoulder (n=20)","Spike (n=6)","Unclustered (n=2)","Wide (n=8)")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# more infection in non-normal clusters?
# ggplot(mdd_all) + 
#   geom_bar(aes(x = cluster, fill = site), position = "fill") +
#   scale_x_discrete("Cluster type", label = c("Double (n=5)","Normal (n=56)","Post-shoulder (n=20)","Spike (n=6)","Unclustered (n=2)","Wide (n=8)")) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# #correlation with infection/carriage? focus on that first, then try to minimize categories in site
ggplot(mdd_all) + 
  geom_bar(aes(x = cluster, fill = pvl), position = "fill") +
  scale_x_discrete("Cluster type", label = c("Double (n=5)","Normal (n=56)","Post-shoulder (n=20)","Spike (n=6)","Unclustered (n=2)","Wide (n=8)")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# pvl positive in double and spike, related to infection/carriage?
ggplot(mdd_all) + 
  geom_bar(aes(x = cluster, fill = lineage), position = "fill") +
  scale_x_discrete("Cluster type", label = c("Double (n=5)","Normal (n=56)","Post-shoulder (n=20)","Spike (n=6)","Unclustered (n=2)","Wide (n=8)")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# correlation with CC30, CC22 and CC5?
ggplot(mdd_all) + 
  geom_bar(aes(x = odd, fill = lineage_short), position = "fill") + # 
  scale_x_discrete("", label = c("Normal (n=56)","Odd (n=41)")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


################################################################################
## Binary GLMs #################################################################
################################################################################

#Variables to test:
data <-dplyr::select(mdd_all, c(cluster, country, institute, year, success, carriage, site, lineage, pvl, # metadata
                                blaZ, mecA, cat,ermA, ermB, ermC, msrA, fusB, fusC, AAC_6_Ie_APH_2_Ia, ANT_4_Ib, ANT_6_Ia, APH_3_IIIa,
                                aad_6, mupA, tetK, tetM, # AMR genes
                                dfrC, dfrG, dfrK, adsA, aur, cap8A, cap8B, cap8C, cap8D, cap8E, cap8F, cap8G, cap8H, cap8I, cap8J,
                                cap8K, cap8L, cap8M, cap8N, cap8O, cap8P, chp, clfA, clfB, clpP,cna, coa, ebp, esaA, esaB, esaC, 
                                essA, essB, essC, esxA, esxB, fnbA, fnbB, geh, hlb, hld, hlgA, hlgB, hly_hla, hysA, icaA, icaB,
                                icaC, icaD, icaR, isdA, isdB, isdC, isdD, isdE, isdF, isdG, lip, lplA1, lukF_PV, lukS_PV, map, 
                                sak, sbi, scn, sdrC, sdrD, sdrE, sea, seb, sec, sed, see, seh, selk, sell, selq, spa, srtB, sspA,
                                sspB, sspC, tsst_1, vWbp, odd, lineage_short)) 

data_carriage<-na.omit(data)
data_pvl<-na.omit(data)

## Univariate GLMs ############################################################# 
glm1 <- glmmTMB::glmmTMB(odd ~ sea + (1|lineage), family = binomial, data = data)
summary(glm1)

glm1 <- glmmTMB::glmmTMB(odd ~ carriage, family = binomial, data = data_carriage)
summary(glm1)

# p > 0.2 for these variables: lineage, lineage_short, mecA, ermA, ANT_4_Ib, APH_3_IIIa, aad_6, cna, sak, scn, sea
table(data$odd, data$sea) # data$odd as rows
plot(data$odd, data$sea)

## Multivariate
#glm2 <- glmmTMB(odd ~ lineage_short + mecA + ermA + ANT_4_Ib + APH_3_IIIa + aad_6 + cna + sak + scn + sea, family=binomial, data=data)
#glm2 <- glmmTMB(odd ~ lineage_short + mecA + ermA + ANT_4_Ib + cna + sak + scn + sea, family=binomial, data=data)
glm2 <- glmmTMB(odd ~ lineage_short + mecA + ermA + ANT_4_Ib + cna + sak + scn + sea, family=binomial, data=data)
summary(glm2)
glmmTMB:::Anova.glmmTMB(glm2, type="II")
check_collinearity(glm2) # vif>10: remove variable, preferably vif<5
AIC(glm2)

step<-stepAIC(glm2,direction='backward', scope=glm1, k = 2) #to perform model selection based on backwards AIC selection instead of p-values
step

#glm2 <- glmmTMB(odd ~ lineage_short + mecA + ermA + scn + sea, family=binomial, data=data)
glm2 <- glmmTMB(odd ~ mecA + ermA + ANT_4_Ib + sea, family=binomial, data=data)
glmmTMB:::Anova.glmmTMB(glm2, type="II")
summary(glm2)

################################################################################
## Multinomial GLMs ############################################################
################################################################################

# # Paramatric test: multinomial logistic regression > not enough data for this
# modeldata <- mlogit.data(mdd2, choice = "cluster", shape = "wide")
# 
# m0 <- mlogit(cluster ~ 1, data = mdd2, reflevel = "normal") #baseline model
# summary(m0)
# m1 <- mlogit(cluster ~ 1 | lineage, data = modeldata, reflevel = "normal")
# summary(m1)
# # model with success, carriage, site or pvl does not explain variability in the data
# # model with lineage does explain variability in the data > attributed to CC30 in cluster wide
# 
# ## Chi-square test
# 
# fisher.test(mdd1$cluster, mdd1$site, simulate.p.value=TRUE) #large dataset, simulated p-value for lineage and site
# #fisher.test(mdd2$cluster, mdd2$lineage, simulate.p.value=TRUE)
# 
# # CrossTable(mdd1$lineage, mdd1$cluster)
# # chisq.test(mdd1$lineage, mdd1$cluster)
# # fisher.test(mdd1$cluster, mdd1$lineage, simulate.p.value=TRUE)
# # CrossTable(mdd2$lineage, mdd2$cluster)
# # chisq.test(mdd2$lineage, mdd2$cluster)
# # fisher.test(mdd2$cluster, mdd2$country, simulate.p.value=TRUE)
# 
# table(mdd2$lineage, mdd2$cluster)
# table(mdd3$lineage, mdd3$cluster)
# 
# #CrossTable(mdd3$cluster, mdd3$lineage)
# #chisq.test(mdd3$cluster, mdd3$success)
# fisher.test(mdd3$cluster, mdd3$success, simulate.p.value=FALSE)
# 
# #chisq.test(mdd3$cluster, mdd3$carriage)
# fisher.test(mdd3$cluster, mdd3$carriage, simulate.p.value=FALSE)
# 
# chisq.test(mdd3$cluster, mdd3$lineage)
# fisher.test(mdd3$cluster, mdd3$lineage, simulate.p.value=TRUE)
# 
# # library(rstatix)
# # fisher_test(table(mdd2$lineage, mdd2$cluster), simulate.p.value=TRUE)
# # pairwise_fisher_test(as.matrix(table(mdd2$lineage, mdd2$cluster)), p.adjust.method = "fdr")
# # how to convert table to two-dimensional contingency table? maybe not working as 2xc comparison, too many categories?
# 
# ggplot(mdd3 , aes(cluster, width_peak)) +
#   geom_boxplot(aes(colour = lineage)) +
#   geom_point(aes(colour = lineage)) +
#   facet_wrap(~lineage)
# ggsave("plots/metadata_width_peakCC30.png", width = 10, height = 6)
