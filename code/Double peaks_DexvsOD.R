##Double peaks Dex vs OD##
##Growth curve analysis of MRSA strains 11006 and 11016 with double peaks##
##based on Gwen Knight's survival model script##

##Libraries##

library(reshape2) # for data manipulation
library(ggplot2) # for plotting
library(grofit) # for fitting growth curves
library(plyr) # for data manipulation
library(data.table)
library(reshape2)
library(magrittr)
library(dplyr)
library(MASS)
library(grid)
library(gridExtra)
library(here)
theme_set(theme_bw(base_size=24)) # theme setting for plots: black and white (bw) and font size (24)

##Home
#home <- "D:/ErasmusMC/Proteomics/Analyse dubbele pieken"
#home <- "C:/Users/Val?rie/Desktop/Schijf_EMC"
#setwd(home)
#getwd()
setwd(here())

## Load in grofit functions if package no longer working
source("grofit_functions.R")

## Load in functions for this analysis
source("20_functions_for_heat_curves.R")

##Organize OD data ####
dd <- read.table("data_paper2/ODvsDextrose_exp2_20210715.txt",header = TRUE)
dd$Time <- dd$Time / 3600
ddm <-reshape2::melt(dd, id.vars=c("Time", "rep"))

ddm$glucose <- 0; ddm$inoc <- 0
ddm[,c("inoc_name", "glucose_conc")]<-colsplit(ddm$variable, "(?<=\\p{L})(?=[\\d+$])", c("char", "digit"))

# inoculum size
w<- which(ddm$inoc_name == "A")
ddm[w,"inoc"] = 7
w<- which(ddm$inoc_name == "B")
ddm[w,"inoc"] = 6
w<- which(ddm$inoc_name == "C")
ddm[w,"inoc"] = 5
w<- which(ddm$inoc_name == "D")
ddm[w,"inoc"] = 4
w<- which(ddm$inoc_name == "E")
ddm[w,"inoc"] = 3
w<- which(ddm$inoc_name == "F")
ddm[w,"inoc"] = 2
w<- which(ddm$inoc_name == "G")
ddm[w,"inoc"] = 1
# glucose concentrations  (g/L)
w<- which(ddm$glucose_conc == "1")
ddm[w,"glucose"] = "5"
w<- which(ddm$glucose_conc == "2")
ddm[w,"glucose"] = "2.5"
w<- which(ddm$glucose_conc == "3")
ddm[w,"glucose"] = "1.25"
w<- which(ddm$glucose_conc == "4")
ddm[w,"glucose"] = "0"


ddm$value <- ddm$value - 0.094 ##medium blank constant

# Check plots
ddm$inoc <- as.factor(ddm$inoc)

ggplot(ddm, aes(x=Time,y = value)) + geom_line(aes(col = inoc)) + facet_wrap(rep~glucose)

### save data
write.csv(ddm, "data/11016_glucose.csv")
