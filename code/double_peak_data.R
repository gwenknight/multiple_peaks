##Double peaks##
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


## Load in grofit functions if package no longer working
source("code/grofit_functions.R")

##Organize CS heat flow data - in muW ####
dd <- read.table("data_paper2/CS_11006_11016_T0_internalbaseline.txt",header = TRUE)
dd$Time <- dd$Time / 3600
ddm <-reshape2::melt(dd, id.vars=c("Time"))

# Cumulative heat curve in joules
# Heat flow in muW. 1 W = 1 Joule / Second. 
# Assume that heat curve output is x muW over that interval of 1/6hrs
intervals = (1/4)/60/60 # time interval between readings in seconds (every 10 min)
ddm$value_J <-  ddm$value * intervals # convert to Joules (W = J/s => J = W*s)
ddm <- ddply(ddm,.(variable),transform,csum=cumsum(value_J))

# strain names
ddm$strain <- 0; ddm$inoc <- 0; ddm$rep <- 0;
ddm[,c("strain_label","inoc_name")]<-colsplit(ddm$variable, "(?<=\\p{L})(?=[\\d+$])", c("char", "digit"))

w<- which(ddm$strain_label == "B")
ddm[w,"strain"] = "11006A"
ddm[w,"rep"] = "A"
w<- which(ddm$strain_label == "C")
ddm[w,"strain"] = "11006B"
ddm[w,"rep"] = "B"
w<- which(ddm$strain_label == "D")
ddm[w,"strain"] = "11016A"
ddm[w,"rep"] = "A"
w<- which(ddm$strain_label == "E")
ddm[w,"strain"] = "11016B"
ddm[w,"rep"] = "B"

# inoculum size
w<- which(ddm$inoc_name == "2")
ddm[w,"inoc"] = 7
w<- which(ddm$inoc_name == "3")
ddm[w,"inoc"] = 6
w<- which(ddm$inoc_name == "4")
ddm[w,"inoc"] = 5
w<- which(ddm$inoc_name == "5")
ddm[w,"inoc"] = 4
w<- which(ddm$inoc_name == "6")
ddm[w,"inoc"] = 3
w<- which(ddm$inoc_name == "7")
ddm[w,"inoc"] = 2
w<- which(ddm$inoc_name == "8")
ddm[w,"inoc"] = 1



# Check plots
ddm$inoc <- as.factor(ddm$inoc)

ggplot(ddm, aes(x=Time,y = value)) + geom_line(aes(col = inoc)) + facet_wrap(~strain)
ggsave("data_paper2/plots/CS_T0_curves.png")
ggplot(ddm, aes(x=Time,y = csum)) + geom_line(aes(col = inoc)) + facet_wrap(~strain)
ggsave("data_paper2/plots/CS_T0_csumcurves.png")

# OUTPUT
ddm <- subset(ddm, select = -c(strain_label,inoc_name)) # remove label
write.csv(ddm, "data_paper2/ddm_CS.csv")

####################################Organize OD data ####
dd <- read.table("data_paper2/OD_11006_11016_T0.txt",header = TRUE)
dd$Time <- dd$Time / 3600
ddm <-reshape2::melt(dd, id.vars=c("Time"))

ddm$strain <- 0; ddm$inoc <- 0
ddm[,c("strain_label","inoc_name")]<-colsplit(ddm$variable, "(?<=\\p{L})(?=[\\d+$])", c("char", "digit"))

# strain names
w<- which(ddm$strain_label == "A")
ddm[w,"strain"] = "11006A"
w<- which(ddm$strain_label == "B")
ddm[w,"strain"] = "11006B"
w<- which(ddm$strain_label == "D")
ddm[w,"strain"] = "11016A"
w<- which(ddm$strain_label == "E")
ddm[w,"strain"] = "11016B"
# inoculum size
w<- which(ddm$inoc_name == "1")
ddm[w,"inoc"] = 7
w<- which(ddm$inoc_name == "2")
ddm[w,"inoc"] = 6
w<- which(ddm$inoc_name == "3")
ddm[w,"inoc"] = 5
w<- which(ddm$inoc_name == "4")
ddm[w,"inoc"] = 4
w<- which(ddm$inoc_name == "5")
ddm[w,"inoc"] = 3
w<- which(ddm$inoc_name == "6")
ddm[w,"inoc"] = 2
w<- which(ddm$inoc_name == "7")
ddm[w,"inoc"] = 1

ddm$value <- ddm$value - 0.098 ##constant

# Check plots
ddm$inoc <- as.factor(ddm$inoc)

ggplot(ddm, aes(x=Time,y = value)) + geom_line(aes(col = inoc)) + facet_wrap(~strain)
ggsave("data_paper2/plots/OD_T0_curves.png")

#odplot <- ggplot(ddms, aes(x=Time,y = value)) + geom_line(aes(group = inoc))
#odplot
#sacplot <- grid.arrange(odplot, heatplot, nrow = 1, top="SAC042W")
#ggsave("output/sacplot.png", sacplot, width = 12, height = 5)

# OUTPUT
ddm <- subset(ddm, select = -c(strain_label,inoc_name)) # remove label
write.csv(ddm, "data_paper2/ddm_OD.csv")


# ##Model fitting CS heat flow data ####
# 
# #Grab CS data
# ddm <- as.data.table(read.csv("ddm_CS.csv")[,-1])
# 
# #Create matrices
# u <- as.character(unique(ddm$strain)) #unique strains
# q <- unique(ddm$inoc) #unique inocula
# 
# param_n <- matrix(0, length(u)*length(q)*1, 2); # Strain x inoc x number of experimental conditions x number of calculated variables
# param <- matrix(0, length(u)*length(q)*1, 5); #Strain x inoc x number of experimental conditions x number of calculated variables
# index <- 1
# 
# # Calculate parameters and fill matrices
# 
# for (jj in 1:length(u)){ #for each strain
#     for (ll in 1:length(q)){# for each inoculum 
#       print(c(jj, ll))
#       
#       #select data for unique combination of u, r and q
#       w <- intersect(which(ddm$inoc == q[ll]),which(ddm$strain == u[jj]))
#       #w <- intersect(which(ddm$inoc == 7), which(ddm$strain == "11006A"))
#       
#       if(length(w) > 0){
#         data1 <- ddm[w,]
#         
#         #Fit growth curve
#         gc_fit <- gcFitSpline(data1$Time, data1$csum)
#         #Calculate growth parameters
#         wmax <- which.max(data1$value_J)
#         time_max_heat_flow <- as.numeric(data1[wmax, "Time"])
#         value_max_heat_flow <- as.numeric(data1[wmax, "value"])
#         print(c(u[jj], time_max_heat_flow, value_max_heat_flow))
#         
#         #Plot data1, csum and fit
#         #data1
#         datam <- reshape2::melt(data1[,c("Time", "value_J", "csum", "strain")], id.vars = c("Time", "strain"))
#         gg <- ggplot(datam, aes(x=Time, y=value)) + geom_point() + facet_wrap(~strain, scales = "free")
#         gg
#         #model fit
#         gc_df <- as.data.frame(cbind(gc_fit$fit.time, gc_fit$fit.data))
#         colnames(gc_df) <- c("Time", "csum")
#         #gc_df$value_J <- c(gc_df[1, "csum"], diff(gc_df$csum))
#         gc_dfm <- reshape2::melt(gc_df, id.vars = "Time")
#         #add fit to plot
#         gg <- gg + geom_line(data = gc_dfm, aes(x=Time, y=value), col = "red")
#         gg
#         ##save
#         #ggsave(paste0("output/",u[jj],"_rep_", r[ii],"_inoc_", q[ll], "_model_fit_CS.pdf"))
#         
#         #Plot heat flow, reference line time peak
#         ggplot(data1, aes(x=Time, y=value)) + geom_point() +
#           geom_vline(xintercept = time_max_heat_flow, col = "black") + geom_text(aes(time_max_heat_flow, 0, label = round(time_max_heat_flow, digits = 2))) +
#           scale_x_continuous("Time (h)", limits = c(0, 25), minor_breaks = seq(0, 25, 5)) + 
#           scale_y_continuous(expression(paste("Heatflow (", mu, "W)")), limits = c(0,125), breaks = seq(0, 125, 25)) +
#           ggtitle(paste("Strain",u[jj]," inoc",q[ll]))
#         ggsave(paste0("",u[jj],"_inoc_", q[ll], "_heatflow_timemax_CS.png"))
#         
#         #Save parameters
#         s <- summary(gc_fit)
#         #u_d <- levels(droplevels(u))
#         #param_n[index,] <- c(unlist(strsplit(u_d[jj],"")), r[ii], q[ll])
#         param_n[index,] <- c(u[jj], q[ll])
#         #param[index,] <- c(time_max_heat_flow, value_max_heat_flow, s$mu.spline, s$lambda.spline,s$integral.spline)
#         param[index,] <- c(time_max_heat_flow, value_max_heat_flow, gc_fit$parameters$mu, gc_fit$parameters$lambda, gc_fit$parameters$integral)
#         index <- index + 1
#         
#       }
#     }
#   }
# 
# 
# param <- as.data.frame(param)
# param_n <- as.data.frame(param_n)
# colnames(param_n) <- c("strain","inoc")
# colnames(param) <- c("tmax", "ymax", "expgr","lag","auc")
# 
# param_orig_cs <- param
# param$strain <- param_n$strain
# param$inoc <- param_n$inoc
# 
# #w<-which(param[,1] == 0)
# #if(length(w) > 0){param <- param[-w,]} # remove 0 rows only if there are some
# 
# dim(param) 
# 
# #Output
# write.csv(param, "param_CS.csv")
# 
# ##Model fitting OD data ####
# 
# #Grab OD data
# ddm <- as.data.table(read.csv("ddm_OD.csv")[,-1])
# #ggplot(ddm, aes(x=Time,y = value)) + geom_line(aes(group = inoc)) + facet_wrap(~strain)
# 
# #Create matrices
# u <- as.character(unique(ddm$strain)) #unique strains
# q <- unique(ddm$inoc) #unique inocula
# 
# param_n <- matrix(0, length(u)*length(q)*1, 2); # Strain x rep x inoc x number of experimental conditions x number of calculated variables
# param <- matrix(0, length(u)*length(q)*1, 5); #Strain x rep x inoc x number of experimental conditions x number of calculated variables
# index <- 1
# 
# # Calculate parameters and fill matrices
# 
# for (jj in 1:length(u)){ #for each strain
#     for (ll in 1:length(q)){# for each inoculum 
#       print(c(jj, ll))
#       
#       #select data for unique combination of u, r and q
#       w <- intersect(which(ddm$inoc == q[ll]),which(ddm$strain == u[jj]))
#       #w <- intersect(which(ddm$inoc == 7), which(ddm$strain == "11006A"))
#       
#       if(length(w) > 0){
#         data1 <- ddm[w,]
#         
#         #Fit growth curve
#         gc_fit <- gcFitSpline(data1$Time, data1$value)
#         #plot.gcFitSpline(gc_fit)
#         #dev.copy(png, "output/_modelspline_OD.png")
#         
#         #Calculate growth parameters
#         y.spl     <- smooth.spline(data1$Time, data1$value)
#         dydt.spl   <- predict(y.spl, data1$Time, deriv = 1)
#         inde      <- which.max(dydt.spl$y)         
#         t.max      <- dydt.spl$x[inde]
#         dydt.max   <- max(dydt.spl$y)
#         y.max      <- y.spl$y[inde]
#         mu.spl     <- dydt.max;
#         b.spl      <- y.max-dydt.max*t.max
#         
#         #Plot data1, splines, tangent and ablines
#         gc_df <- as.data.frame(cbind(gc_fit$fit.time, gc_fit$fit.data))
#         colnames(gc_df) <- c("Time", "value")
#         
#         ggplot(data1, aes(x=Time, y=value)) + geom_point() +
#           geom_line(data = gc_df, aes(x=Time, y=value), col = "red") +
#           geom_abline(intercept = b.spl, slope = mu.spl, col = "red") +
#           geom_vline(xintercept = t.max, col = "black") + geom_text(aes(t.max, 0, label = round(t.max, digits = 2))) +
#           scale_x_continuous("Time (h)", limits = c(0, 25), minor_breaks = seq(0, 25, 5)) + 
#           scale_y_continuous("OD at 600nm", limits = c(0, 1), breaks = seq(0.25, 1, 0.25)) +
#           ggtitle(paste("Strain",u[jj]," inoc",q[ll]))
#         ggsave(paste0("",u[jj],"_inoc_", q[ll], "_model_fit_OD.png"))
#         
#         #Save parameters
#         s <- summary(gc_fit)
#         t.max <- as.numeric(t.max)
#         y.max <- as.numeric(y.max)
#         param_n[index,] <- c(u[jj], q[ll])
#         param[index,] <- c(t.max, y.max, gc_fit$parameters$mu, gc_fit$parameters$lambda, gc_fit$parameters$integral)
#         #param[index,] <- c(s$mu.spline, s$lambda.spline,s$integral.spline)
#         index <- index + 1
#         
#       }
#     }
#   }
# 
# 
# 
# param <- as.data.frame(param)
# param_n <- as.data.frame(param_n)
# colnames(param_n) <- c("strain","inoc")
# colnames(param) <- c("tmax", "ymax", "expgr","lag","auc")
# 
# param_orig_od <- param
# param$strain <- param_n$strain
# param$inoc <- param_n$inoc
# 
# dim(param) 
# 
# #Output
# write.csv(param, "param_OD.csv")
