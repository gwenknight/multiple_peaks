##### Double peaks exploration 


ddm <- read.csv("data_paper2/ddm_CS.csv")

###******** UNITS / DATA *************#######################################################################################################################################
#ddm$value_J = ddm$value

###******* MODEL FITTING *************#######################################################################################################################################
## Fit growth curves to each set of data to determine the underlying parameters 

## What are the strains?
u <- as.character(unique(ddm$strain,stringsasFactors = FALSE)) # strains
## How many replicates? 
r <- unique(ddm$rep) # replicates
# What are the inoculums? 
q <- unique(ddm$inoc)

# Where the parameters for each strain are stored
param <- matrix(0, length(u)*length(q)*3, 18); # number of strains x number of replicates x number of experimental conditions
index <- 1 # for counting 
max <- c() # for calibration

# Double curves
doubles <- matrix(0,0,14)

## Run thru each strain/rep/drying time/ inoculum: fit a growth curve to each
for(jj in 1:length(u)){ # for each strain
  r <- unique(ddm %>% filter(strain == u[jj]) %>% dplyr::select(rep))[,1] # which replicates for this strain (all have different names)
  for(ii in 1:length(r)){ # for each replicate
    #for(kk in c(1,2)){ #each of the experimental conditions
    for(ll in 1:length(q)){ #each of the inocula
      
      strain <- u[jj];
      replicate <- r[ii]
      #   condition <- drying_times[kk]
      inocl <- q[ll]
      
      wi <- intersect(which(ddm$strain == strain),which(ddm$rep == replicate)) # if fit to each replicate
      #  wj <- intersect(wi, which(ddm$drytime == condition))
      w <- intersect(wi, which(ddm$inoc == as.numeric(inocl)))
      
      if(length(w) > 0){ # if this replicate exists for this strain (i.e. there is data)
        data1 <- ddm[w,] # Grab data
        
        print(c(strain, replicate, inocl)) # output so can track how it is working
        
        p <- cut_extract(data1, "Time", "value_J", paste(strain, replicate, inocl,sep="_"),plot = 0, plot_where = "data_paper2/plots/", early_cut = 0) ### NEW function: runs on any timeseries: gives baseline and cut parameter analysis
        
        ## Required parameters
        param[index,] <- c(strain, replicate, inocl, p$param)
        index <- index + 1 # counting for storing matrix - next row each iteration
        
        ## Double
        if(p$param[10] == 1){
          if(!any(is.na(p$double_param))){
            curve <- "norm"
            doubles <- rbind(doubles, cbind(strain, replicate, inocl, curve, p$double_param))
          }
          if(!any(is.na(p$double_param_logn))){
            curve <- "log_norm"
            doubles <- rbind(doubles, cbind(strain, replicate, inocl, curve, p$double_param_logn))
          }
          if(!any(is.na(p$double_param_norm3))){
            curve <- "triple_norm"
            doubles <- rbind(doubles, cbind(strain, replicate, inocl, curve, p$double_param_norm3))
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
colnames(param) <- c("strain_name","rep","inocl",
                     "t_m_h_flow", "v_m_h_flow", "exp_gr","lag","auc",
                     "odd_peaks","odd_width","width_peak","odd_shoulder","odd_double","shoulder_point_t","shoulder_point_v",
                     "cut_exp", "timepeak", "valpeak")

w<-which(param$lag!=0); param <- param[w,] # remove 0 values
dim(param)
write.csv(param,"data_paper2/param_orig.csv")

### how many had double peaks? 
param %>% filter(odd_double > 0) 
dim(doubles)
ggplot(doubles %>% dplyr::select(strain:curve3) %>% dplyr::select(-value) %>% group_by(strain, replicate, inocl, time, curve) %>% pivot_longer(cols = "fit":"curve3"), 
       aes(x=time, y = value, group = name)) +  geom_point(data = doubles %>% dplyr::select(strain:value) %>% group_by(strain, replicate, inocl, time), aes(x=time, y = value), alpha = 0.3) +
  geom_line(aes(col = name)) + facet_grid(inocl ~ curve + strain + replicate , scales = "free") 

ggsave("data_paper2/plots/curve_fit_3.pdf")

