# EXPLORE BEHAVIOUR FOR ONE STRAIN 
# DATASET
strain <- "11223" # 11161
replicate <-11.2  #10.2
condition <- 0
inocl <- 5


#11050             6.2t168-inoc4/3
#11283             2.1t168-inoc5
#11179             10.3t0-inoc5/4 s


odd_shoulder <- 0
odd_width <- 0
odd_peak <- 0

####### DATA
data <- ddm
wi <- intersect(which(data$strain == strain),which(data$rep == replicate)) # if fit to each replicate
w <- intersect(wi, which(data$inoc == as.numeric(inocl)))

data_orig <- data[w,]  # just get the data for this experiment (time, value, )
data1 <- subset(data_orig, drytime == condition) 

gc_fit <- gcFitSpline(data1$Time, data1$value_J)
# parameters from this fit
s <- summary(gc_fit)

plot(gc_fit$raw.time, gc_fit$raw.data) 
lines(gc_fit$fit.time, gc_fit$fit.data)


##### WIDTH
wmax <- which.max( data1$value_J[-c(1:6)]) + 6 # take off funny initial behaviour < 2hr (take off in which.max and then add on 6 for position as taken 6 off)
time_max_heat_flow <- as.numeric(data1[wmax,"Time"])
value_max_heat_flow <- as.numeric(data1[wmax,"value_J"])

## ODD 
# (1) Is the peak broad? 
interval_peak <- time_max_heat_flow + c(-2.5, 2.5) # time interval around peak time
interval_value <- as.numeric(unlist(data1[c( which(round(data1$Time,4) == round(interval_peak[1],4)), which(round(data1$Time,4) == round(interval_peak[2],4))),"value_J"]))

max_level <- 0;
max_level <- max(max_level,round(100*interval_value / value_max_heat_flow,0)) # asymmetric peak then take highest

if(max_level >= thresh_wide){odd_width <- 1}


plot(data1$Time,data1$value_J, type = 'l', xlab = "Time", ylab = "Value_J")
points(time_max_heat_flow, value_max_heat_flow, col = 'red', pch = 19)
points(interval_peak, interval_value, pch = 19)
lines(interval_peak, c(0.9*value_max_heat_flow,0.9*value_max_heat_flow), col = "red")
lines(c(interval_peak[1],interval_peak[1]), c(-0.020,0.1), col = "red", lty = "dashed")
lines(c(interval_peak[2],interval_peak[2]), c(-0.020,0.1), col = "red", lty = "dashed")



##### PEAKS
# (2) Are there multiple peaks? 
# GIVES WHERE ALL PEAKS ARE (greater or equal to (m) 5 points around them)
peaks_index = find_peaks(data1$value_J, m = 5)
if(length(peaks_index) > 1){
  peaks_index_orig = find_peaks(data1$value_J, m = 5)
  # REMOVE - early ones
  we<-which(data1[peaks_index,"Time"]>3) 
  # REMOVE - late ones
  w<-intersect(we,which(data1[peaks_index,"Time"]< 0.90*max(data1$Time))) # remove > 95% of time
  wl<-intersect(w,which(data1[peaks_index,"Time"] > 0.6*max(data1$Time))) # which in the odd 70-95% of the time range
  if(length(wl)>0){ # if a late peak check it is big 
    if(data1[peaks_index[wl],"value_J"] < 0.45*max(data1[peaks_index,"value_J"])){ # if not bigger than 40% of peak
      w <-setdiff(w,wl)    # then remove
    }
  }
  
  peaks_index <- peaks_index[w] # Keep the ones not at the beginning / end
  
  # Sort by height - want to compare to and keep the tallest (first now in peaks_index)
  o <- order(data1[peaks_index,"value_J"], decreasing = "TRUE")
  peaks_index <- peaks_index[o]
}

# When are the peaks? 
time_peaks <- as.numeric(unlist(data1[peaks_index,"Time"]))
time_peaks_diff <- time_peaks - time_peaks[1] # how far apart are they?
# If multiple far apart then issue: double peaks
keep_time_far_apart <- which(abs(time_peaks_diff) > 5) 

if(length(keep_time_far_apart) >= 1){odd_peak <- 1} # if multiple peaks
if(any(data1[peaks_index,"value_J"] < 0.001)){odd_peak <- 1} # or if peak low

# If close and same height (90% of tallest) then odd (peak decline plateau decline OK)
close_peaks <- which(abs(time_peaks_diff) <= 5)
close_peaks_i <- 1
if(length(close_peaks)>1){
  for(i in 2:length(close_peaks)){
    ifelse(data1[peaks_index[close_peaks[i]],"value_J"]/data1[peaks_index[close_peaks[1]],"value_J"]> 0.9,
           close_peaks_i <- c(close_peaks_i,i),"")
  }}

if(length(close_peaks_i) > 1){odd_peak <- 1} # if multiple close time and height peaks


plot(data1$Time,data1$value_J, type = 'l',xlab = "Time", ylab = "Value_J")
points(unlist(data1[peaks_index[c(1)],"Time"]), unlist(data1[peaks_index[c(1)],"value_J"]), col = 'black', pch = 19)
points(unlist(data1[peaks_index[-c(1)],"Time"]), unlist(data1[peaks_index[-c(1)],"value_J"]), col = 'red', pch = 19)
lines(c(time_peaks[1]-5,time_peaks[1]+5),c(data1[peaks_index[1],"value_J"],data1[peaks_index[1],"value_J"]))

# (2) Are there multiple peaks? 
# GIVES WHERE ALL PEAKS ARE (greater or equal to (m) 5 points around them)
peaks_index = find_peaks(data1$value_J, m = 5)
# REMOVE - early ones
w<-which(data1[peaks_index,"Time"]>3) 
# REMOVE - late ones
w<-intersect(w,which(data1[peaks_index,"Time"]< 0.95*max(data1$Time))) 
peaks_index <- peaks_index[w] # Keep the ones not at the beginning / end
# Sort by height - want to compare to and keep the tallest (first now in peaks_index)
o <- order(data1[peaks_index,"value_J"], decreasing = "TRUE")
peaks_index <- peaks_index[o]

# When are the peaks? 
time_peaks <- as.numeric(unlist(data1[peaks_index,"Time"]))
time_peaks_diff <- time_peaks - time_peaks[1] # how far apart are they?
# If multiple far apart then issue: double peaks
keep_time_far_apart <- which(abs(time_peaks_diff) > 5) 

if(length(keep_time_far_apart) >= 1){odd_peak <- 1} # if multiple peaks
if(any(data1[peaks_index,"value_J"] < 0.001)){odd_peak <- 1} # or if peak low

# If close and same height (90% of tallest) then odd (peak decline plateau decline OK)
close_peaks <- which(abs(time_peaks_diff) <= 5)
close_peaks_i <- 1
if(length(close_peaks)>1){
  for(i in 2:length(close_peaks)){
    ifelse(data1[peaks_index[close_peaks[i]],"value_J"]/data1[peaks_index[close_peaks[1]],"value_J"]> 0.9,
           close_peaks_i <- c(close_peaks_i,i),"")
  }}

if(length(close_peaks_i) > 1){odd_peak <- 1} # if multiple close time and height peaks


plot(data1$Time,data1$value_J, type = 'l',xlab = "Time", ylab = "Value_J")
points(unlist(data1[peaks_index[c(1)],"Time"]), unlist(data1[peaks_index[c(1)],"value_J"]), col = 'black', pch = 19)
points(unlist(data1[peaks_index[c(2,3)],"Time"]), unlist(data1[peaks_index[c(2,3)],"value_J"]), col = 'red', pch = 19)
lines(c(time_peaks[1]-5,time_peaks[1]+5),c(0.9*data1[peaks_index[1],"value_J"],0.9*data1[peaks_index[1],"value_J"]), col="red", lty = "dashed")
lines(c(time_peaks[3],time_peaks[3]),c(0.9*data1[peaks_index[1],"value_J"],data1[peaks_index[3],"value_J"]), col="red", lty = "dashed")
lines(c(time_peaks[2],time_peaks[2]),c(0.9*data1[peaks_index[1],"value_J"],data1[peaks_index[2],"value_J"]), col="red", lty = "dashed")
lines(c(time_peaks[1],time_peaks[1]),c(0.9*data1[peaks_index[1],"value_J"],data1[peaks_index[1],"value_J"]))




#OLD
# ### SHOULDER
# time_endline <- max(min(data1$Time)+0.5, time_peaks[1] - 10) # 10 hrs from first peak or at least 30min in to recording data
# where.end <- which.min(abs(data1$Time - time_endline)) # What time exactly is this? 
# value_endline <- data1[where.end, "value_J"] # Find value at this time point
# ## Start of line is where? (top point)
# time_startline <- time_peaks[1] # highest peak
# value_startline <- as.numeric(data1[peaks_index[1],"value_J"])
# 
# # Draw straight line, assuming peak time = time 0. 
# grad = as.numeric((value_endline - value_startline)/(time_endline - time_startline)) # Gradient of line
# times_line <- as.numeric(unlist(data1[peaks_index[1]:where.end,"Time"])) # The times for the line (x values)
# 
# # Predicted straight line
# pred_points_fit <- grad*(times_line - time_startline) + value_startline # remove times_startline as taking this to be time zero (so can use value_startline as y intercept)
# 
# ### Visualise where peaks are if needed
# plot(data1$Time,data1$value_J, type = 'l', xlab = "Time", ylab = "Value_J")
# points(unlist(data1[peaks_index,"Time"]), unlist(data1[peaks_index,"value_J"]), col = 'black', pch = 19)
# ## and can plot line against this if needed too 
# lines(times_line, pred_points_fit, col= "blue")
# 
# #### How far is the predicted straight line from the data? 
# dist <- pred_points_fit - data1[peaks_index[1]:where.end,"value_J"]
# 
# ## Move the line to fit as close to the data as possible 
# # Want there to be some wobble on the line - if the predicted line never goes below actual data line then 
# # the distance will always be positive. Want there to be some later negative values or at least a close fit initially
# # Hence the removal of the first few terms ("-c(1:10)")- near the peak is always some wobble
# 
# # While there are not enough negative points: i.e. not enough overlap between the straight line and the curve
# # And while the time of the endline is not too close to peak
# # move the end of the line (bottom point) nearer to the peak time (can't go closer than 3 hrs)
# while(all(round(dist[-c(1:10)],10) >= 0) & time_endline < (time_peaks[1]-3)){ 
#   
#   time_endline <- as.numeric(time_endline + 1) # 10 hrs from first peak # SUBTRACT ONE EACH TIME
#   where.end <- which.min(abs(data1$Time - time_endline))
#   value_endline <- as.numeric(data1[where.end, "value_J"]) # Find value nearest this time point
#   
#   time_startline <- time_peaks[1]
#   value_startline <- as.numeric(data1[peaks_index[1],"value_J"])
#   
#   # straight line 
#   grad = (value_endline - value_startline)/(time_endline - time_startline)
#   times_line <- as.numeric(unlist(data1[peaks_index[1]:where.end,"Time"]))
#   
#   pred_points_fit <- grad*(times_line - time_startline) + value_startline # remove times_startline as taking this to be time zero (so can use value_startline as y intercept)
#   
#   dist <- pred_points_fit - as.numeric(unlist(data1[peaks_index[1]:where.end,"value_J"]))
#   
# }
# 
# ## No look at the distance from the straight line to the data: 
# ## if there are multiple peaks then have the bump then the data crosses the straight line multiple times
# ## = a shoulder 
# fpsq <- find_peaks(as.numeric(unlist(dist)), m = 5) # Are there more than 1 peaks? 
# if(length(fpsq)>1){odd_shoulder <- 1} # if yes = shoulder, otherwise just a smooth curve 
# 
# 
# odd_shoulder
# 





### SHOULDER 222222
time_endline <- -40 #max(min(data1$Time)+0.5, time_peaks[1] - 10) # 10 hrs from first peak or at least 30min in to recording data
value_endline <- 0
#where.end <- which.min(abs(data1$Time - time_endline)) # What time exactly is this? 
#value_endline <- data1[where.end, "value_J"] # Find value at this time point
## Start of line is where? (top point)
time_startline <- time_peaks[1] # highest peak
value_startline <- as.numeric(data1[peaks_index[1],"value_J"])

# Draw straight line, assuming peak time = time 0. 
grad = as.numeric((value_endline - value_startline)/(time_endline - time_startline)) # Gradient of line
times_line <- as.numeric(unlist(data1[peaks_index[1]:1,"Time"])) # The times for the line (x values)

# Predicted straight line
pred_points_fit <- grad*(times_line - time_startline) + value_startline # remove times_startline as taking this to be time zero (so can use value_startline as y intercept)

### Visualise where peaks are if needed
plot(data1$Time,data1$value_J, type = 'l', xlab = "Time", ylab = "Value_J", xlim = c(-40,25))
points(unlist(data1[peaks_index,"Time"]), unlist(data1[peaks_index,"value_J"]), col = 'black', pch = 19)
## and can plot line against this if needed too 
lines(times_line, pred_points_fit, col= "blue")

### Run through lines. 
odd_shoulder<- 0
shoulder_point <- time_peaks[1] # shoulder before this: this is the max peak 
### Visualise where peaks are if needed

plot(data1$Time,data1$value_J, type = 'l', xlab = "Time", ylab = "Value_J", xlim = c(0,25))
points(data1$Time,data1$value_J)
points(unlist(data1[peaks_index,"Time"]), unlist(data1[peaks_index,"value_J"]), col = 'black', pch = 19)

for(i in -40:5){
  print(i)
  time_endline <- i #max(min(data1$Time)+0.5, time_peaks[1] - 10) # 10 hrs from first peak or at least 30min in to recording data
  value_endline <- 0
  #where.end <- which.min(abs(data1$Time - time_endline)) # What time exactly is this? 
  #value_endline <- data1[where.end, "value_J"] # Find value at this time point
  ## Start of line is where? (top point)
  time_startline <- time_peaks[1] # highest peak
  value_startline <- as.numeric(data1[peaks_index[1],"value_J"])
  
  # Draw straight line, assuming peak time = time 0. 
  grad = as.numeric((value_endline - value_startline)/(time_endline - time_startline)) # Gradient of line
  # Only want exponential growth line
  exp_start <- which.min(abs(data1$Time-s$lambda.spline+1))
  times_line <- as.numeric(unlist(data1[peaks_index[1]:exp_start,"Time"])) # The times for the line (x values)
  
  # Predicted straight line
  pred_points_fit <- grad*(times_line - time_startline) + value_startline # remove times_startline as taking this to be time zero (so can use value_startline as y intercept)
  
  
  ## and can plot line against this if needed too 
  lines(times_line, pred_points_fit, col= "blue")
  
  #### How far is the predicted straight line from the data? 
  #if(length(pred_points_fit) > 28){leng = 28}else{leng = length(pred_points_fit)}
  dist <- as.numeric(unlist(c(pred_points_fit - data1[peaks_index[1]:exp_start,"value_J"])))
  
  # if(sum(dist>0) > 0 && sum(dist<0)>0){ # need dist to cross multiple times
  #   fpsq <- find_peaks(as.numeric(unlist(dist)), m = 5) # Are there more than 1 peaks? 
  #   
  #   if(length(fpsq)>1 ){ # if there are more than one peak
  #     if(max(-dist[fpsq]) > 0.0001) # if the curve moves a long way above the straight line 
  #       if(abs(sum(sign(dist[fpsq]))) != length(fpsq)){ # Want each peak to be on a different side of the line
  #         if((shoulder_point-1) > max(times_line[fpsq])){ # if the crossing is more than an hour from peak
  #           shoulder_point = min(shoulder_point-1,max(times_line[fpsq]))
  #           print(c(i,fpsq,times_line[fpsq])); odd_shoulder <- 1
  #         }}}} # if yes = shoulder, otherwise just a smooth curve 
  
  ## Crossing points
  #if(sum(diff(sign(dist)) != 0) > 2){ # if cross more than twice then shoulder
  if(length(which(abs(diff(sign(dist)))==2)) > 2){
    if(max(-dist) > 0.001){ # far enough away to be a proper shoulder
      if(max(dist) > 0.001){ # far enough away on both sides
        #if(times_line[which.max(-dist)] < (time_peaks[1] - 2)){ # and not just another peak
        shoulder_values = c(max(dist[which(dist>0)]), max(dist[which(dist<0)]))
        ws <- c() ####FINISH HERE
        for(i in 1:length(shoulder_values)){
          ws <- which(round(data1$Time,5) == round(shoulder_values[1],5)); 
          
          fpsq <- find_peaks(as.numeric(unlist(dist)), m = 3) # Are there more than 1 peaks?
          shoulder_point1 = min(time_peaks[1]-1,max(times_line[fpsq]))
          ws <- which(round(data1$Time,5) == round(shoulder_point1,5)); 
          shoulder_point_v1 <- data1[ws,"value_J"]
          if((shoulder_point_v1 / data1[peaks_index[1],"value_J"]) > 0.5){ # if the shoulder is greater than half way up
            if((shoulder_point_v1 / data1[peaks_index[1],"value_J"]) < 0.96){ # if the shoulder is not too near peak
              shoulder_point = shoulder_point1 
              shoulder_point_v = shoulder_point_v1
              print(shoulder_point_v / data1[peaks_index,"value_J"])
              odd_shoulder <- 1; print(c("new",i)) # Then odd shoulder     
            }
          }
        }
        #}
      }
    }
    
    
  }
  
}

ws <- which(data1$Time == shoulder_point)

# Run above for the i that is the one that gives the shoulder then plot this... 
plot(data1$Time,data1$value_J, type = 'l', xlab = "Time", ylab = "Value_J")#, xlim = c(0,10))#,ylim = c(-0.03,0.1))
points(unlist(data1[peaks_index,"Time"]), unlist(data1[peaks_index,"value_J"]), col = 'black', pch = 19)
lines(times_line[1:28], dist[1:28],col = "grey")
points(times_line[fpsq], data1[c(peaks_index[1]:1)[fpsq],"value_J"], col = 'red', pch = 19)
points(times_line[which.max(-dist)], data1[c(peaks_index[1]:1)[which.max(-dist)],"value_J"], col = 'green', pch = 19)
points(shoulder_point, shoulder_point_v, col = 'pink', pch = 19)
## and can plot line against this if needed too 
lines(times_line, pred_points_fit, col= "blue")


dist[fpsq]

# #### differential
# plot((data1$Time),(data1$value_J),ylim = c(-0.01,0.1))
# points(data1[ws,"Time"],data1[ws,"value_J"], col = "red")
# points(diff(diff(data1$value_J)))
# lines(diff(diff(data1$value_J)))
# #lines(diff(data1$value_J), col = "red")
# lines(seq(0,100,1),rep(0,101))
# 
# w <- which(diff(data1$value_J) < 0.00000000000000000000000000000000000001) 
# w[which(w <= peaks_index)]
# peaks_index_orig

#### If no shoulder but multiple peaks, want to grab time of first peak
if(length(time_peaks_diff) > 1 & odd_shoulder == 0){ # if multiple peaks but no shoulder
  shoulder_point1 = min(time_peaks)
  shoulder_point_v1 = data1[which(data1$Time == shoulder_point),"value_J"]
  if((shoulder_point_v1 / data1[peaks_index[1],"value_J"]) > 0.5){ # shoulder needs to be high still! 
    shoulder_point = shoulder_point1; shoulder_point_v = shoulder_point_v1 # Update shoulder values to ones that are > half way
  }
}

#### Past peak shoulder
###### SHOULDER - PAST PEAK
# (4b) Is there a shoulder after the peak? 
# GIVES distance from line to curve post peak - if far from line then there is a "shoulder"
# First, draw straight line from peak
# endline = bottom of line near end time 
time_endline <- 40 #max(min(ts[,Time])+0.5, time_peaks[1] - 10) # 10 hrs from first peak or at least 30min in to recording data
value_endline <- 0
## Start of line is where? (top point)
time_startline <- time_peaks[1] # highest peak
value_startline <- as.numeric(data1[peaks_index[1],"value_J"])

# Draw straight line, assuming peak time = time 0. 
grad = as.numeric((value_endline - value_startline)/(time_endline - time_startline)) # Gradient of line
# Only want exponential growth line
exp_End <- which.max(abs(data1$Time))
time_step = median(diff(data1$Time)) ## should be constant but some have variation so take normal step size, but should be constant to work in the below
times_line <- as.numeric(unlist(data1[peaks_index[1]:exp_End,"Time"]))# The times for the line (x values)

# Predicted straight line
pred_points_fit <- grad*(times_line - time_startline) + value_startline # remove times_startline as taking this to be time zero (so can use value_startline as y intercept)

## Start of line is where? (top point)
time_startline <- time_peaks[1] # highest peak
value_startline <- as.numeric(data1[peaks_index[1],"value_J"])

# Draw straight line, assuming peak time = time 0. 
grad = as.numeric((value_endline - value_startline)/(time_endline - time_startline)) # Gradient of line

### Visualise where peaks are if needed
# plot(unlist(data1[,"Time"]),unlist(data1[,"value_J"]), type = 'l', xlab = "Time", ylab = "value", xlim = c(0,40))
# points(time_startline, value_startline, col = 'black', pch = 19)
# #and can plot line against this if needed too
# lines(times_line, pred_points_fit, col= "blue")

### Run through lines. 
odd_shoulder_past<- 0
shoulder_point_past <- time_startline # shoulder before this

### Visualise where peaks are if needed
plot(unlist(data1[,"Time"]),unlist(data1[,"value_J"]), type = 'l', xlab = "Time", ylab = "value", xlim = c(0,25))
points(time_startline, value_startline, col = 'black', pch = 19)

shoulder_point_past_keep <- max(times_line) # initialise as furthest point

for(i in round((time_peaks[1]+5),0):40){
  time_endline <- i #max(min(data1$Time)+0.5, time_peaks[1] - 10) # 10 hrs from first peak or at least 30min in to recording data
  value_endline <- 0
  #where.end <- which.min(abs(data1$Time - time_endline)) # What time exactly is this? 
  #value_endline <- data1[where.end, "value_J"] # Find value at this time point
  ## Start of line is where? (top point)
  time_startline <- time_peaks[1] # highest peak
  value_startline <- as.numeric(data1[peaks_index[1],"value_J"])
  
  # Draw straight line, assuming peak time = time 0. 
  grad = as.numeric((value_endline - value_startline)/(time_endline - time_startline)) # Gradient of line
  # Want to end of curve
  exp_End <- which.max(abs(data1$Time))
  times_line <- as.numeric(unlist(data1[peaks_index[1]:exp_End,"Time"])) # The times for the line (x values)
  
  # Predicted straight line
  pred_points_fit <- grad*(times_line - time_startline) + value_startline # remove times_startline as taking this to be time zero (so can use value_startline as y intercept)
  
  ## and can plot line against this if needed too 
  lines(times_line, pred_points_fit, col= "blue")
  
  #### How far is the predicted straight line from the data? 
  #if(length(pred_points_fit) > 28){leng = 28}else{leng = length(pred_points_fit)}
  integer_times <- data1 %>% filter(Time >= (min(times_line)-0.02), Time <= (max(times_line)+0.02))
  dist <- as.numeric(c(pred_points_fit - (unlist(integer_times[,"value_J"]))))
  
  print(c("new",i,max(-dist),max(dist),length(which(abs(diff(sign(dist)))==2))))
  # Crossing points
  if(length(which(abs(diff(sign(dist)))==2)) > 2){# if cross more than twice then shoulder
    if(max(-dist) > 0.001 & max(dist) > 0.01){# || max(abs(dist)>0.05)){
      #if(max(-dist) > 0.001){ # far enough away to be a proper shoulder
      #  if(max(dist) > 0.01){
      #if(times_line[which.max(-dist)] < (time_peaks[1] - 2)){ # and not just another peak
      fpsq <- find_peaks(as.numeric(unlist(abs(dist))), m = 3) # Are there more than 1 peaks?
      if(length(fpsq) != 0){ # if no peaks then not a clear shoulder
        ws <- which(times_line[fpsq] > time_peaks[1]+1.01) # want shoulder to be more than hour from peak 
        if(length(ws) != 0 ){
          shoulder_point1_past = min(times_line[fpsq][ws]) # And then the next closest one
          ws <- which.min(abs(unlist(data1[,"Time"]) - shoulder_point1_past)) ## instead of which(round(data1[,"Time"],5) == round(shoulder_point1_past,5)); 
          shoulder_point_v1_past <- data1[ws,"value_J"]
          height <- shoulder_point_v1_past / data1[peaks_index[1],"value_J"]
          if(height > 0.55 && height < 0.99){ # if the shoulder is greater than 51% up but not too close
            shoulder_point_past = shoulder_point1_past; shoulder_point_v_past = shoulder_point_v1_past# Update shoulder values to ones that are > half way
            odd_shoulder_past <- 1; print(c("new",i,max(-dist),max(dist))) # Then odd shoulder  
            lines(times_line, pred_points_fit, col= "red") # overlay with red line on plot
            shoulder_point_past_keep = min(shoulder_point_past_keep, shoulder_point_past) # keep the earliest
          }
        }
      }
    }
  }
  
}

odd_shoulder_past
if(odd_shoulder_past == 0){shoulder_point_past <- 0; shoulder_point_v_past <- 0 # no shoulder_point if no shoulder!
}else{ws <- which.min(abs(unlist(data1[,"Time"]) - shoulder_point_past_keep)); shoulder_point_past <- shoulder_point_past_keep; shoulder_point_v_past <- as.numeric(data1[ws,"value_J"])}

points(shoulder_point_past, shoulder_point_v_past, col = "red")

# ws <- which(data1$Time == shoulder_point)
# data1$Time[ws]
# 
# ###### 
# 
# try <- ddm_cut %>% filter(strain == "11280", rep == "2.3", drytime == 168, inoc == 5) %>% filter(Time > 5)
# g1 <- ggplot(try,aes(x=Time, y= value_J)) + geom_line() + geom_point()
# st <- c()
# for(i in 2:(-2 + dim(try)[1])){
#   try1 <- try[max(1,i-4):i,]
#   try2 <- try[(i+1):dim(try)[1],]
#   lm.1 <- lm(try1$value_J ~ try1$Time)
#   lm.2 <- lm(try2$value_J ~ try2$Time)
#   
#   g1 + geom_line(data = try1, aes(x=try1$Time, y = lm.1$coefficients[1] + lm.1$coefficients[2]*try1$Time), col ="red") + 
#     geom_line(data = try2, aes(x=try2$Time, y = lm.2$coefficients[1] + lm.2$coefficients[2]*try2$Time), col ="blue")
#   
#   st <- rbind(st, c(i,c(max(try1$Time),lm.1$coefficients[2],lm.2$coefficients[2])))
#   print(c(i,c(max(try1$Time),lm.1$coefficients[2],lm.2$coefficients[2])))
#   
#   
# }
# 
# st <- as.data.frame(st)
# colnames(st) <- c("i","maxtime","f_ang","s_ang")
# st$d <- 10000
# st$d[1:(dim(st)[1]-1)] = abs(diff(st$s_ang))
# timepeak = st[which.min(st$d),"maxtime"]
# valpeak = try[which(try$Time == timepeak),"value_J"]
# 
# # ggplot(st, aes(x=maxtime, y = f_ang)) + geom_line() + geom_point() + 
# #   geom_line(aes(x=maxtime, y = s_ang), col = "red")  + geom_point(aes(x=maxtime, y = s_ang), col = "red")  + 
# #   geom_vline(xintercept = 7, lty = "dashed") 
# 
# ggplot(try,aes(x=Time, y= value_J)) + geom_line() + 
#   geom_point(data = try[which(try$Time == timepeak),c("Time","value_J")])



# #fitting smoothing splines using smooth.spline(X,Y,df=...)
# fit1<-smooth.spline(try$Time,try$value_J,df=5) #16 degrees of freedom
# #Plotting both cubic and Smoothing Splines 
# plot(try$Time,try$value_J,col="grey")
# #adding cutpoints
# lines(fit1,col="red",lwd=2)




# ##### CUT
# strain <- "11050"
# replicate <- "6.2"
# condition <- "168"
# inocl <- "5"
# data <- ddm_cut
# timepeak <- 0
# valpeak <- 0
# 
# ## which rows of the data are this strain and replicate?
# wi <- intersect(which(data$strain == strain),which(data$rep == replicate)) # if fit to each replicate
# wj <- intersect(wi, which(data$drytime == condition))
# w <- intersect(wj, which(data$inoc == as.numeric(inocl)))
# 
# if(length(w) > 0){ # if this replicate exists for this strain (i.e. there is data)
#   data1 <- data[w,] # %>% filter(Time > 5)# just get the data for this experiment (strain, time, value, drying time)
#   
#   ## First time to peak 
#   peaks_index = find_peaks(data1$value_J, m = 3)
#   if(length(peaks_index) > 0){
#     w <- which(data1[peaks_index,"Time"] < 3)
#     if(length(w) > 0){peaks_index <- peaks_index[-w]}
#     if(length(peaks_index) > 0){ # if any later than 3 yrs
#       peaks_index = min(peaks_index)
#       timepeak = data1[peaks_index,"Time"]
#       valpeak = data1[peaks_index,"value_J"]
#     }
#     
#   }else{
#     valpeak = 0; timepeak = 0
#     if(dim(data1)[1]>4){
#       
#       
#       st <- c()
#       for(i in 2:(-2 + dim(data1)[1])){
#         data11 <- data1[max(1,i-4):i,]
#         data12 <- data1[(i+1):dim(data1)[1],]
#         lm.1 <- lm(data11$value_J ~ data11$Time)
#         lm.2 <- lm(data12$value_J ~ data12$Time)
#         
#         st <- rbind(st, c(i,c(max(data11$Time),lm.1$coefficients[2],lm.2$coefficients[2])))
#       }
#       
#       st <- as.data.frame(st)
#       colnames(st) <- c("i","maxtime","f_ang","s_ang")
#       st$d <- 10000
#       st$d[1:(dim(st)[1]-1)] = abs(diff(st$s_ang))
#       st_upper <- st%>% filter(maxtime > 0.50*max(st$maxtime))
#       w1 <- which.min(st_upper$d)
#       w2 <- which.min(st_upper$d[-w1])
#       timepeak = min(st_upper[w1,"maxtime"], st_upper[-w1,"maxtime"][w2])
#       valpeak = data1[which(data1$Time == timepeak),"value_J"]
#     }else{valpeak = 0; timepeak = 0}
#   }
#   
#   ## Check not too low: if cut point already good enough
#   if(valpeak < 0.65*max(data1$value_J)){valpeak <- data1[dim(data1)[1],"value_J"]; timepeak =data1[dim(data1)[1],"Time"] }
#   
#   ggplot(data1,aes(x=Time, y= value_J)) + geom_line() + 
#     geom_point(data = data1[which(round(data1$Time,2) == as.numeric(round(timepeak,2))),c("Time","value_J")],col="red") + 
#     ggtitle(paste0(strain,"_", replicate,"_", inocl,"_",condition))
# }

