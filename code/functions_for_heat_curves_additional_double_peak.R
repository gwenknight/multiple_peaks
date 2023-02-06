#### Functions to extract key parameters. 

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


#### Function updated to give initial fit, cut and then subsequent fit: builds on above code

### Input: times series data
### Output: time series data up to cut point: time of end of exponential growth
### Output: lag time / exponential growth / time of end of exponential growth 

#### Function

### Input: times series data
### Output: time series data up to cut point: time of end of exponential growth
### Output: lag time / exponential growth / time of end of exponential growth 

cut_extract_dp <- function(ts, Time, value, name4output, thresh_wide = 86, plot = 0, plot_where = "plots/", early_cut = 3.5){
  ## ts = timeseries
  ## Time = name of time column
  ## value = name of value column
  ## name4output = strain, replicate, condition, inocl = labels for output
  ## thresh_wide = 86: % what is a wide peak? 
  ## plot = 0:  don't plot. 1 to plot
  ## plot_where: location for files to output to  
  ## early_cut: for heat flow, remove the first 3.5hrs of data 
  
  ## is this strain replicate ODD? Set all ODD indicators to zero initially
  odd_peak <- 0 # 0 = one peak
  odd_width <- 0 # 0 = not a wide peak
  odd_shoulder <- 0 # 0 = no shoulder
  odd_shoulder_past <- 0 # 0 = no shoulder after peak
  max_level <- 0 # height of peak
  odd_double <- 0 # double curve fits this data
  
  ## (1) Fit spline
  ## This gives lag time and exponential growth rate cumulative
  gc_fit <- gcFitSpline(ts[,Time], ts[,value])
  # parameters from this fit
  s <- summary(gc_fit)
  
  ## (2) What is the maximum heat flow and when?
  wmax <- which.max( unlist(ts[,value])[-c(1:6)]) + 6 # take off funny initial behaviour < 2hr (take off in which.max and then add on 6 for position as taken 6 off)
  time_max_heat_flow <- as.numeric(ts[wmax,Time])
  value_max_heat_flow <- as.numeric(ts[wmax,value])
  
  ## ODD characteristic determination
  ## (3) Looking at peaks
  # Is the peak broad? 
  interval_peak <- time_max_heat_flow + c(-2.5, 2.5) # time interval around peak time
  interval_value <- as.numeric(unlist(ts[c( which(round(ts[,Time],4) == round(interval_peak[1],4)), which(round(ts[,Time],4) == round(interval_peak[2],4))),value]))
  
  max_level <- 0;
  max_level <- max(max_level,round(100*interval_value / value_max_heat_flow,0)) # asymmetric peak then take highest
  
  # ODD? 
  if(max_level >= thresh_wide){odd_width <- 1}
  
  # Are there multiple peaks? 
  # GIVES WHERE ALL PEAKS ARE (greater or equal to (m) 5 points around them)
  peaks_index = as.numeric(find_peaks(unlist(ts[,value]), m = 5))
  if(length(peaks_index) > 0){ # MAY BE NO PEAK
    if(length(peaks_index) > 1){
      # REMOVE - early ones
      we<-which(ts[peaks_index,Time]>4) 
      # REMOVE - late ones
      w<-intersect(we,which(ts[peaks_index,Time]< 0.90*max(ts[,Time]))) # remove > 95% of time
      wl<-intersect(w,which(ts[peaks_index,Time] > 0.6*max(ts[,Time]))) # which in the odd 70-95% of the time range
      if(length(wl)>0){ # if a late peak check it is big 
        for(gg in 1:length(wl)){
          if(ts[peaks_index[wl[gg]],value] < 0.45*max(ts[,value])){ # if not bigger than 40% of maximum value in timeseries
            w <-setdiff(w,wl[gg]) }}}   # then remove
      peaks_index <- peaks_index[w] # Keep the ones not at the beginning / end
      
      # Check height ok - only want places with more than 45% of maximum (remove those little noisy bumps)
      w <- which(ts[peaks_index,value] > 0.45*max(ts[,value]))
      peaks_index <- peaks_index[w]
      
      # Sort by height - want to compare to and keep the tallest (first now in peaks_index)
      o <- order(ts[peaks_index,value], decreasing = "TRUE")
      peaks_index <- peaks_index[o]
    } else{ # if only one, check really a high point: greater than 45% of max of data
      if(ts[peaks_index,value] < 0.45*max(ts[,value])){
        peaks_index <- NA} # if too small then remove
    }
    
    if(is.numeric(peaks_index)){  # If there remain peaks
      
      # When are the peaks? 
      time_peaks <- as.numeric(unlist(ts[peaks_index,Time]))
      time_peaks_diff <- time_peaks - time_peaks[1] # how far apart are they?
      # If multiple far apart then issue: double peaks
      keep_time_far_apart <- which(abs(time_peaks_diff) > 5) 
      
      if(length(keep_time_far_apart) >= 1){odd_peak <- 1} # if multiple peaks
      #if(any(ts[peaks_index,value] < 0.001)){odd_peak <- 1} # or if peak low
      
      # If close and same height (90% of tallest) then odd (peak decline plateau decline OK)
      close_peaks <- which(abs(time_peaks_diff) <= 5)
      close_peaks_i <- 1
      if(length(close_peaks)>1){
        for(i in 2:length(close_peaks)){
          ifelse(ts[peaks_index[close_peaks[i]],value]/ts[peaks_index[close_peaks[1]],value]> 0.9,
                 close_peaks_i <- c(close_peaks_i,i),"")
        }}
      
      #if(length(close_peaks_i) > 1){odd_peak <- 1} # if multiple close time and height peaks
      
      # Only keep those peaks that are far apart
      time_peaks <- time_peaks[c(1,keep_time_far_apart)]
      peaks_index <- peaks_index[c(1,keep_time_far_apart)]
    }
  }
  
  ###### SHOULDER
  # (4) Is there a shoulder? 
  # GIVES distance from line to curve post peak - if far from line then there is a "shoulder"
  # First, draw straight line from peak
  # endline = bottom of line near time zero
  time_endline <- -40 #max(min(ts[,Time])+0.5, time_peaks[1] - 10) # 10 hrs from first peak or at least 30min in to recording data
  value_endline <- 0
  ## Start of line is where? (top point)
  time_startline <- ifelse(length(peaks_index)==0,  max(ts[,Time]), ifelse(is.na(peaks_index),max(ts[,Time]),time_peaks[1])) # highest peak or last point
  value_startline <- ifelse(length(peaks_index)==0, as.numeric(ts[which(unlist(ts[,Time])==time_startline),value]), 
                            ifelse(is.na(peaks_index),as.numeric(ts[which(unlist(ts[,Time])==time_startline),value]),as.numeric(ts[peaks_index[1],value]))) 
  
  # Draw straight line, assuming peak time = time 0. 
  grad = as.numeric((value_endline - value_startline)/(time_endline - time_startline)) # Gradient of line
  # Only want exponential growth line
  exp_start <- as.numeric(unlist(ts[which.min(unlist(abs(ts[,Time]-(s$lambda.spline-1)))),Time]))
  time_step = median(diff(ts$Time)) ## should be constant but some have variation so take normal step size, but should be constant to work in the below
  times_line <- seq(time_startline,exp_start,by = -time_step) # The times for the line (x values)
  
  
  # Predicted straight line
  pred_points_fit <- grad*(times_line - time_startline) + value_startline # remove times_startline as taking this to be time zero (so can use value_startline as y intercept)
  
  ### Visualise where peaks are if needed
  # plot(unlist(ts[,Time]),unlist(ts[,value]), type = 'l', xlab = "Time", ylab = "value", xlim = c(-40,25))
  # points(time_startline, value_startline, col = 'black', pch = 19)
  # and can plot line against this if needed too
  # lines(times_line, pred_points_fit, col= "blue")
  
  ### Run through lines. 
  odd_shoulder<- 0
  shoulder_point <- time_startline # shoulder before this
  
  ### Visualise where peaks are if needed
  # plot(unlist(ts[,Time]),unlist(ts[,value]), type = 'l', xlab = "Time", ylab = "value", xlim = c(0,25))
  # points(time_startline, value_startline, col = 'black', pch = 19)
  
  for(i in -40:5){
    time_endline <- i #max(min(ts[,Time])+0.5, time_peaks[1] - 10) # 10 hrs from first peak or at least 30min in to recording data
    value_endline <- 0
    #where.end <- which.min(abs(ts[,Time] - time_endline)) # What time exactly is this? 
    #value_endline <- ts[where.end, value] # Find value at this time point
    
    # Draw straight line, assuming peak time = time 0. 
    grad = as.numeric((value_endline - value_startline)/(time_endline - time_startline)) # Gradient of line
    times_line <- seq(time_startline, exp_start,by = -time_step) # The times for the line (x values) #c(time_startline:exp_start) # The times for the line (x values)
    
    # Predicted straight line
    pred_points_fit <- grad*(times_line - time_startline) + value_startline # remove times_startline as taking this to be time zero (so can use value_startline as y intercept)
    
    ## and can plot line against this if needed too 
    #lines(times_line, pred_points_fit, col= "blue")
    
    #### How far is the predicted straight line from the data? 
    #if(length(pred_points_fit) > 28){leng = 28}else{leng = length(pred_points_fit)}
    integer_times <- ts %>% filter(Time >= (min(times_line)-0.02), Time <= (max(times_line)+0.02))
    dist <- as.numeric(c(pred_points_fit - rev(unlist(integer_times[,value]))))
    
    # Crossing points
    if(length(which(abs(diff(sign(dist)))==2)) > 2){# if cross more than twice then shoulder
      if(max(-dist) > 0.001){ # far enough away to be a proper shoulder
        if(max(dist) > 0.001){
          #if(times_line[which.max(-dist)] < (time_peaks[1] - 2)){ # and not just another peak
          fpsq <- find_peaks(as.numeric(unlist(dist)), m = 3) # Are there more than 1 peaks?
          if(length(fpsq) != 0){ # if no peaks then not a clear shoulder
            shoulder_point1 = min(time_peaks[1]-1,max(times_line[fpsq]))
            ws <- which(round(ts[,Time],5) == round(shoulder_point1,5)); shoulder_point_v1 <- ts[ws,value]
            height <- shoulder_point_v1 / ts[peaks_index[1],value]
            if(height > 0.55 && height < 0.96){ # if the shoulder is greater than half way up but not too close
              shoulder_point = shoulder_point1; shoulder_point_v = shoulder_point_v1 # Update shoulder values to ones that are > half way
              odd_shoulder <- 1;# print(c("new",i)) # Then odd shoulder  
            }
          }
        }
      }
    }
    
  }
  
  if(odd_shoulder == 0){shoulder_point <- 0; shoulder_point_v <- 0 # no shoulder_point if no shoulder!
  }else{ws <- which(round(ts[,Time],5) == round(shoulder_point,5)); shoulder_point_v <- as.numeric(ts[ws,value])}
  
  #### If no shoulder but multiple peaks, want to grab time of first peak
  if(length(peaks_index)!=0){
    if(length(time_peaks_diff) > 1 & odd_shoulder == 0){ # if multiple peaks but no shoulder
      shoulder_point1 = min(time_peaks)
      shoulder_point_v1 = ts[which(round(ts[,Time],5) == round(shoulder_point1,5)),value]
      if((shoulder_point_v1 / ts[peaks_index[1],value]) > 0.5){ # shoulder needs to be high still! 
        shoulder_point = shoulder_point1; shoulder_point_v = shoulder_point_v1 # Update shoulder values to ones that are > half way
      }
    }
  }
  
  
  ###### SHOULDER - PAST PEAK
  # (4b) Is there a shoulder after the peak? 
  # GIVES distance from line to curve post peak - if far from line then there is a "shoulder"
  # First, draw straight line from peak
  # endline = bottom of line near end time 
  time_endline <- 40 #max(min(ts[,Time])+0.5, time_peaks[1] - 10) # 10 hrs from first peak or at least 30min in to recording data
  value_endline <- 0
  ## Start of line is where? (top point)
  time_startline <- time_peaks[1] # highest peak
  value_startline <- as.numeric(ts[peaks_index[1],value])
  
  # Draw straight line, assuming peak time = time 0. 
  grad = as.numeric((value_endline - value_startline)/(time_endline - time_startline)) # Gradient of line
  # Only want exponential growth line
  exp_End <- which.max(abs(ts$Time))
  time_step = median(diff(ts$Time)) ## should be constant but some have variation so take normal step size, but should be constant to work in the below
  times_line <- as.numeric(unlist(ts[peaks_index[1]:exp_End,Time]))# The times for the line (x values)
  
  # Predicted straight line
  pred_points_fit <- grad*(times_line - time_startline) + value_startline # remove times_startline as taking this to be time zero (so can use value_startline as y intercept)
  
  ## Start of line is where? (top point)
  time_startline <- time_peaks[1] # highest peak
  value_startline <- as.numeric(ts[peaks_index[1],value])
  
  # Draw straight line, assuming peak time = time 0. 
  grad = as.numeric((value_endline - value_startline)/(time_endline - time_startline)) # Gradient of line
  
  ### Visualise where peaks are if needed
  # plot(unlist(ts[,Time]),unlist(ts[,value]), type = 'l', xlab = "Time", ylab = "value", xlim = c(0,40))
  # points(time_startline, value_startline, col = 'black', pch = 19)
  # #and can plot line against this if needed too
  # lines(times_line, pred_points_fit, col= "blue")
  
  ### Run through lines. 
  odd_shoulder_past<- 0
  shoulder_point_past <- time_startline # shoulder before this
  
  ### Visualise where peaks are if needed
  #plot(unlist(ts[,Time]),unlist(ts[,value]), type = 'l', xlab = "Time", ylab = "value", xlim = c(0,25))
  #points(time_startline, value_startline, col = 'black', pch = 19)
  
  shoulder_point_past_keep <- max(times_line) # initialise as furthest point
  
  for(i in round((time_peaks[1]+5),0):40){
    time_endline <- i #max(min(ts$Time)+0.5, time_peaks[1] - 10) # 10 hrs from first peak or at least 30min in to recording data
    value_endline <- 0
    #where.end <- which.min(abs(ts$Time - time_endline)) # What time exactly is this? 
    #value_endline <- ts[where.end, "value_J"] # Find value at this time point
    ## Start of line is where? (top point)
    time_startline <- time_peaks[1] # highest peak
    value_startline <- as.numeric(ts[peaks_index[1],value])
    
    # Draw straight line, assuming peak time = time 0. 
    grad = as.numeric((value_endline - value_startline)/(time_endline - time_startline)) # Gradient of line
    # Want to end of curve
    exp_End <- which.max(abs(ts$Time))
    times_line <- as.numeric(unlist(ts[peaks_index[1]:exp_End,Time])) # The times for the line (x values)
    
    # Predicted straight line
    pred_points_fit <- grad*(times_line - time_startline) + value_startline # remove times_startline as taking this to be time zero (so can use value_startline as y intercept)
    
    ## and can plot line against this if needed too 
    #lines(times_line, pred_points_fit, col= "blue")
    
    #### How far is the predicted straight line from the data? 
    #if(length(pred_points_fit) > 28){leng = 28}else{leng = length(pred_points_fit)}
    integer_times <- ts %>% filter(Time >= (min(times_line)-0.02), Time <= (max(times_line)+0.02))
    dist <- as.numeric(c(pred_points_fit - (unlist(integer_times[,value]))))
    
    
    # Crossing points
    if(length(which(abs(diff(sign(dist)))==2)) > 2){# if cross more than twice then shoulder
      if(max(-dist) > 0.001){ # far enough away to be a proper shoulder
        if(max(dist) > 0.01){
          #if(times_line[which.max(-dist)] < (time_peaks[1] - 2)){ # and not just another peak
          fpsq <- find_peaks(as.numeric(unlist(abs(dist))), m = 3) # Are there more than 1 peaks?
          if(length(fpsq) != 0){ # if no peaks then not a clear shoulder
            ws <- which(times_line[fpsq] > time_peaks[1]+1.01) # want shoulder to be more than hour from peak 
            if(length(ws) != 0 ){
              shoulder_point1_past = min(times_line[fpsq][ws]) # And then the next closest one
              ws <- which.min(abs(unlist(ts[,Time]) - shoulder_point1_past)) ## instead of which(round(ts[,Time],5) == round(shoulder_point1_past,5)); 
              shoulder_point_v1_past <- ts[ws,value]
              height <- shoulder_point_v1_past / ts[peaks_index[1],value]
              if(height > 0.55 && height < 0.96){ # if the shoulder is greater than 55% up but not too close
                shoulder_point_past = shoulder_point1_past; shoulder_point_v_past = shoulder_point_v1_past# Update shoulder values to ones that are > half way
                odd_shoulder_past <- 1; #print(c("new",i)) # Then odd shoulder  
                shoulder_point_past_keep = min(shoulder_point_past_keep, shoulder_point_past) # keep the earliest
              }
            }
          }
        }
      }
    }
    
  }
  
  if(odd_shoulder_past == 0){shoulder_point_past <- 0; shoulder_point_v_past <- 0 # no shoulder_point if no shoulder!
  }else{ws <- which.min(abs(unlist(ts[,Time]) - shoulder_point_past_keep)); shoulder_point_past <- shoulder_point_past_keep; shoulder_point_v_past <- as.numeric(ts[ws,value])}
  
  ##### (5) Double curves? Only use to count how many have this - fit not explored here
  # x <- ts[,Time]
  # y <- ts[,value]
  # 
  # # Guesses for the parameters in the curve
  # startl=list(a=unlist(ts[peaks_index,value])[1]*1.5, 
  #             b=diff(interval_peak)/3, 
  #             c=unlist(ts[peaks_index,Time])[1],
  #             d=unlist(ts[peaks_index,value])[1]*1.5, 
  #             e=diff(interval_peak), 
  #             f=unlist(ts[peaks_index,Time])[1])
  # # Set indicator/output parameters to zero
  # fit1 <- 0; fit2 <- 0; fit3 <- 0;
  # plot_p <- NA; plot_p2 <- NA; plot_p3 <- NA;
  # p <- 0
  # # Try to fit two normal curves using non-linear least squares methods 
  # try(fit1 <- nls(y~(a/b)*exp(-(x-c)^2/(2*b^2))+(d/e)*exp(-(x-f)^2/(2*e^2)),start = startl, algorithm="port"), silent = TRUE)
  # try(fit2 <- nls(y~peak1/(x*sig*sqrt(2*pi))*exp(-(log(x) - mu)^2/(2*sig^2)) + peak2/(x*sig2*sqrt(2*pi))*exp(-(log(x) - mu2)^2/(2*sig2^2)), 
  #                 start = list(peak1 = 800, peak2 = 900, sig = 0.1, mu = 1.7, sig2 = 0.1, mu2 = 1.9)), silent = TRUE)
  # try(fit3 <- nls(y~(a/(2*b))*exp(-(x-c)^2/(2*b^2))+(a/(2*b))*exp(-(x-c/2)^2/(2*b^2)) +(d/e)*exp(-(x-f)^2/(2*e^2)),start = startl, algorithm="port"), silent = TRUE)
  # 
  # # If manage to fit
  # if(length(fit1) > 1){
  #   # print that can 
  #   print(paste("Double curve fit Norm", name4output,sep = " "))
  #   odd_double <- 1
  #   # Save the predictions of the fit 
  #   pred_p <- predict(fit1)
  #   # Save the coefficients of the fit
  #   p <- coef(fit1)
  #   # Pull out the actual curves using the parameters in p
  #   g1p <- (p["a"]/p["b"])*exp(-(x-p["c"])^2/(2*p["b"]^2)) 
  #   g2p <- (p["d"]/p["e"])*exp(-(x-p["f"])^2/(2*p["e"]^2))
  #   if(all(g1p > 0) && all(g2p>0)){
  #     # Store the curves
  #     plot_p <- as.data.frame(cbind(x,y,pred_p, g1p, g2p,0))
  #     colnames(plot_p) <- c("time","value","fit","curve1","curve2","curve3")
  #     # Store the parameters
  #     plot_p$a <- p["a"]; plot_p$b <- p["b"]; plot_p$c <- p["c"]
  #     plot_p$d <- p["d"]; plot_p$e <- p["e"]; plot_p$f <- p["f"]
  #   }
  # }
  # 
  # 
  # # If manage to fit
  # if(length(fit2) > 1){
  #   # print that can 
  #   print(paste("Double curve fit Lognorm", name4output,sep = " "))
  #   odd_double <- 1
  #   # Save the predictions of the fit 
  #   pred_p <- predict(fit2)
  #   # Save the coefficients of the fit
  #   p <- coef(fit2)
  #   # Pull out the actual curves using the parameters in p
  #   sig = p["sig"]; mu = p["mu"]; peak1 <- p["peak1"]
  #   sig2 <- p["sig2"]; mu2 <- p["mu2"]; peak2 <- p["peak2"]
  #   g1p <- peak1/(x*sig*sqrt(2*pi))*exp(-(log(x) - mu)^2/(2*sig^2))
  #   g2p <- peak2/(x*sig2*sqrt(2*pi))*exp(-(log(x) - mu2)^2/(2*sig2^2))
  #   if(all(g1p > 0) && all(g2p>0)){
  #     # Store the curves
  #     plot_p2 <- as.data.frame(cbind(x,y,pred_p, g1p, g2p,0))
  #     colnames(plot_p2) <- c("time","value","fit","curve1","curve2","curve3")
  #     # Store the parameters
  #     plot_p2$a <- p["sig"]; plot_p2$b <- p["mu"]; plot_p2$c <- p["peak1"]
  #     plot_p2$d <- p["sig2"]; plot_p2$e <- p["mu2"]; plot_p2$f <- p["peak2"]
  #   }
  # }
  # 
  # # If manage to fit
  # if(length(fit3) > 1){
  #   # print that can 
  #   print(paste("Double curve fit Norm x3", name4output,sep = " "))
  #   odd_double <- 1
  #   # Save the predictions of the fit 
  #   pred_p <- predict(fit3)
  #   # Save the coefficients of the fit
  #   p <- coef(fit3)
  #   # Pull out the actual curves using the parameters in p
  #   g1p <- (p["a"]/(2*p["b"]))*exp(-(x-p["c"])^2/(2*p["b"]^2)) 
  #   g2p <- (p["a"]/(2*p["b"]))*exp(-(x-p["c"]/2)^2/(2*p["b"]^2)) 
  #   g3p <- (p["d"]/p["e"])*exp(-(x-p["f"])^2/(2*p["e"]^2))
  #   if(all(g1p > 0) && all(g2p>0) && all(g3p>0)){
  #     # Store the curves
  #     plot_p3 <- as.data.frame(cbind(x,y,pred_p, g1p, g2p, g3p))
  #     colnames(plot_p3) <- c("time","value","fit","curve1","curve2","curve3")
  #     # Store the parameters
  #     plot_p3$a <- p["a"]; plot_p3$b <- p["b"]; plot_p3$c <- p["c"]
  #     plot_p3$d <- p["d"]; plot_p3$e <- p["e"]; plot_p3$f <- p["f"]
  #   }
  # }
  # 
  # ## Plot ts, cumulative and fit
  # # if functions takes in a command to plot
  # if(plot == 1){
  #   
  #   ### Fit to value : could add cumulative plot in here but not done
  #   datam <- reshape2::melt(ts[,c(Time,value)], id.vars = Time)
  #   gg <- ggplot(datam, aes(x=Time,y=value)) + geom_point() 
  #   ## model fit:
  #   gc_df <- as.data.frame(cbind(gc_fit$fit.time, gc_fit$fit.data))
  #   colnames(gc_df) <- c(Time,value)
  #   #gc_df[,value]<- c(gc_df[1,"csum"],diff(gc_df$csum)) # CHANGE TO CUMULATIVE? 
  #   #gc_dfm <- reshape2::melt(gc_df, id.vars = Time)
  #   ## add fit to data plot
  #   gg <- gg + geom_line(data = gc_df, aes(x=Time,y=value), col= "red")
  #   
  #   ## save
  #   ggsave(paste0(plot_where, name4output, "_model_fit.pdf"))
  # } 
  # 
  ## Build vectors of required parameters to output
  param_o   <- c(time_max_heat_flow, value_max_heat_flow, 
                 s$mu.spline, s$lambda.spline,s$integral.spline, 
                 odd_peak, odd_width, max_level, odd_shoulder, odd_shoulder_past, odd_double, 
                 shoulder_point, as.numeric(shoulder_point_v), shoulder_point_past, as.numeric(shoulder_point_v_past))
  
  #print(param_o)  
  
  
  # (6) Now cut at shoulder or first peak
  if(shoulder_point > 0){
    cut_point_t <- shoulder_point; cut_point_v <- shoulder_point_v}else{
      cut_point_t <- time_max_heat_flow; cut_point_v <- value_max_heat_flow;
    }
  
  # NEW Time series generated: cut up to first peak or shoulder 
  #print(paste0("Early cut = ",early_cut))
  ts_cut <- ts %>% 
    filter(Time > early_cut) %>% # cut off first 3hrs 
    mutate(cutpart = ifelse(Time <= cut_point_t,1,0)) %>%
    filter(cutpart == 1) # Trim off extra parts
  
  ### Look at this cut data
  ## Currently the time to use is: peak growth 
  timepeak = cut_point_t
  valpeak = cut_point_v[1]
  
  ## (7) What if there is a peak in this?
  peaks_index = find_peaks(ts_cut[,value], m = 3)
  
  ## If there is a peak then reassign peak 
  if(length(peaks_index) > 0){
    w_early <- which(ts_cut[peaks_index,Time] < 5) # Remove early ones
    if(length(w_early) > 0){peaks_index <- peaks_index[-w_early]}
    w_high <- which(ts_cut[peaks_index,value] < 0.6*max(ts_cut[,value])) # Remove low ones
    if(length(w_high) > 0){peaks_index <- peaks_index[-w_high]}
    if(length(peaks_index) > 0){ # if any later than 3 
      peaks_index = min(peaks_index)
      timepeak = ts_cut[peaks_index,Time]
      valpeak = ts_cut[peaks_index,value]
    }
  }
  
  # Reassign up to peak 
  ts_cut = ts[which(ts[,Time] <= as.numeric(timepeak)),]
  
  ### (8) Check if the slope changes substantially in this period: may be a shoulder or a plateau near peak 
  if(dim(ts_cut)[1]>4){ # If enough data
    st <- c()
    
    #plot(ts_cut[,Time], ts_cut[,value], "l")
    #points_cut(ts_cut[,Time], ts_cut[,value])
    
    for(i in 2:(-1 + dim(ts_cut)[1])){ # fit a linear model to the data in segments_cut
      ts_cut1 <- ts_cut[max(1,i-6):i,]
      ts_cut2 <- ts_cut[(i):dim(ts_cut)[1],]
      lm.1 <- lm(unlist(ts_cut1[,value]) ~ unlist(ts_cut1[,Time]))
      lm.2 <- lm(unlist(ts_cut2[,value]) ~ unlist(ts_cut2[,Time]))
      
      #lines(ts_cut1[,Time], lm.1$coefficients[1] + lm.1$coefficients[2]*ts_cut1[,Time], col = "blue")
      #lines(ts_cut2[,Time], lm.2$coefficients[1] + lm.2$coefficients[2]*ts_cut2[,Time], col = "red")
      
      st <- rbind(st, c(i,c(max(ts_cut1[,value]), max(ts_cut1[,Time]),lm.1$coefficients[2],lm.2$coefficients[2])))
    }
    
    st <- as.data.frame(st)
    colnames(st) <- c("i","maxval","maxtime","f_ang","s_ang") # first angle, second angle
    st$d <- 0; st$da <- 1000
    st$d[1:(dim(st)[1]-1)] = diff(st$s_ang) # look at change in second angle: want to know when substantial change
    st$da[1:(dim(st)[1]-1)] = abs(diff(st$s_ang)) # look at change in second angle: want to know when substantial change
    
    #plot(st$maxtime, st$d) 
    #plot(ts_cut[,Time], ts_cut[,value])
    
    ## Check if shoulder
    st_upper <- st%>% filter(maxtime > 0.50*max(st$maxtime)) # but want to be past halfway
    w1 <- which.min(st_upper$da) # as absolute this detects_cut a plateau 
    w2 <- which.min(st_upper$da[-w1])
    timepeak_s = min(st_upper[w1,"maxtime"], st_upper[-w1,"maxtime"][w2]) # earliest 
    valpeak_s = as.numeric(ts_cut[which(ts_cut[,Time] == timepeak_s),value])
    
    ## Check not too low: if cut point already good enough. 
    # If OK then cut at shoulder value
    if(valpeak_s > 0.65*max(ts_cut[,value])){
      valpeak <- valpeak_s; timepeak = timepeak_s}else{
        
        ## Check if end peak should be moved forward at all (i.e. a peak: want end of exponential growth )
        st_upper <- st %>% filter(maxval > 0.80*max(st$maxval)) # want to be near the end 
        #if(dim(st_upper)[1] > 5){ # have at least 4 d values to look at
        w1 <- which.max(st_upper$d)
        w2 <- which.max(st_upper$d[-w1])
        if(min(w1,w2) != 1){ # if its_cut not just a slope down - if there is a plateau near the top? i.e. index now just the earliest point
          timepeak = st_upper[min(w1,w2),"maxtime"]
          valpeak = as.numeric(ts_cut[which(ts_cut[,Time] == timepeak),value])}
      }
  }
  
  if(plot == 1){
    g1 <- ggplot(ts_cut,aes(x=Time, y= value)) + geom_line() + 
      geom_point(data = ts_cut[which(round(ts_cut[,Time],2) == as.numeric(round(timepeak,2))),c(Time,value)],col="red") + 
      geom_point(data= ts_cut[1,], aes(x=shoulder_point, y = shoulder_point_v), col = "black") + 
      ggtitle(paste0(u[jj],"_",r[ii], "_",drying_times[kk],"_",q[ll]))
    dir.create(file.path(here(), paste0(plot_where,"/shoulder_curves")),showWarnings = FALSE)
    ggsave(paste0(plot_where,"shoulder_curves/cutpoint_highlighted_",name4output,".pdf")) 
  }
  
  # sames as ts_cut
  tstopeak = ts[which(ts[,Time] <= as.numeric(timepeak)),c(Time,value)]
  
  ## (9) Growth Curve 
  ## Ths gives lag time and exponential growth rate 
  gc_fit <- gcFitSpline(tstopeak[,Time], tstopeak[,value])
  # parameters from this fit
  s_peak <- summary(gc_fit)
  
  ## (10) Record all minor peaks 
  peaks_index_minor = as.numeric(find_peaks(unlist(ts[,value]), m = 5))
  if(length(peaks_index_minor) > 0){ # MAY BE NO PEAK
    if(length(peaks_index_minor) > 1){
      # REMOVE - early ones
      we<-which(ts[peaks_index_minor,Time]>4) 
      # REMOVE - late ones
      w<-intersect(we,which(ts[peaks_index_minor,Time]< 0.90*max(ts[,Time]))) # remove > 90% of time
      wl<-intersect(w,which(ts[peaks_index_minor,Time] > 0.6*max(ts[,Time]))) # which in the odd 70-95% of the time range
      if(length(wl)>0){ # if a late peak check it is big 
        for(gg in 1:length(wl)){
          if(ts[peaks_index_minor[wl[gg]],value] < 0.45*max(ts[,value])){ # if not bigger than 40% of maximum value in timeseries
            w <-setdiff(w,wl[gg]) }}}   # then remove
      peaks_index_minor <- peaks_index_minor[w] # Keep the ones not at the beginning / end
      
      # Check height ok - only want places with more than 45% of maximum (remove those little noisy bumps)
      w <- which(ts[peaks_index_minor,value] > 0.45*max(ts[,value]))
      peaks_index_minor <- peaks_index_minor[w]
      
      # # Sort by height - want to compare to and keep the tallest (first now in peaks_index)
      # o <- order(ts[peaks_index_minor,value], decreasing = "TRUE")
      # peaks_index <- peaks_index[o]
    } else{ # if only one, check really a high point: greater than 45% of max of data
      if(ts[peaks_index_minor,value] < 0.45*max(ts[,value])){
        peaks_index_minor <- NA} # if too small then remove
    }
    
    if(is.numeric(peaks_index_minor)){  # If there remain peaks
      
      # When are the peaks? 
      time_peaks <- as.numeric(unlist(ts[peaks_index_minor,Time]))
      time_peaks_diff <- time_peaks - time_peaks[1] # how far apart are they?
      # If multiple far apart then issue: double peaks
      keep_time_far_apart <- which(abs(time_peaks_diff) > 1) 
      
      #if(length(keep_time_far_apart) >= 1){odd_peak <- 1} # if multiple peaks
      #if(any(ts[peaks_index,value] < 0.001)){odd_peak <- 1} # or if peak low
      
      # If close and same height (90% of tallest) then odd (peak decline plateau decline OK)
      close_peaks <- which(abs(time_peaks_diff) <= 1)
      close_peaks_i <- 1
      if(length(close_peaks)>1){
        for(i in 2:length(close_peaks)){# first always 0 as same peak
          ifelse(ts[peaks_index[close_peaks[i]],value]/ts[peaks_index[close_peaks[1]],value]> 0.9,
                 close_peaks_i <- c(close_peaks_i,i),"")
        }}
      
      # if(length(close_peaks_i) >= 1){odd_peak <- 1} # if multiple close time and height peaks
      
      # Only keep those peaks that are far apart
      time_peaks <- time_peaks[c(1,keep_time_far_apart)]
      peaks_index_minor <- peaks_index_minor[c(1,keep_time_far_apart)]
      val_peaks <- unlist(ts[peaks_index_minor, value])
    }
  }
  
  ## Minor peaks 
  minor_peaks_timing <- c(time_peaks, matrix(0,1,10-length(time_peaks)))
  minor_peaks_height <- c(val_peaks, matrix(0,1,10-length(val_peaks)))
  gaps <- (time_peaks - time_peaks[1])[-1]
  minor_peaks_gaps <- c(gaps, matrix(0,1,9-length(gaps)))
  minor_peaks <- c(minor_peaks_timing, minor_peaks_height, minor_peaks_gaps)
  
  ## (11) AUC
  # From lag time to 10 hrs after
  lag_time <- as.numeric(which.min(abs(unlist(ts[,Time]) - s$lambda.spline)))
  ten_after = ifelse(max(ts[,Time])<100,10,10*60*60) # units in current data either hours or seconds (latter is 10 hrs of seconds)
  ten_after_lag = as.numeric(which.min(abs(unlist(ts[,Time]) - (s$lambda.spline + ten_after))))
  gc_fit_auc <- gcFitSpline(ts[c(lag_time:ten_after_lag),Time], ts[c(lag_time:ten_after_lag),value])
  # parameters from this fit
  s_auc <- summary(gc_fit_auc)
  
  ## Build vectors of required parameters to output
  param_o[5] <- s_auc$integral.spline # replace original auc with just that up to 10hrs from lag
  param_o_new <- c(param_o, as.numeric(s_peak$mu.spline), as.numeric(timepeak), as.numeric(valpeak), minor_peaks)
  
  return(list(param = param_o_new))
  
}

#### Function to determine clusters

### Input: parameter information and timeseries
### Output: clustering

cluster <- function(ts, parameters, name = "clusters", plot_where = "plots/"){
  
  # Add cluster column
  ts$cluster <- ""; parameters$cluster <- ""
  
  # Loop thru each condition 
  drytimes <- unique(ts$drytime)
  inocs <- unique(ts$inoc)
  
  # to bind to
  ts_end <- c()
  parameters_end <- c()
  
  for(j in 1:length(drytimes)){
    for(k in 1:length(inocs)){
      
      print(c(j,k))
      # subset data
      sub_ts <- ts %>% filter(drytime == drytimes[j], inoc == inocs[k])
      sub_parm <- parameters %>% filter(drytime == drytimes[j], inoc == inocs[k])
      
      ####**** (Cluster 1) clear double peaks
      ####* Are there any? 
      dp <- sub_parm %>% filter(odd_peaks == 1) 
      if(dim(dp)[1]>0){
        double_peak_curves_analysis <- as.data.frame(dp %>% group_by(strain) %>% dplyr::mutate(n_odd_peaks = n()) )
        
        nd <- left_join(as.data.frame(sub_ts), double_peak_curves_analysis[,c("strain","rep","odd_peaks","n_odd_peaks")], by = c("strain", "rep"))
        g <- ggplot(nd, aes(x=Time, y = value_J, group = rep)) + geom_line(aes(col = factor(n_odd_peaks))) + facet_wrap(~strain) + 
          scale_color_manual("Number of\nreps with\nodd peaks", breaks = c(3,2,1,"NA"), labels = c(3,2,1,0), values = c("red","orange","blue","grey"))
        ggsave(paste0(plot_where,drytimes[j],"_",inocs[k],"_",name,"_double_peaks.pdf"))
        
        # which are in cluster 1
        peaks_double <- double_peak_curves_analysis %>% filter(n_odd_peaks > 1) %>% ungroup() %>% summarise(unique(strain)) # double peaks in all replicates
        sub_ts[which(sub_ts$strain %in% unlist(peaks_double)), "cluster"] = "double"
        sub_parm[which(sub_parm$strain %in% unlist(peaks_double)), "cluster"] = "double"
      }
      
      ####**** (Cluster 2) clear "normal peaks"
      normal_peak_curves_analysis <- as.data.frame(sub_parm %>% mutate(odd = odd_peaks + odd_width + odd_shoulder + odd_shoulder_past) %>% filter(odd == 0) %>% 
                                                     group_by(strain) %>% dplyr::mutate(n_normal = n()))
      nd <- left_join(as.data.frame(sub_ts), normal_peak_curves_analysis[,c("strain","rep","n_normal")], by = c("strain", "rep"))
      g <- ggplot(nd, aes(x=Time, y = value_J, group = rep)) + geom_line(aes(col = factor(n_normal))) + facet_wrap(~strain)+ 
        scale_color_manual("Number of\nreps with\nnormal peaks", breaks = c(3,2,1,"NA"), labels = c(3,2,1,0), values = c("red","orange","blue","grey"))
      ggsave(paste0(plot_where,drytimes[j],"_",inocs[k],"_",name,"_normal_peaks.pdf"))
      
      peaks_normal <- normal_peak_curves_analysis %>% filter(n_normal > 1) %>% ungroup() %>% summarise(unique(strain)) # double peaks in two or all replicates
      sub_ts[which(sub_ts$strain %in% unlist(peaks_normal)), "cluster"] = "normal"
      sub_parm[which(sub_parm$strain %in% unlist(peaks_normal)), "cluster"] = "normal"
      
      ####**** (Cluster 3) Some have a single rep that is a "peak" and then rest are "shoulders" => basically multiple peaks = spike
      spiked_analysis <- as.data.frame(sub_parm %>% mutate(odd_spike = odd_peaks + odd_shoulder) %>% filter(odd_spike > 0) %>% group_by(strain) %>% dplyr::mutate(n_spike = n())) # how many have spikes in all? 
      nd <- left_join(as.data.frame(sub_ts), spiked_analysis[,c("strain","rep","n_spike")], by = c("strain", "rep"))
      g <- ggplot(nd, aes(x=Time, y = value_J, group = rep)) + geom_line(aes(col = factor(n_spike))) + facet_wrap(~strain)+ 
        scale_color_manual("Number of\nreps with\nspike peaks", breaks = c(3,2,1,"NA"), labels = c(3,2,1,0), values = c("red","orange","blue","grey"))
      ggsave(paste0(plot_where,drytimes[j],"_",inocs[k],"_",name,"_spike_peaks.pdf"))
      
      peaks_spike <- spiked_analysis %>% filter(n_spike > 1) %>% filter(!strain %in% unlist(peaks_double)) %>% ungroup() %>% summarise(unique(strain)) # 2 or more have spike
      sub_ts[which(sub_ts$strain %in% unlist(peaks_spike)), "cluster"] = "spike"
      sub_parm[which(sub_parm$strain %in% unlist(peaks_spike)), "cluster"] = "spike"
      
      ## Those with assignation
      clustered = unlist(sub_ts %>% filter(!cluster == "") %>% summarise(unique(strain)))
      
      ####***** (Cluster 4) Clear post shoulders
      post_should_analysis <- as.data.frame(sub_parm %>% filter(odd_shoulder_past == 1) %>% group_by(strain, inoc, drytime) %>% dplyr::mutate(n_spast= n())) # how many have spikes in all? 
      nd <- left_join(as.data.frame(sub_ts), post_should_analysis[,c("strain","rep","n_spast")], by = c("strain", "rep"))
      g <- ggplot(nd, aes(x=Time, y = value_J, group = rep)) + geom_line(aes(col = factor(n_spast))) + facet_wrap(~strain)+ 
        scale_color_manual("Number of\nreps with\npast peaks", breaks = c(3,2,1,"NA"), labels = c(3,2,1,0), values = c("red","orange","blue","grey"))
      ggsave(paste0(plot_where,drytimes[j],"_",inocs[k],"_",name,"_shoulder_past.pdf")) 
      
      peaks_post_shoulder <- post_should_analysis %>% filter(n_spast > 1) %>% filter(!strain %in% clustered) %>% ungroup() %>% summarise(unique(strain)) # 2 or more have post shoulders
      sub_ts[which(sub_ts$strain %in% unlist(peaks_post_shoulder)), "cluster"] = "post_shoulder"
      sub_parm[which(sub_parm$strain %in% unlist(peaks_post_shoulder)), "cluster"] = "post_shoulder"
      
      ## Those with assignation (update)
      clustered = unlist(sub_ts %>% filter(!cluster == "") %>% summarise(unique(strain)))
      
      # (x) Clear pre shoulders - not used atm
      #pre_should_analysis <- as.data.frame(sub_parm %>% filter(odd_shoulder == 1) %>% group_by(strain) %>% mutate(n_spre= n())) # how many have spikes in all? 
      #ggplot(pre_should_analysis %>% group_by(strain) %>% slice(1), aes(n_spre)) + geom_histogram(binwidth = 1)
      #nd <- left_join(as.data.frame(sub_ts), pre_should_analysis[,c("strain","rep","n_spre")], by = c("strain", "rep"))
      #g <- ggplot(nd, aes(x=Time, y = value_J, group = rep)) + geom_line(aes(col = factor(n_spre))) + facet_wrap(~strain)
      #ggsave(paste0(plot_where,name,"_shoulder_pre.pdf"))
      
      # (Cluster 5) Wide
      wide_analysis <- as.data.frame(sub_parm %>% filter(odd_width == 1) %>% group_by(strain, inoc, drytime) %>% dplyr::mutate(n_wide = n())) # how many have spikes in all? 
      nd <- left_join(as.data.frame(sub_ts), wide_analysis[,c("strain","rep","n_wide")], by = c("strain", "rep"))
      g <- ggplot(nd, aes(x=Time, y = value_J, group = rep)) + geom_line(aes(col = factor(n_wide))) + facet_wrap(~strain) + 
        scale_color_manual("Number of\nreps with\nwide peaks", breaks = c(3,2,1,"NA"), labels = c(3,2,1,0), values = c("red","orange","blue","grey"))
      ggsave(paste0(plot_where,drytimes[j],"_",inocs[k],"_",name,"_wide.pdf"))
      
      peaks_wide <- wide_analysis %>% filter(n_wide > 1) %>% filter(!strain %in% clustered) %>% ungroup() %>% summarise(unique(strain)) # 2 or more are wide
      sub_ts[which(sub_ts$strain %in% unlist(peaks_wide)), "cluster"] = "wide"
      sub_parm[which(sub_parm$strain %in% unlist(peaks_wide)), "cluster"] = "wide"
      
      ## Those with assignation (update)
      clustered = unlist(sub_ts %>% filter(!cluster == "") %>% summarise(unique(strain)))
      
      ### Those not clustered
      odd_normal_peak_curves_analysis <- as.data.frame(sub_parm %>% mutate(odd = odd_peaks + odd_width + odd_shoulder + odd_shoulder_past, 
                                                                           odd_label = paste(odd_peaks,odd_width,odd_shoulder,odd_shoulder_past)))
      nd <- left_join(as.data.frame(sub_ts %>% filter(cluster =="")), odd_normal_peak_curves_analysis[,c("strain","rep","odd","odd_label")], by = c("strain", "rep"))
      if(dim(nd)[1]!=0){
        g <- ggplot(nd, aes(x=Time, y = value_J, group = rep)) + geom_line(aes(col = factor(odd_label))) + facet_wrap(~strain) + 
          scale_color_discrete("peak / width / shoulder / shoulder past")
        ggsave(paste0(plot_where,drytimes[j],"_",inocs[k],"_",name,"_not_clustered.pdf"))
      }
      
      # Wide assign to shoulder past - analysis of 0t / 5inoc suggest most had one past shoulder / one wide and often wide == really past shoulder
      make_shoulder_past <- odd_normal_peak_curves_analysis %>% filter(odd_label == "0 1 0 0") %>% filter(!strain %in% clustered) %>% ungroup() %>% summarise(unique(strain))
      sub_ts[which(sub_ts$strain %in% unlist(make_shoulder_past)), "cluster"] = "post_shoulder"
      sub_parm[which(sub_parm$strain %in% unlist(make_shoulder_past)), "cluster"] = "post_shoulder"
      
      ### Group level plots
      cbp = c("#D55E00","#E69F00","#0072B2","#56B4E9","#F0E442","#999999")
      g <- ggplot(sub_ts, aes(x=Time, y = value_J, group = rep)) + geom_line(aes(col = factor(cluster))) + facet_wrap(~strain) + 
        scale_color_manual("Cluster", breaks = c("double","spike","normal","post_shoulder","wide",""), labels = c("double","spike","normal","post_shoulder","wide","none"), 
                           values = cbp) + 
        ggtitle(paste0("Drytime = ",drytimes[j],", Inoculum = ",inocs[k]))
      ggsave(paste0(plot_where,drytimes[j],"_",inocs[k],"_",name,"_grouped.pdf"))
      
      g <- ggplot(sub_ts %>% ungroup() %>% arrange(desc(cluster)), aes(x=Time, y = value_J, group = interaction(strain,rep))) + geom_line(aes(col = factor(cluster))) + 
        facet_grid(cluster~.) +
        scale_color_manual("Cluster", breaks = c("double","spike","normal","post_shoulder","wide",""), labels = c("double","spike","normal","post_shoulder","wide","none"), 
                           values = cbp)
      ggsave(paste0(plot_where,drytimes[j],"_",inocs[k],"_",name,"_grouped_rows.pdf"))
      
      
      
      ### Bind together for now 
      parameters_end <- rbind(parameters_end, sub_parm)
      ts_end <- rbind(ts_end, sub_ts)
      
      print(dim(ts_end))
      
    }
  }
  
  
  ## Return with assignation
  return(list(ts = ts_end, parameters = parameters_end))
}




##### Function to explore, normalise and extract non-growth output from OD vs CS data comparison 
data_exploration_od_cs <- function(datae, inoc_name, normal_strain = "11257"){
  # datae = time series to explore
  # inoc_name = which inoculum is being explored
  # normal_strain = strain which has the normal curve - take max of this for the normalisation 
  
  # remove all data before the time point that they all have which is the max of the minimum times to avoid odd completion curves prior to start and after
  cutoff_time_dn = max(datae %>% dplyr::group_by(rep, exp, inoc) %>% dplyr::summarise(minn = min(Time)) %>% ungroup() %>% dplyr::select(minn))
  cutoff_time_up = min(datae %>% dplyr::group_by(rep, exp, inoc) %>% dplyr::summarise(maxx = max(Time)) %>% ungroup() %>% dplyr::select(maxx))
  
  # Plot raw data
  ggplot(datae, aes(x=Time, y = value, group = interaction(inoc, exp, strain, rep))) + geom_line(aes(col = factor(rep))) + 
    facet_grid(exp~strain, scales = "free")
  ggsave(paste0("plots/ODvsC2_", inoc_name,"_rawdata.pdf"))
  
  data_od <- datae %>% filter(Time > cutoff_time_dn, Time < cutoff_time_up) %>% 
    dplyr::select(Time, rep, exp, value, strain, inoc) %>% 
    filter(inoc == 5) %>% # only 10^5 for this analysis
    group_by(rep, exp, strain) %>%
    mutate(ma_value = zoo::rollapply(value, 5, mean,fill = NA),
           differ = c(0,diff(ma_value)),
           compara = ifelse(exp == "CS", value, differ)) %>%
    ungroup() 
  
  g2 <- ggplot(data_od, aes(x=Time, y = ma_value, group = interaction(rep, exp, strain))) + 
    geom_line(aes(col = factor(rep))) + 
    facet_grid(exp~strain,scales = "free") + 
    scale_x_continuous("Time (h)") 
  ggsave(paste0("plots/ODvsCS_", inoc_name,"_data_smoothed.pdf"))
  
  g3a <- ggplot(data_od, aes(x=Time, y = compara, group = interaction(rep, exp, strain))) + 
    geom_line(aes(col = factor(rep))) + 
    facet_grid(exp~strain, scales = "free") + 
    scale_x_continuous("Time (h)") + 
    scale_y_continuous("Raw data") + 
    scale_color_discrete("Replicate")
  ggsave(paste0("plots/ODvsCS_",inoc_name,"_data_compare.pdf"))
  
  
  ##### Normalisation 
  # Try to plot together? Need to normalise... 
  ## First tried scale by max. Normalise - but then don't get comparison between strains
  # data_od <- data_od %>% group_by(strain, rep, exp) %>% mutate(max_v = max(compara, na.rm = TRUE),
  #                                                              compara_norm = compara / max_v)
  
  # Second find strain that is normal - divide by max of that then subtract
  #strains <- c("11016","11051", "11210", "11257")
  #clusters <- c("spike", "double", "spike", "normal") # categorizations from 2_cluster.R 
  
  max_vals_norm <- data_od %>% filter(strain == normal_strain) %>% group_by(rep, exp )%>% summarise(max_norms = max(compara, na.rm = TRUE))
  
  data_od <- left_join(data_od, max_vals_norm) %>% mutate(compara_norm = compara / max_norms)
  
  ggplot(data_od, aes(x=Time, y = compara_norm, group = interaction(rep, exp, strain))) + 
    geom_line(aes(col = factor(rep), lty = exp), lwd = 1) + 
    facet_wrap(~strain, scales = "free", ncol = 2) + 
    scale_x_continuous("Time (h)") + 
    scale_color_discrete("Replicate") + 
    scale_linetype("Data")
  ggsave(paste0("plots/ODvsCS_",inoc_name,"_data_compare_norm.pdf"))
  
  ##### Difference in output 
  # Subtract normalised data? Need to complete: measured at different time points
  data_od_normd <- data_od %>% ungroup() %>% filter(Time > cutoff_time_dn, Time < cutoff_time_up) %>% # make sure cover same time for all 
    complete(rep, strain, exp, Time) %>% mutate(compara_norm_inp = na_ma(compara_norm, k = 4, weighting = "linear", maxgap = 10)) %>% # fill in all time points and then linear imputation between (tried exponential and simple but get more odd bumps)
    group_by(strain, rep, exp) %>% dplyr::select(Time, strain,rep, exp, compara_norm_inp) %>% # Take imputed values
    pivot_wider(id_cols = c(strain, Time, rep), names_from = exp, values_from = compara_norm_inp) %>% mutate(nongrowth_only = CS - OD) %>% # look for difference between OD and heat output
    pivot_longer(cols = c("CS","OD"), names_to = "exp", values_to ="imput_val")
  
  data_od_normd_ana <- left_join(data_od, data_od_normd)
  
  ggplot(data_od_normd_ana, aes(x=Time, y = nongrowth_only, group = interaction(exp,rep, strain))) + 
    geom_line(aes(col = factor(rep)), lwd = 1) + 
    geom_line(aes(y = imput_val,col = interaction(rep))) + 
    facet_wrap(~strain, scales = "free", ncol = 2) + 
    scale_x_continuous("Time (h)") 
  ggsave(paste0("plots/ODvsCS_",inoc_name,"_data_nongrowth_togplot.pdf"))
  
  g1 <- ggplot(data_od_normd_ana, aes(x=Time, y = nongrowth_only, group = interaction(rep, strain))) + 
    geom_line(aes(col = interaction(rep))) + 
    facet_wrap(~strain, ncol = 4) + 
    scale_x_continuous("Time (h)") + 
    scale_y_continuous("Non growth only") + 
    scale_color_discrete("Replicate") + 
    geom_hline(yintercept = c(0,0.5), alpha = 0.2) #+ 
  #geom_vline(xintercept = c(25000,30000))
  
  g2 <- ggplot(data_od_normd_ana, aes(x=Time, group = interaction(exp,rep, strain))) + 
    geom_line(aes(y = imput_val,col = interaction(rep), linetype = exp)) + 
    facet_wrap(~strain, ncol = 4) + 
    scale_x_continuous("Time (h)") + 
    scale_y_continuous("Normalised measure") + 
    scale_color_discrete("Data type and\nreplicate") + 
    scale_linetype("Data type")
  
  g3a / g2 / g1 
  ggsave(paste0("plots/ODvsCS_",inoc_name,"_data_nongrowth_tog_grid.pdf"), width = 20, height = 20)
  
}


#### Simple data extraction for non-heat curves

charac_extract <- function(ts, Time, value, ts_orig_1, ts_orig_2, value_orig, name4output){
  ## ts = timeseries
  ## Time = name of time column
  ## value = name of value column
  ## name4output = strain, replicate, condition, inocl = labels for output
  
  ## (1) Fit spline
  ## This gives lag time and exponential growth rate cumulative
  gc_fit <- gcFitSpline(ts[,Time], ts[,value])
  # parameters from this fit
  s <- summary(gc_fit)
  
  ## (2) What is the maximum heat flow and when?
  wmax <- which.max( unlist(ts[,value])[-c(1:6)]) + 6 # take off funny initial behaviour < 2hr (take off in which.max and then add on 6 for position as taken 6 off)
  time_max_heat_flow <- as.numeric(ts[wmax,Time])
  value_max_heat_flow <- as.numeric(ts[wmax,value])
  
  ## (3) What is the minimum heat flow and when?
  wmin <- which.min( unlist(ts[,value])[-c(1:6)]) + 6 # take off funny initial behaviour < 2hr (take off in which.min and then add on 6 for position as taken 6 off)
  time_min_heat_flow <- as.numeric(ts[wmin,Time])
  value_min_heat_flow <- as.numeric(ts[wmin,value])
  
  ## (4) AUC - calculate original and then subtract for non-growth heat
  ## Orig 1
  gc_fit <- gcFitSpline(ts_orig_1[,Time], ts_orig_1[,value_orig])
  # parameters from this fit
  s1 <- summary(gc_fit)
  
  # From lag time to 10 hrs after
  lag_time <- as.numeric(which.min(abs(unlist(ts_orig_1[,Time]) - s1$lambda.spline)))
  ten_after = ifelse(max(ts_orig_1[,Time])<100,10,10*60*60) # units in current data either hours or seconds (latter is 10 hrs of seconds)
  ten_after_lag = as.numeric(which.min(abs(unlist(ts_orig_1[,Time]) - (s1$lambda.spline + ten_after))))
  gc_fit_auc <- gcFitSpline(ts_orig_1[c(lag_time:ten_after_lag),Time], ts_orig_1[c(lag_time:ten_after_lag),value_orig])
  # parameters from this fit
  s_auc_1 <- summary(gc_fit_auc)
  
  ## Orig 2
  if(sum(ts_orig_2[,value_orig])!=0){
    gc_fit <- gcFitSpline(ts_orig_2[,Time], ts_orig_2[,value_orig])
    # parameters from this fit
    s2 <- summary(gc_fit)
    
    # From lag time to 10 hrs after
    lag_time <- as.numeric(which.min(abs(unlist(ts_orig_2[,Time]) - s2$lambda.spline)))
    ten_after = ifelse(max(ts_orig_2[,Time])<100,10,10*60*60) # units in current data either hours or seconds (latter is 10 hrs of seconds)
    ten_after_lag = as.numeric(which.min(abs(unlist(ts_orig_2[,Time]) - (s2$lambda.spline + ten_after))))
    gc_fit_auc <- gcFitSpline(ts_orig_2[c(lag_time:ten_after_lag),Time], ts_orig_2[c(lag_time:ten_after_lag),value_orig])
    # parameters from this fit
    s_auc_2 <- summary(gc_fit_auc)
    
    ### Auc non-growth 
    s_auc <- s_auc_2$integral.spline - s_auc_1$integral.spline}else{s_auc <- s_auc_1$integral.spline} # if running on just heat flow then second dataset = 0
  
  ### Simple parameter extract 
  para_simple <- c(time_max_heat_flow, value_max_heat_flow, time_min_heat_flow, value_min_heat_flow, s_auc, s$lambda.spline)
  
  return(para_simple)
}
