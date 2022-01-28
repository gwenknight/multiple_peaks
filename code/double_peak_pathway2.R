#### Function

#### From Quentin: optim

startl=list(a=unlist(ts[peaks_index,value])[1]*1.5, 
            b=diff(interval_peak)/3, 
            c=unlist(ts[peaks_index,Time])[1],
            d=unlist(ts[peaks_index,value])[1]*1.5, 
            e=diff(interval_peak), 
            f=unlist(ts[peaks_index,Time])[1])

startl = list(a = 150, b = 0.2, 
              c = 0.1, d = 0.1, e = 0.1, f = 0.1)

startl = list(peak1 = 792, peak2 = 1200, sig = 0.4, mu = 2.5, sig2 = 0.1, mu2 = 1.9)

#wrapper function which takes as input parameters and data
#calculates whatever needed with params (in this case, bacteria growth curve), and returns likelihood
fit_0 = function(p, data_0){
  #p = parameters 
  peak1 = as.numeric(p["a"])
  sig = as.numeric(p["b"])
  mu = as.numeric(p["c"])
  peak2 = as.numeric(p["d"])
  sig2 = as.numeric(p["e"])
  mu2 = as.numeric(p["f"])
  # data_0 = data to fit
  #quick_model is just a deSolve wrapper that models bacterial growth based on mu and Bmax. I also feed it the starting value for bacteria, and the times vector.
  x = data_0[,"Time"]
  y <- data_0[,"differ"] - as.numeric(data_0[1,"differ"])
  two_logs = peak1/(x*sig*sqrt(2*pi))*exp(-(log(x) - mu)^2/(2*sig^2)) + peak2/(x*sig2*sqrt(2*pi))*exp(-(log(x) - mu2)^2/(2*sig2^2))
  #then I just calculate least squares between model and data
  sum(((y - two_logs))^2)
  #sum(((y - two_logs)/two_logs)^2)
}

fit_0(startl, ts)

#optim needs as input starting values for your parameters, the function to optimise, and any other argument that function needs (in this case, the data: data_0)
growth_par = optim(startl, fit_0, data_0 = ts)

p <- as.data.frame(t(growth_par$par))

two_logs = p$peak1/(x*p$sig*sqrt(2*pi))*exp(-(log(x) - p$mu)^2/(2*p$sig^2)) + p$peak2/(x*p$sig2*sqrt(2*pi))*exp(-(log(x) - p$mu2)^2/(2*p$sig2^2))
g1l <-  p$peak1/(x*p$sig*sqrt(2*pi))*exp(-(log(x) - p$mu)^2/(2*p$sig^2)) 
g2l <-  p$peak2/(x*p$sig2*sqrt(2*pi))*exp(-(log(x) - p$mu2)^2/(2*p$sig2^2))
two_logs = g1l + g2l
plot(x,y)
lines(x, two_logs)
lines(x, g1l, type = "l", col = "red")
lines(x, g2l, type = "l", col = "blue")


p$peak1 = 0.005
p$sig = 0.4
p$mu = 2.5
p$peak2 = 0.005
p$sig2 = 0.22
p$mu2 = 2.8

growth_par = optim(p, fit_0, data_0 = ts)

p <- as.data.frame(t(growth_par$par))

g1l <-  p$peak1/(x*p$sig*sqrt(2*pi))*exp(-(log(x) - p$mu)^2/(2*p$sig^2)) 
g2l <-  p$peak2/(x*p$sig2*sqrt(2*pi))*exp(-(log(x) - p$mu2)^2/(2*p$sig2^2))
two_logs = g1l + g2l
plot(x,y, ylim = c(0,0.01))
lines(x, two_logs)
lines(x, g1l, type = "l", col = "red")
lines(x, g2l, type = "l", col = "blue")




### Input: times series data
### Output: time series data up to cut point: time of end of exponential growth
### Output: lag time / exponential growth / time of end of exponential growth 

cut_extract_p2 <- function(ts, Time, value, name4output, thresh_wide = 90, plot = 0, plot_where = "plots/", early_cut = 3.5){
  ## ts = timeseries
  ## Time = name of time column
  ## value = name of value column
  ## name4output = strain, replicate, condition, inocl = labels for output
  ## thresh_wide = 90: % what is a wide peak? 
  ## plot = 0:  don't plot. 1 to plot
  ## plot_where: location for files to output to  
  ## early_cut: for heat flow, remove the first 3.5hrs of data 
  
  ## is this strain replicate ODD? Set all ODD indicators to zero initially
  odd_peak <- 0 # 0 = one peak
  odd_width <- 0 # 0 = not a wide peak
  odd_shoulder <- 0 # 0 = no shoulder
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
      
      if(length(close_peaks_i) > 1){odd_peak <- 1} # if multiple close time and height peaks
      
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
          if(length(fpsq) != 0){ # of no preaks then not a clear shoulder
            shoulder_point1 = min(time_peaks[1]-1,max(times_line[fpsq]))
            ws <- which(round(ts[,Time],5) == round(shoulder_point1,5)); shoulder_point_v1 <- ts[ws,value]
            height <- shoulder_point_v1 / ts[peaks_index[1],value]
            if(height > 0.5 && height < 0.96){ # if the shoulder is greater than half way up but not too close
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
  if(length(peaks_index)!=0 && !is.na(peaks_index)){
    if(length(time_peaks_diff) > 1 & odd_shoulder == 0){ # if multiple peaks but no shoulder
      shoulder_point1 = min(time_peaks)
      shoulder_point_v1 = ts[which(round(ts[,Time],5) == round(shoulder_point1,5)),value]
      if((shoulder_point_v1 / ts[peaks_index[1],value]) > 0.5){ # shoulder needs to be high still! 
        shoulder_point = shoulder_point1; shoulder_point_v = shoulder_point_v1 # Update shoulder values to ones that are > half way
      }
    }
  }
  
  
  ##### (5) Double curves? Only use to count how many have this - fit not explored here
  x <- ts[,Time]
  y <- ts[,value]
  
  # Guesses for the parameters in the curve
  startl=list(a=unlist(ts[peaks_index,value])[1]*1.5, 
              b=diff(interval_peak)/3, 
              c=unlist(ts[peaks_index,Time])[1],
              d=unlist(ts[peaks_index,value])[1]*1.5, 
              e=diff(interval_peak), 
              f=unlist(ts[peaks_index,Time])[1])
  # Set indicator/output parameters to zero
  fit1 <- 0; fit2 <- 0; fit3 <- 0;
  plot_p <- NA; plot_p2 <- NA; plot_p3 <- NA; plot_p4 <- NA;
  p <- 0
  # Try to fit two normal curves using non-linear least squares methods 
  try(fit1 <- nls(y~(a/b)*exp(-(x-c)^2/(2*b^2))+(d/e)*exp(-(x-f)^2/(2*e^2)),start = startl, algorithm="port"), silent = TRUE)
  try(fit2 <- nls(y~peak1/(x*sig*sqrt(2*pi))*exp(-(log(x) - mu)^2/(2*sig^2)) + peak2/(x*sig2*sqrt(2*pi))*exp(-(log(x) - mu2)^2/(2*sig2^2)), 
                  start = list(peak1 = 0.03, peak2 = 0.03, sig = 0.2, mu = 2.4, sig2 = 0.1, mu2 = 2.8)), silent = TRUE)
  try(fit3 <- nls(y~(a/(2*b))*exp(-(x-c)^2/(2*b^2))+(a/(2*b))*exp(-(x-c/2)^2/(2*b^2)) +(d/e)*exp(-(x-f)^2/(2*e^2)),start = startl, algorithm="port"), silent = TRUE)
  
  try(fit4 <- t(optim(startl, fit_0, data_0 = ts)$par))
  

  
  # If manage to fit
  if(length(fit1) > 1){
    # print that can 
    print(paste("Double curve fit Norm", name4output,sep = " "))
    odd_double <- 1
    # Save the predictions of the fit 
    pred_p <- predict(fit1)
    # Save the coefficients of the fit
    p <- coef(fit1)
    # Pull out the actual curves using the parameters in p
    g1p <- (p["a"]/p["b"])*exp(-(x-p["c"])^2/(2*p["b"]^2)) 
    g2p <- (p["d"]/p["e"])*exp(-(x-p["f"])^2/(2*p["e"]^2))
    if(all(g1p > 0) && all(g2p>0)){
      # Store the curves
      plot_p <- as.data.frame(cbind(x,y,pred_p, g1p, g2p,0))
      colnames(plot_p) <- c("time","value","fit","curve1","curve2","curve3")
      # Store the parameters
      plot_p$a <- p["a"]; plot_p$b <- p["b"]; plot_p$c <- p["c"]
      plot_p$d <- p["d"]; plot_p$e <- p["e"]; plot_p$f <- p["f"]
    }
  }
  
  
  # If manage to fit
  if(length(fit2) > 1){
    # print that can 
    print(paste("Double curve fit Lognorm", name4output,sep = " "))
    odd_double <- 1
    # Save the predictions of the fit 
    pred_p <- predict(fit2)
    # Save the coefficients of the fit
    p <- coef(fit2)
    # Pull out the actual curves using the parameters in p
    sig = p["sig"]; mu = p["mu"]; peak1 <- p["peak1"]
    sig2 <- p["sig2"]; mu2 <- p["mu2"]; peak2 <- p["peak2"]
    g1p <- peak1/(x*sig*sqrt(2*pi))*exp(-(log(x) - mu)^2/(2*sig^2))
    g2p <- peak2/(x*sig2*sqrt(2*pi))*exp(-(log(x) - mu2)^2/(2*sig2^2))
    if(all(g1p > 0) && all(g2p>0)){
      # Store the curves
      plot_p2 <- as.data.frame(cbind(x,y,pred_p, g1p, g2p,0))
      colnames(plot_p2) <- c("time","value","fit","curve1","curve2","curve3")
      # Store the parameters
      plot_p2$a <- p["sig"]; plot_p2$b <- p["mu"]; plot_p2$c <- p["peak1"]
      plot_p2$d <- p["sig2"]; plot_p2$e <- p["mu2"]; plot_p2$f <- p["peak2"]
    }
  }
  
  # If manage to fit
  if(length(fit3) > 1){
    # print that can 
    print(paste("Double curve fit Norm x3", name4output,sep = " "))
    odd_double <- 1
    # Save the predictions of the fit 
    pred_p <- predict(fit3)
    # Save the coefficients of the fit
    p <- coef(fit3)
    # Pull out the actual curves using the parameters in p
    g1p <- (p["a"]/(2*p["b"]))*exp(-(x-p["c"])^2/(2*p["b"]^2)) 
    g2p <- (p["a"]/(2*p["b"]))*exp(-(x-p["c"]/2)^2/(2*p["b"]^2)) 
    g3p <- (p["d"]/p["e"])*exp(-(x-p["f"])^2/(2*p["e"]^2))
    if(all(g1p > 0) && all(g2p>0) && all(g3p>0)){
      # Store the curves
      plot_p3 <- as.data.frame(cbind(x,y,pred_p, g1p, g2p, g3p))
      colnames(plot_p3) <- c("time","value","fit","curve1","curve2","curve3")
      # Store the parameters
      plot_p3$a <- p["a"]; plot_p3$b <- p["b"]; plot_p3$c <- p["c"]
      plot_p3$d <- p["d"]; plot_p3$e <- p["e"]; plot_p3$f <- p["f"]
    }
  }
  
  # If manage to fit
  if(length(fit4) > 1){
    print(c("fit4",fit4))
    # print that can 
    print(paste("Double curve optim Log normal x3", name4output,sep = " "))
    odd_double <- 1
    # Save the predictions of the fit 
    p <- as.data.frame(fit4)
    pred_p <- p$a/(x*p$b*sqrt(2*pi))*exp(-(log(x) - p$c)^2/(2*p$b^2)) + p$d/(x*p$e*sqrt(2*pi))*exp(-(log(x) - p$f)^2/(2*p$e^2))
    # Pull out the actual curves using the parameters in p
    g1l <-  p$a/(x*p$b*sqrt(2*pi))*exp(-(log(x) - p$c)^2/(2*p$b^2)) 
    g2l <-  p$d/(x*p$e*sqrt(2*pi))*exp(-(log(x) - p$f)^2/(2*p$e^2))
    
    if(all(g1l > 0) && all(g2l>0)){
      # Store the curves
      plot_p4 <- as.data.frame(cbind(x,y,pred_p, g1l, g2l,0))
      colnames(plot_p4) <- c("time","value","fit","curve1","curve2","curve3")
         
      # Store the parameters
      plot_p4$a <- p$peak1; plot_p4$b <- p$sig;  plot_p4$c <- p$mu
      plot_p4$d <- p$peak2; plot_p4$e <- p$sig2; plot_p4$f <- p$mu2
    }
  }
  
  
  ## Plot ts, cumulative and fit
  # if functions takes in a command to plot
  if(plot == 1){
    
    ### Fit to value : could add cumulative plot in here but not done
    datam <- reshape2::melt(ts[,c(Time,value)], id.vars = Time)
    gg <- ggplot(datam, aes(x=Time,y=value)) + geom_point() 
    ## model fit:
    gc_df <- as.data.frame(cbind(gc_fit$fit.time, gc_fit$fit.data))
    colnames(gc_df) <- c(Time,value)
    #gc_df[,value]<- c(gc_df[1,"csum"],diff(gc_df$csum)) # CHANGE TO CUMULATIVE? 
    #gc_dfm <- reshape2::melt(gc_df, id.vars = Time)
    ## add fit to data plot
    gg <- gg + geom_line(data = gc_df, aes(x=Time,y=value), col= "red")
    
    ## save
    ggsave(paste0(plot_where, name4output, "_model_fit.pdf"))
  } 
  
  ## Build vectors of required parameters to output
  param_o   <- c(time_max_heat_flow, value_max_heat_flow, 
                 s$mu.spline, s$lambda.spline,s$integral.spline, 
                 odd_peak, odd_width, max_level, odd_shoulder, odd_double, shoulder_point, as.numeric(shoulder_point_v))
  
  #print(param_o)  
  
  
  # (6) Now cut at shoulder or first peak
  if(shoulder_point > 0){
    cut_point_t <- shoulder_point; cut_point_v <- shoulder_point_v}else{
      cut_point_t <- time_max_heat_flow; cut_point_v <- value_max_heat_flow;
    }
  
  # NEW Time series generated: cut up to first peak or shoulder 
  #print(paste0("Early cut = ",early_cut))
  ts <- ts %>% 
    filter(Time > early_cut) %>% # cut off first 3hrs 
    mutate(cutpart = ifelse(Time <= cut_point_t,1,0)) %>%
    filter(cutpart == 1) # Trim off extra parts
  
  ### Look at this cut data
  ## Currently the time to use is: peak growth 
  timepeak = cut_point_t
  valpeak = cut_point_v[1]
  
  ## (7) What if there is a peak in this?
  peaks_index = find_peaks(ts[,value], m = 3)
  
  ## If there is a peak then reassign peak 
  if(length(peaks_index) > 0){
    w_early <- which(ts[peaks_index,Time] < 5) # Remove early ones
    if(length(w_early) > 0){peaks_index <- peaks_index[-w_early]}
    w_high <- which(ts[peaks_index,value] < 0.6*max(ts[,value])) # Remove low ones
    if(length(w_high) > 0){peaks_index <- peaks_index[-w_high]}
    if(length(peaks_index) > 0){ # if any later than 3 
      peaks_index = min(peaks_index)
      timepeak = ts[peaks_index,Time]
      valpeak = ts[peaks_index,value]
    }
  }
  
  # Reassign up to peak 
  ts = ts[which(ts[,Time] <= as.numeric(timepeak)),]
  
  ### (8) Check if the slope changes substantially in this period: may be a shoulder or a plateau near peak 
  if(dim(ts)[1]>4){ # If enough data
    st <- c()
    
    #plot(ts[,Time], ts[,value], "l")
    #points(ts[,Time], ts[,value])
    
    for(i in 2:(-1 + dim(ts)[1])){ # fit a linear model to the data in segments
      ts1 <- ts[max(1,i-6):i,]
      ts2 <- ts[(i):dim(ts)[1],]
      lm.1 <- lm(unlist(ts1[,value]) ~ unlist(ts1[,Time]))
      lm.2 <- lm(unlist(ts2[,value]) ~ unlist(ts2[,Time]))
      
      #lines(ts1[,Time], lm.1$coefficients[1] + lm.1$coefficients[2]*ts1[,Time], col = "blue")
      #lines(ts2[,Time], lm.2$coefficients[1] + lm.2$coefficients[2]*ts2[,Time], col = "red")
      
      st <- rbind(st, c(i,c(max(ts1[,value]), max(ts1[,Time]),lm.1$coefficients[2],lm.2$coefficients[2])))
    }
    
    st <- as.data.frame(st)
    colnames(st) <- c("i","maxval","maxtime","f_ang","s_ang") # first angle, second angle
    st$d <- 0; st$da <- 1000
    st$d[1:(dim(st)[1]-1)] = diff(st$s_ang) # look at change in second angle: want to know when substantial change
    st$da[1:(dim(st)[1]-1)] = abs(diff(st$s_ang)) # look at change in second angle: want to know when substantial change
    
    #plot(st$maxtime, st$d) 
    #plot(ts[,Time], ts[,value])
    
    ## Check if shoulder
    st_upper <- st%>% filter(maxtime > 0.50*max(st$maxtime)) # but want to be past halfway
    w1 <- which.min(st_upper$da) # as absolute this detects a plateau 
    w2 <- which.min(st_upper$da[-w1])
    timepeak_s = min(st_upper[w1,"maxtime"], st_upper[-w1,"maxtime"][w2]) # earliest 
    valpeak_s = as.numeric(ts[which(ts[,Time] == timepeak_s),value])
    
    ## Check not too low: if cut point already good enough. 
    # If OK then cut at shoulder value
    if(valpeak_s > 0.65*max(ts[,value])){
      valpeak <- valpeak_s; timepeak = timepeak_s}else{
        
        ## Check if end peak should be moved forward at all (i.e. a peak: want end of exponential growth )
        st_upper <- st %>% filter(maxval > 0.80*max(st$maxval)) # want to be near the end 
        #if(dim(st_upper)[1] > 5){ # have at least 4 d values to look at
        w1 <- which.max(st_upper$d)
        w2 <- which.max(st_upper$d[-w1])
        if(min(w1,w2) != 1){ # if its not just a slope down - if there is a plateau near the top? i.e. index now just the earliest point
          timepeak = st_upper[min(w1,w2),"maxtime"]
          valpeak = as.numeric(ts[which(ts[,Time] == timepeak),value])}
      }
  }
  
  if(plot == 1){
    g1 <- ggplot(ts,aes(x=Time, y= value)) + geom_line() + 
      geom_point(data = ts[which(round(ts[,Time],2) == as.numeric(round(timepeak,2))),c(Time,value)],col="red") + 
      geom_point(data= ts[1,], aes(x=shoulder_point, y = shoulder_point_v), col = "black") + 
      ggtitle(paste0(u[jj],"_",r[ii], "_",drying_times[kk],"_",q[ll]))
    dir.create(file.path(here(), paste0(plot_where,"/shoulder_curves")),showWarnings = FALSE)
    ggsave(paste0(plot_where,"shoulder_curves/cutpoint_highlighted_",name4output,".pdf")) 
  }
  
  tstopeak = ts[which(ts[,Time] <= as.numeric(timepeak)),c(Time,value)]
  
  ## (9) Growth Curve 
  ## Ths gives lag time and exponential growth rate 
  gc_fit <- gcFitSpline(tstopeak[,Time], tstopeak[,value])
  # parameters from this fit
  s <- summary(gc_fit)
  
  
  ## Build vectors of required parameters to output
  param_o   <- c(param_o, as.numeric(s$mu.spline), as.numeric(timepeak), as.numeric(valpeak))
  
  return(list(param = param_o, double_param = plot_p, double_param_logn = plot_p2, double_param_norm3 = plot_p3, double_param_op_lognorm2 = plot_p4))
  
}



###******* MODEL FITTING *************#######################################################################################################################################
## Fit to differentiated OD data

# strain names
data_od$rep <- "0";
data_od[,c("strain_label","inoc_name")]<-colsplit(data_od$variable, "(?<=\\p{L})(?=[\\d+$])", c("char", "digit"))

w<- which(data_od$strain_label == "A")
data_od[w,"rep"] = "A"
w<- which(data_od$strain_label == "B")
data_od[w,"rep"] = "B"
w<- which(data_od$strain_label == "D")
data_od[w,"rep"] = "A"
w<- which(data_od$strain_label == "E")
data_od[w,"rep"] = "B"

## What are the strains?
u <- as.character(unique(data_od$strain,stringsasFactors = FALSE)) # strains
## How many replicates? 
r <- unique(data_od$rep) # replicates
# What are the inoculums? 
q <- unique(data_od$inoc)

# Where the parameters for each strain are stored
param <- matrix(0, length(u)*length(q)*length(r), 18); # number of strains x number of replicates x number of experimental conditions
index <- 1 # for counting 
max <- c() # for calibration

## Run thru each strain/rep/drying time/ inoculum: fit a growth curve to each
for(jj in 1:length(u)){ # for each strain
  r <- unique(data_od %>% filter(strain == u[jj]) %>% dplyr::select(rep))[,1] # which replicates for this strain (all have different names)
  for(ii in 1:length(r)){ # for each replicate
    #for(kk in c(1,2)){ #each of the experimental conditions
    for(ll in 1:length(q)){ #each of the inocula
      
      strain <- u[jj];
      replicate <- as.character(r[ii])
      #   condition <- drying_times[kk]
      inocl <- q[ll]
      
      wi <- intersect(which(data_od$strain == strain),which(data_od$rep == replicate)) # if fit to each replicate
      #  wj <- intersect(wi, which(data_od$drytime == condition))
      w <- intersect(wi, which(data_od$inoc == as.numeric(inocl)))
      
      if(length(w) > 0){ # if this replicate exists for this strain (i.e. there is data)
        data1 <- data_od[w,] # Grab data
        wna<-which(is.na(data1$differ)) # remove NA at beginning
        data1 <- data1[-wna,]
        
        print(c(strain, replicate, inocl)) # output so can track how it is working
        
        p <- cut_extract_p2(data1, "Time", "differ", paste(strain, replicate, inocl,sep="_"),plot = 0, plot_where = "data_paper2/plots/path_1_", early_cut = 0) ### NEW function: runs on any timeseries: gives baseline and cut parameter analysis
        
        # print(p$param) # was getting errors generating 
        
        ## Required parameters
        param[index,] <- c(strain, replicate, inocl, p$param)
        index <- index + 1 # counting for storing matrix - next row each iteration
      }
    }
  }
}

# The original set 
param_orig <- param

