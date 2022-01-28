##### GROFIT FUNCTIONS FROM NO LONGER SUPPORT PACKAGE:
## ***** Citation: ***** #################
## Matthias Kahm, Guido Hasenbrink, Hella Lichtenberg-Frate, Jost Ludwig, Maik Kschischo (2010). 
## grofit: Fitting Biological Growth Curves with R. Journal of Statistical Software, 33(7), 1-21. 
## URL http://www.jstatsoft.org/v33/i07/. 
## **********************************

### Fit a Spline to curve: output relevant parameters
gcFitSpline <-
  function(time,data, gcID ="undefined", control=grofit.control())
  {
    
    # /// check input parameters
    if (is(control)!="grofit.control") stop("control must be of class grofit.control!")
    if (!control$fit.opt%in%c("s","b")) stop("Fit option is not set for a spline fit. See grofit.control()")
    
    # /// conversion to handle even data.frame inputs
    time <- as.vector(as.numeric(as.matrix(time)))
    data <- as.vector(as.numeric(as.matrix(data)))
    
    # /// check length of input data
    if (length(time)!=length(data)) stop("gcFitSpline: length of input vectors differ!")
    
    # /// determine which values are not valid
    #bad.values <-  (is.na(time))|(time<0)|(is.na(data))|(data<0)|(!is.numeric(time))|(!is.numeric(data))
    bad.values <-  (is.na(time))|(time<0)|(is.na(data))|(!is.numeric(time))|(!is.numeric(data)) # OK with negative data
    
    # /// remove bad values or stop program
    if (TRUE%in%bad.values)
    {
      if (control$neg.nan.act==FALSE)
      {
        time    <- time[!bad.values]
        data    <- data[!bad.values]
      }
      else{
        stop("Bad values in gcFitSpline")
      }
    }
    
    if (length(data)<5){
      cat("gcFitSpline: There is not enough valid data. Must have at least 5!")
      gcFitSpline <- list(raw.time = time, raw.data = data, gcID = gcID, fit.time = NA, fit.data = NA, parameters = list(A= NA, mu=NA, lambda=NA, integral=NA), parametersLowess=list(A=NA, mu=NA, lambda=NA), spline = NA, reliable=NULL, fitFlag=FALSE, control = control)
      class(gcFitSpline) <- "gcFitSpline"
      return(gcFitSpline)
    }
    else
    {
      # /// apply transformation
      if (control$log.x.gc==TRUE){time <- log(1+time)}
      if (control$log.y.gc==TRUE){ data   <- log(1+data)}
      
      #will be used as start value in nls
      halftime <- (min(time)+max(time))/2
      
      # spline fit and computation of the maximum derivative
      try(y.spl <- smooth.spline(time,data,spar=control$smooth.gc))
      if (is.null(y.spl)==TRUE){
        warning("Spline could not be fitted to data!")
        if (is.null(control$smooth.gc)==TRUE){
          cat("This might be caused by usage of smoothing parameter NULL\n")
          fit.nonpara        <- list(raw.x = time, raw.y = data, fit.x = NA, fit.y = NA, parameters = list(A= NA, mu=NA, lambda=NA, integral=NA), spline = NA, parametersLowess=list(A= NA, mu=NA, lambda=NA), spline = NA, reliable=NULL, fitFlag=FALSE, control = control)
          class(gcFitSpline) <- "gcFitSpline"
          return(gcFitSpline)
        }		    	
      }
      # spline fit
      dydt.spl   <- predict(y.spl, time, deriv = 1) # derivative of smoothed spline
      index      <- which.max(dydt.spl$y)          #index of maximum derivative
      t.max      <- dydt.spl$x[index] # time at maximum derivative
      dydt.max   <- max(dydt.spl$y) # Maximum derivative
      y.max      <- y.spl$y[index] # y value at max deriv
      mu.spl     <- dydt.max; # slope at maximum derivative = dydt.max
      b.spl      <- y.max-dydt.max*t.max           #intercept: y value at end of lag time (y = mx + b => y - mx = b) 
      lambda.spl <- -b.spl/mu.spl # 
      integral   <- low.integrate(y.spl$x,y.spl$y) #the integral under the curve = AUC
      
      #time_end_lag <- 
      
      # lowess fit
      low        <- lowess(time,data,f=0.25)
      y.low      <- low$y
      x.low      <- low$x
      dydt.low   <- diff(y.low)/diff(time)
      mu.low     <- max(dydt.low)
      index      <- which.max(dydt.low)            #index of maximum derivative
      t.max      <- x.low[index]
      y.max      <- y.low[index]
      b.low      <- y.max-mu.low*t.max             #intercept
      lambda.low <- (-1)*b.low/mu.low
      
      
    }
    
    gcFitSpline        <- list(raw.time = time, raw.data = data, gcID = gcID, 
                               fit.time = y.spl$x, fit.data = y.spl$y, 
                               parameters = list(A= max(y.spl$y), mu=mu.spl, lambda=lambda.spl, integral=integral, t.max = t.max), 
                               parametersLowess=list(A= max(y.low), mu=mu.low, lambda=lambda.low), 
                               spline = y.spl, reliable=NULL, fitFlag=TRUE, control = control)
    
    class(gcFitSpline) <- "gcFitSpline"
    gcFitSpline
    
  }


grofit.control <-
  function(
    neg.nan.act       = FALSE,
    clean.bootstrap   = TRUE ,
    suppress.messages = FALSE,
    fit.opt     = "b",
    log.x.gc    = FALSE,
    log.y.gc    = FALSE,
    interactive = TRUE,
    nboot.gc    = 0,
    smooth.gc   = NULL,
    model.type  = c("logistic", "richards", "gompertz", "gompertz.exp"),
    have.atleast   = 6,
    parameter      = 9,
    smooth.dr      = NULL,
    log.x.dr       = FALSE,
    log.y.dr       = FALSE,
    nboot.dr      = 0)
  {
    if ((is.character(fit.opt)==FALSE)|(length(fit.opt)!=1))
      stop("value of fit.opt must be character and of one element")
    
    if (is.character(model.type)==FALSE)
      stop("value of model.type must be character")
    
    if ((is.logical(neg.nan.act)==FALSE)|(length(neg.nan.act)!=1))
      stop("value of neg.nan.act must be logical and of one element")
    
    if ((is.logical(clean.bootstrap)==FALSE)|(length(clean.bootstrap)!=1))
      stop("value of clean.bootstrap must be logical and of one element")
    
    if ((is.logical(suppress.messages)==FALSE)|(length(suppress.messages)!=1))
      stop("value of suppress.messages must be logical and of one element")
    
    if ((is.logical(log.x.gc)==FALSE)|(length(log.x.gc)!=1))
      stop("value of log.x.gc must be logical and of one element")
    
    if ((is.logical(log.y.gc)==FALSE)|(length(log.y.gc)!=1))
      stop("value of log.y.gc must be logical and of one element")
    
    if ((is.logical(interactive)==FALSE)|(length(interactive)!=1))
      stop("value of interactive must be logical and of one element")
    
    if ((is.logical(log.x.dr)==FALSE)|(length(log.x.dr)!=1))
      stop("value of log.x.dr must be logical and of one element")
    
    if ((is.logical(log.y.dr)==FALSE)|(length(log.y.dr)!=1))
      stop("value of log.y.dr must be logical and of one element")
    
    if ((is.numeric(nboot.gc)==FALSE)|(length(nboot.gc)!=1)|(nboot.gc<0))
      stop("value of nboot.gc must be numeric (>=0) and of one element")
    
    if ((is.numeric(have.atleast)==FALSE)|(length(have.atleast)!=1)|(have.atleast<6))
      stop("value of have.atleast must be numeric (>=6) and of one element")
    
    if ((is.numeric(parameter)==FALSE)|(length(parameter)!=1))
      stop("value of parameter must be numeric and of one element")
    
    if ((is.numeric(nboot.dr)==FALSE)|(length(nboot.dr)!=1)|(nboot.dr<0))
      stop("value of nboot.dr must be numeric (>=0) and of one element")
    
    if (((is.numeric(smooth.gc)==FALSE) && (is.null(smooth.gc)==FALSE)))
      stop("value of smooth.gc must be numeric or NULL")
    
    if (((is.numeric(smooth.dr)==FALSE) && (is.null(smooth.dr)==FALSE)))
      stop("value of smooth.dr must be numeric or NULL")
    
    
    grofit.control <- list(neg.nan.act=neg.nan.act, clean.bootstrap=clean.bootstrap, suppress.messages=suppress.messages,
                           fit.opt=fit.opt, log.x.gc=log.x.gc, log.y.gc=log.y.gc, interactive=interactive,
                           nboot.gc=round(nboot.gc), smooth.gc=smooth.gc, smooth.dr=smooth.dr,
                           have.atleast=round(have.atleast), parameter=round(parameter), log.x.dr=log.x.dr, log.y.dr=log.y.dr,
                           nboot.dr=round(nboot.dr), model.type=model.type)
    class(grofit.control) <- "grofit.control"
    grofit.control
  }


#### Function to get auc 
low.integrate <-
  function (x,y)
  {
    
    if(is.vector(x)==FALSE || is.vector(y)==FALSE)
      stop("low.integrate: two vectors x and y are needed !")
    if(length(x)!=length(y))
      stop("low.integrate: x and y have to be of same length !")
    
    
    #spline fit and computation of the maximum derivative
    y.spl <- NULL
    try(y.spl <- smooth.spline(x,y))
    if (is.null(y.spl)==TRUE){
      warning("Spline could not be fitted to data!")
      stop("Error in low.integrate")	
    }
    
    f     <- function(t){
      p<-predict(y.spl,t)
      f<-p$y
    }
    
    low.integrate <- integrate(f,min(x),max(x))$value
    
  }

summary.gcFitSpline <-
  function(object,...)
  {
    
    # object of class gcFitSpline
    
    contents.fitted.spline  <- c("mu.spline", "lambda.spline", "A.spline", "integral.spline")
    
    if ((is.na(object$fitFlag)==TRUE)|(object$fitFlag==FALSE)){
      table<-rep(NA,length(contents.fitted.spline))
    }
    else{
      table <- c(object$parameters$mu, object$parameters$lambda,  object$parameters$A, object$parameters$integral)
    }
    
    table               <- data.frame(t(table))
    colnames(table)     <- contents.fitted.spline
    summary.gcFitSpline <- table
    
  }