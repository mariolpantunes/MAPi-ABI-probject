library(zoo)
library(tseries)
library(TTR)
library(e1071)

###############################################################################
# User defined functions
###############################################################################

# Pre processing
# Applying log to minimize seasonally effects
pre.processing <- function(ts.data)
{
    ts.data.processed <- log(ts.data)

    return(ts.data.processed)
}

# Evaluation metrics

# Mean square error (MSE)

mse <- function(preds, trues)
{
    error <- mean((preds - trues)^2, na.rm = T)

    return(error)
}

# Mean absolute percentage error (MAPE)
mape <- function(preds, trues)
{
    error <- mean(abs((preds - trues)/trues),na.rm = T)

    return(error)
}

# Median absolute deviation (MAD)
mad <- function(preds, trues)
{
    error <- mean(abs(preds - trues), na.rm = T)

    return(error)
}

# Embedded Time series
embedded.dataset <- function(data, dimension = 3, tlag = 1)
{
    dataset <- data
    names <- c("fr")
    for (i in 1:(dimension-1))
    {
        dataset.tmp <- lag(data, -(i*tlag))
        dataset <- cbind(dataset, dataset.tmp)
        names <- cbind(names,paste("pr",i,sep = ""))
    }
    colnames(dataset) <- names
    #rows <- (dimension*tlag):(nrow(data))
    #dataset <- dataset[rows,]

    dataset <- na.omit(dataset)

    return(dataset)
}


# Monte Carlo simulation with sliding window
#
# Provided in KDD MapI class
MonteCarlo.estimates <- function(data, nreps=10, train.size, test.size, learner, learner.pars, perf.func, relearn.step=1, starting.points=NULL)
{
  data <- as.data.frame(data)
  n <- nrow(data)
  
  if (is.null(starting.points) || missing(starting.points))
  {
    selection.range <- (train.size+1):(n-test.size+1)
    starting.points <- sort(sample(selection.range,nreps))
  }
  preds <- list()
  
  for(it in seq(along=starting.points))
  {
    start <- starting.points[it]
    cat('Going for sliding window testing from row ',start,'(train.sz=',train.size,',Relearn=',relearn.step,')\n')
    preds[[it]] <- sliding.window.testing(data[1:(start+test.size-1),], train.size, learner, learner.pars,relearn.step, test.pos=start)
      
    if(exists('MC.results'))
      MC.results[it,] <- do.call(perf.func,list(preds[[it]],data[start:(start+test.size-1),],learner.pars))
    else
      MC.results <- data.frame(do.call(perf.func,list(preds[[it]],data[start:(start+test.size-1),],learner.pars)))
  }

  list(Starting.Points=starting.points,Results=MC.results,Predictions=preds)
}

sliding.window.testing <- function(orig.data, window.size, learner, learner.pars, relearn.step=1, test.pos=window.size+1)
{
    init.test <- test.pos
    n <- nrow(orig.data)
    preds <- vector()
    while (test.pos <= n)
    {
        cat('*')
        learner.pars$data <- orig.data[(test.pos-window.size):(test.pos-1),]
        model <- do.call(learner,learner.pars)
        preds <- c(preds,predict(model,orig.data[test.pos:min(n,test.pos+relearn.step-1),]))
        test.pos <- test.pos+relearn.step
    }
    cat('\n')

    return(preds)
}

target.name <- function(form, data)
{
    mt <- terms(form,data=data)
    if ((yvar <- attr(mt, "response")) <= 0)
    {
        stop(paste("Incorrect response variable",yvar))
    }
    as.character(attr(mt, "variables"))[-1][yvar]
}

error.func <- function(preds, test.data, learner.pars)
{
    form <- learner.pars[[1]]  # assumes the first argument is the formula
    list(mse = mse(preds,test.data[,target.name(form,test.data)]),
    mape = mape(preds,test.data[,target.name(form,test.data)]),
    mad = mad(preds,test.data[,target.name(form,test.data)]))
}

###############################################################################
# Code
###############################################################################

# Load data from CSV
raw.data <- read.csv("data.csv", col.names=c('Date','Value'), header = TRUE, stringsAsFactors = FALSE)

# Convert raw data into a TS object
ts.data <- ts(data = raw.data[2], frequency = 12, start = c(1983,1), end = c(2011,11))

# Plot time series
plot.ts(ts.data)

# Pre-processing time series
ts.data.processed <- pre.processing(ts.data)
#x11()
#plot.ts(ts.data.processed)

# Linear approach
cat("Linear approach\n")
# Moving average using differente time windows
cat("Moving average\n")
preds.window.3 <- SMA(ts.data.processed, 3)
cat("Window size = 3 (trimester)\n")
cat("MAPE:",mape(preds.window.3, ts.data.processed),"\n")
cat("MAD :",mad(preds.window.3, ts.data.processed),"\n")

preds.window.6 <- SMA(ts.data.processed, 6)
cat("Window size = 6 (semester)\n")
cat("MAPE:",mape(preds.window.6, ts.data.processed),"\n")
cat("MAD :",mad(preds.window.6, ts.data.processed),"\n")

preds.window.12 <- SMA(ts.data, 12)
cat("Window size = 12 (annually)\n")
cat("MAPE:",mape(preds.window.12, ts.data.processed),"\n")
cat("MAD :",mad(preds.window.12, ts.data.processed),"\n")

# Exponencial moving average 
cat("Exponencial moving average\n")
preds <- EMA(ts.data.processed)
cat("MSE :",mse(preds, ts.data.processed),"\n")
cat("MAPE:",mape(preds, ts.data.processed),"\n")
cat("MAD :",mad(preds, ts.data.processed),"\n")

# Generate simple dataset
dataset <- embedded.dataset(ts.data,13,1)

# Monte carlo simulation
res <- MonteCarlo.estimates(coredata(dataset),nreps=10,train.size=100,test.size=10,'svm',
list(fr ~ .,type='eps-regression', kernel='sigmoid'),'error.func',relearn.step=50)
