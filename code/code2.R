library(zoo)
library(tseries)
library(TTR)
library(kknn)
library(randomForest)
library(e1071)
library(rpart)

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

# Mean absolute error (MAE)
mae <- function(preds, trues)
{
    error <- mean(abs((preds - trues)),na.rm = T)

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
    for (i in 1:(dimension))
    {
        dataset.tmp <- lag(data, -(i*tlag))
        dataset <- cbind(dataset, dataset.tmp)
    }

    names <- c("fr")
    for (i in 1:(length(colnames(dataset))-1))
    {
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
# Provided in KDD MapI class by Luis Torgo
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
    #cat('Going for sliding window testing from row ',start,'(train.sz=',train.size,',Relearn=',relearn.step,')\n')
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
        #cat('*')
        learner.pars$data <- orig.data[(test.pos-window.size):(test.pos-1),]
        model <- do.call(learner,learner.pars)
        preds <- c(preds,predict(model,orig.data[test.pos:min(n,test.pos+relearn.step-1),]))
        test.pos <- test.pos+relearn.step
    }
    #cat('\n')

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
    list(mae = mae(preds,test.data[,target.name(form,test.data)]),
    mape = mape(preds,test.data[,target.name(form,test.data)]),
    mse = mse(preds,test.data[,target.name(form,test.data)]),
    mad = mad(preds,test.data[,target.name(form,test.data)]))
}

###############################################################################
# Code
###############################################################################

# Load data from CSV
raw.data <- read.csv("data2.csv", header = TRUE, stringsAsFactors = FALSE)

# Convert raw data into a TS object
ts.data <- ts(data = raw.data[2:3], frequency = 4, start = c(2000,1), end = c(2011,2))

# Plot time series
plot.ts(ts.data)

# Pre-processing time series
ts.data <- pre.processing(ts.data)
#plot.ts(ts.data.processed)

# Original dataset
# Trimester
cat("Trimester\n")
dataset <- embedded.dataset(ts.data,1,1)

res.svm <- MonteCarlo.estimates(coredata(dataset),nreps=10,train.size=20,test.size=6,'svm',
list(fr ~ .,type='eps-regression', kernel='sigmoid'),'error.func',relearn.step=10)

res.dt <- MonteCarlo.estimates(coredata(dataset),nreps=10,train.size=20,test.size=6, 'rpart',
list(fr ~ .,method="anova"), 'error.func', relearn.step=10)

res.rf <- MonteCarlo.estimates(coredata(dataset),nreps=10,train.size=20,test.size=6,'randomForest',
list(fr ~ .),'error.func',relearn.step=10)

cat("DT Results\n")
print(colMeans(res.dt$Results))
cat("SVM Results\n")
print(colMeans(res.svm$Results))
cat("RF Results\n")
print(colMeans(res.rf$Results))
cat("\n")

# Semester
cat("Semester\n")
dataset <- embedded.dataset(ts.data,2,1)

res.svm <- MonteCarlo.estimates(coredata(dataset),nreps=10,train.size=20,test.size=6,'svm',
list(fr ~ .,type='eps-regression', kernel='sigmoid'),'error.func',relearn.step=10)

res.dt <- MonteCarlo.estimates(coredata(dataset),nreps=10,train.size=20,test.size=6, 'rpart',
list(fr ~ .,method="anova"), 'error.func', relearn.step=10)

res.rf <- MonteCarlo.estimates(coredata(dataset),nreps=10,train.size=20,test.size=6,'randomForest',
list(fr ~ .),'error.func',relearn.step=10)

cat("DT Results\n")
print(colMeans(res.dt$Results))
cat("SVM Results\n")
print(colMeans(res.svm$Results))
cat("RF Results\n")
print(colMeans(res.rf$Results))
cat("\n")

# Anual
cat("Anual\n")
dataset <- embedded.dataset(ts.data,4,1)

res.svm <- MonteCarlo.estimates(coredata(dataset),nreps=10,train.size=15,test.size=5,'svm',
list(fr ~ .,type='eps-regression', kernel='sigmoid'),'error.func',relearn.step=10)

res.dt <- MonteCarlo.estimates(coredata(dataset),nreps=10,train.size=15,test.size=5, 'rpart',
list(fr ~ .,method="anova"), 'error.func', relearn.step=10)

res.rf <- MonteCarlo.estimates(coredata(dataset),nreps=10,train.size=15,test.size=5,'randomForest',
list(fr ~ .),'error.func',relearn.step=10)

cat("DT Results\n")
print(colMeans(res.dt$Results))
cat("SVM Results\n")
print(colMeans(res.svm$Results))
cat("RF Results\n")
print(colMeans(res.rf$Results))
cat("\n\n")

# Homologous period
cat("Homologous period\n")

cat("Trimester\n")
dataset <- embedded.dataset(ts.data,1,4)

res.svm <- MonteCarlo.estimates(coredata(dataset),nreps=10,train.size=20,test.size=6,'svm',
list(fr ~ .,type='eps-regression', kernel='sigmoid'),'error.func',relearn.step=10)

res.dt <- MonteCarlo.estimates(coredata(dataset),nreps=10,train.size=20,test.size=6, 'rpart',
list(fr ~ .,method="anova"), 'error.func', relearn.step=10)

res.rf <- MonteCarlo.estimates(coredata(dataset),nreps=10,train.size=20,test.size=6,'randomForest',
list(fr ~ .),'error.func',relearn.step=10)

cat("DT Results\n")
print(colMeans(res.dt$Results))
cat("SVM Results\n")
print(colMeans(res.svm$Results))
cat("RF Results\n")
print(colMeans(res.rf$Results))
cat("\n")

# Semester
cat("Semester\n")
dataset <- embedded.dataset(ts.data,2,4)

res.svm <- MonteCarlo.estimates(coredata(dataset),nreps=10,train.size=20,test.size=6,'svm',
list(fr ~ .,type='eps-regression', kernel='sigmoid'),'error.func',relearn.step=10)

res.dt <- MonteCarlo.estimates(coredata(dataset),nreps=10,train.size=20,test.size=6, 'rpart',
list(fr ~ .,method="anova"), 'error.func', relearn.step=10)

res.rf <- MonteCarlo.estimates(coredata(dataset),nreps=10,train.size=20,test.size=6,'randomForest',
list(fr ~ .),'error.func',relearn.step=10)

cat("DT Results\n")
print(colMeans(res.dt$Results))
cat("SVM Results\n")
print(colMeans(res.svm$Results))
cat("RF Results\n")
print(colMeans(res.rf$Results))
cat("\n")

# Anual
cat("Anual\n")
dataset <- embedded.dataset(ts.data,4,4)

res.svm <- MonteCarlo.estimates(coredata(dataset),nreps=10,train.size=15,test.size=5,'svm',
list(fr ~ .,type='eps-regression', kernel='sigmoid'),'error.func',relearn.step=10)

res.dt <- MonteCarlo.estimates(coredata(dataset),nreps=10,train.size=15,test.size=5, 'rpart',
list(fr ~ .,method="anova"), 'error.func', relearn.step=10)

res.rf <- MonteCarlo.estimates(coredata(dataset),nreps=10,train.size=15,test.size=5,'randomForest',
list(fr ~ .),'error.func',relearn.step=10)

cat("DT Results\n")
print(colMeans(res.dt$Results))
cat("SVM Results\n")
print(colMeans(res.svm$Results))
cat("RF Results\n")
print(colMeans(res.rf$Results))
cat("\n\n")
