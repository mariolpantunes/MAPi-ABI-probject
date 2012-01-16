library(zoo)
library(tseries)
library(forecast)

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
    rows <- (dimension*tlag):(nrow(data))
    dataset <- dataset[rows,]

    return(dataset)
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
cat("MAPE:",mape(preds, ts.data.processed),"\n")
cat("MAD :",mad(preds, ts.data.processed),"\n")

# Generate simple dataset
simple.dataset <- embedded.dataset(ts.data)


