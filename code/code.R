library(zoo)
library(tseries)

# Load data from CSV
raw.data <- read.csv("data.csv", header = TRUE, stringsAsFactors = FALSE)

# Convert raw data into a TS object
ts.data <- ts(data = raw.data[2], frequency = 12, start = c(1983,1), end = c(2011,11))
