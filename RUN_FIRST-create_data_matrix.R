rm(list = ls())

# Install and load required packages
if ("lubridate" %in% rownames(installed.packages()) == F){
  install.packages("lubridate")
}
if ("RAMP" %in% rownames(installed.packages()) == F){
  install.packages("RAMP")
}
library(lubridate)
library(RAMP)

# Set base directory
base.dir <- 'https://raw.github.com/wsdaniels/COmodeling/main/'

# GET RESPONSE DATA
response <- read.csv(paste0(base.dir, "MSEA_V8JMOPITT_weeklyanomalies_WEDCEN_nofill.csv"))
response$time <- ymd(response$time)

# Placeholder for missing CO values. This will get used later
missing.val <- -9999


# GET PREDICTOR DATA
nino <- read.csv(paste0(base.dir, "nino34_weekly_avg.csv"))
aao  <- read.csv(paste0(base.dir, "aao_weekly_avg.csv"))
tsa  <- read.csv(paste0(base.dir, "tsa_weekly_avg.csv"))
dmi  <- read.csv(paste0(base.dir, "dmi_weekly_avg.csv"))
olr.msea  <- read.csv(paste0(base.dir,  "msea_olr.csv"))

# remove partial first entry in olr
olr.msea <- olr.msea[2:nrow(olr.msea),]

# PUT PREDICTOR DATA INTO A LIST
predictors <- list("nino" = nino,
                   "dmi" = dmi,
                   "tsa" = tsa,
                   "aao" = aao,
                   "olr.msea" = olr.msea)

# Fix column alignment
# NOTE: aao doesn't need the correction for some reason
for (i in 1:3){
  this.var <- predictors[[i]]
  this.var[,2] <- this.var[,1]
  this.var[,1] <- rownames(this.var)
  row.names(this.var) <- NULL
  predictors[[i]] <- this.var
}

# Convert time variable to lubridate datetime
for (i in 1:length(predictors)){
  this.var <- predictors[[i]]
  this.var$time <- ymd(this.var$time)
  predictors[[i]] <- this.var
}

# Clean up
rm(aao,dmi,nino,tsa,olr.msea,this.var,i)

# ALIGN START OF PREDICTOR TIME SERIES
# Get the latest start date - this will be used to align the starts
start.date <- max(as_date(sapply(predictors, function(X) X$time[1])))

for (i in 1:length(predictors)){
  this.var <- predictors[[i]]
  this.var <- this.var[!(this.var$time < start.date), ]
  predictors[[i]] <- this.var
}


# ALIGN END OF PREDICTOR TIME SERIES
terminal.date <- response$time[length(response$time)]

for (i in 1:length(predictors)){
  this.var <- predictors[[i]]
  this.var <- this.var[!(this.var$time > terminal.date), ]
  predictors[[i]] <- this.var
}

# Clean up
rm(this.var, terminal.date, i)

#------------------------------------------------
# ALL END DATES ARE ALLIGNED AT THIS POINT
# ALL TIME VARIABLES ARE THE SAME AT THIS POINT
#------------------------------------------------

# COMPUTE OFFSETS FROM FINAL RESPONSE OBSERVATION
# This will be used for lag calculations later
for (i in 1:length(predictors)){
  predictors[[i]]$offset <- seq(nrow(predictors[[i]])-1, 0)
}
response$offset <- seq(nrow(response)-1, 0)

# COMPUTE SMOOTHED CURVES
# Here we compute gradually smoother gaussian kernels
# Smoothed indices are used for longer lags
x.dist <- response$time[2] - response$time[1]
max.mult <- 8
min.mult <- 1
mult.seq <- seq(min.mult, max.mult, length.out = 8)

for (i in 1:length(predictors)){
  for (j in 1:length(mult.seq)){
    
    this.gaussian.kernel <- ksmooth(x = predictors[[i]]$time,
                                    y = predictors[[i]]$anomaly,
                                    kernel = "normal",
                                    bandwidth = mult.seq[j]*x.dist)
    
    
    predictors[[i]][ length(predictors[[i]]) + 1 ] <- this.gaussian.kernel$y
    names(predictors[[i]])[ length(predictors[[i]]) ] <- paste0("gaussian.kernel.",
                                                                mult.seq[j])
    
  }
}

rm(this.gaussian.kernel, i, j, x.dist)

# CREATE MASKS FOR THE MONTHS WE WANT TO EXPLAIN THE RESPONSE
#---------------------------------------------------------------------
# Period over which to explain response data (in months)
# This corresponds to fire season in MSEA
study.period <- c(9, 12)

# Create vector of months to keep
if (study.period[1] <= study.period[2]){
  months.to.keep <- study.period[1] : study.period[2]
} else {
  months.to.keep <- c(study.period[1] : 12, 1 : study.period[2])
}

# Create month mask
month.mask <- month(response$time) %in% months.to.keep

# Apply month mask
response <- response[month.mask,]

rm(month.mask, months.to.keep, study.period)
#---------------------------------------------------------------------

# Remove NA from response
# NOTE: this must be done AFTER setting up the offsets!!!!!!
response <- response[!(response$num_obs == missing.val), ]


# SET LAG LIMITS
min.lag <- 1
max.lag <- 52
lag.vals <- list(nino = min.lag:max.lag,
                 dmi = min.lag:max.lag,
                 tsa = min.lag:max.lag,
                 aao = min.lag:max.lag,
                 olr.msea = min.lag:max.lag)


# DEFINE SMOOTHING PARAMETERS
start.smoothing <- 4
lag.dividers <- seq(start.smoothing, max.lag, length.out = length(mult.seq)+1)
lag.dividers <- round(lag.dividers)

# COUNT TOTAL NUMBER OF LAGS TO CONSIDER
n.lag <- 0
for (i in 1:length(lag.vals)){
  n.lag <- n.lag + length(lag.vals[[i]])
}

# BUILD DATA MATRIX TO BE USED IN RAMP
data.matrix <- data.frame(matrix(NA, ncol = n.lag, nrow = nrow(response)))

# FILL DATA MATRIX
it <- 1
for (i in 1:length(predictors)){
  these.lags <- lag.vals[[i]]
  if (length(these.lags) > 0){
    for (j in 1:length(these.lags)){
      
      var.offsets <- predictors[[i]]$offset
      required.offsets <- response$offset + these.lags[j]
      to.keep <- var.offsets %in% required.offsets
      
      if (these.lags[j] < start.smoothing){
        this.var <- predictors[[i]]$anomaly[to.keep]
      } else {
        
        for (k in seq( 1, length(lag.dividers)-1 )){
          if (these.lags[j] >= lag.dividers[k] & these.lags[j] < lag.dividers[k+1]){
            this.string <- paste0("gaussian.kernel.", mult.seq[k])
            this.var <- predictors[[i]][, this.string][to.keep]
            break
          }
        }
        
      }
      data.matrix[,it] <- this.var
      var.name <- paste0(names(predictors)[i], "_", these.lags[j])
      names(data.matrix)[it] <- var.name
      it <- it + 1
    }
  }
}

rm(i, it, j, n.lag, required.offsets, these.lags, this.var, to.keep,
   var.name, var.offsets, lag.vals)


write.csv(data.matrix, "data_matrix.csv")
write.csv(response,    "response.csv")
