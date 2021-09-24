rm(list = ls())

library(lubridate)
library(RAMP)
library(foreach)
library(doParallel)

# Set parallel computing parameters
# run.in.parallel == F runs in serial
run.in.parallel <- T
num.cores.to.use <- 7

# Set base directory
base.dir <- 'https://raw.github.com/wsdaniels/COmodeling/main/'

# Read in data matrix and response variable
# See "create_data_matrix.R" for code to generate these files
data.matrix <- read.csv(paste0(base.dir, "data_matrix.csv"))
response <-    read.csv(paste0(base.dir, "response.csv"))

# Clean up data
data.matrix <- data.matrix[,2:ncol(data.matrix)]
response <- response[,2:ncol(response)]

# Prep data for use in RAMP algorithm
# Scale predictor variables
X <- scale(data.matrix)
y <- response$anomaly_co

# Set up RAMP parameter sequences (gamma and ebic)
lseq <- function(from, to, length.out) {
  exp(seq(log(from), log(to), length.out = length.out))
}
gamma.vals <- lseq(1.001, 18, length.out = 140)[1:88]
ebic.vals <-  round(1000*(1-lseq(1, 0.001, length.out = 16)[1:11]))/1000

# Set up matrices to hold resulting models and corresponding EBIC values
form.mat <- matrix(NA, nrow = length(gamma.vals), ncol = length(ebic.vals))
crit.mat <- matrix(NA, nrow = length(gamma.vals), ncol = length(ebic.vals))

# BEGIN GRID SEARCH OVER RAMP PARAMETERS
for (a in 1:length(ebic.vals)){
  
  start.time <- Sys.time()
  print(paste0("Working on ebic gamma ", a, " / ", length(ebic.vals)))
  
  if (run.in.parallel){
    cl <- makeCluster(num.cores.to.use)
    registerDoParallel(cl)
    output <- foreach(b = 1:length(gamma.vals)) %dopar% {
      
      ramp.fit <- RAMP::RAMP(X = X, y = y,
                             penalty = "MCP",
                             tune = "EBIC",
                             n.lambda = 500,
                             ebic.gamma = ebic.vals[a],
                             gamma = gamma.vals[b])   
      
      # CREATE LM FORMULA FROM RAMP OUTPUT
      #------------------------------------------------------------------
      main.terms <- ramp.fit$mainInd
      main.terms <- names(data.matrix)[main.terms]
      
      interactions <- ramp.fit$interInd
      if (!is.null(interactions)){
        for (i in 1:length(interactions)){
          
          this.term <- interactions[i]
          these.subterms <- as.integer(strsplit(this.term, "X")[[1]][2:3])
          these.subterms <- names(data.matrix)[these.subterms]
          
          if (these.subterms[1] != these.subterms[2]){
            this.term <- paste(these.subterms, collapse = ":")
          } else {
            this.term <- paste0("I(", these.subterms[1], "^2)")
          }
          
          interactions[i] <- this.term
        }
      }
      
      if (!is.null(interactions)){
        model.string <- paste( paste(main.terms, collapse = " + "),
                               paste(interactions, collapse = " + "),
                               sep = " + ")
      } else {
        model.string <- paste(main.terms, collapse = " + ")
      }
      
      model.string <- paste0("co ~ ", model.string)
      
      list(formula = model.string,
           crit = ramp.fit$cri.list[ramp.fit$cri.loc])
    }
    
    stopCluster(cl)
    
  } else {
    output <- vector(mode = "list", length = length(gamma.vals))
    
    for (b in 1:length(gamma.vals)){
      
      ramp.fit <- RAMP::RAMP(X = X, y = y,
                             penalty = "MCP",
                             tune = "EBIC",
                             n.lambda = 500,
                             ebic.gamma = ebic.vals[a],
                             gamma = gamma.vals[b])   
      
      # CREATE LM FORMULA FROM RAMP OUTPUT
      #------------------------------------------------------------------
      main.terms <- ramp.fit$mainInd
      main.terms <- names(data.matrix)[main.terms]
      
      interactions <- ramp.fit$interInd
      if (!is.null(interactions)){
        for (i in 1:length(interactions)){
          
          this.term <- interactions[i]
          these.subterms <- as.integer(strsplit(this.term, "X")[[1]][2:3])
          these.subterms <- names(data.matrix)[these.subterms]
          
          if (these.subterms[1] != these.subterms[2]){
            this.term <- paste(these.subterms, collapse = ":")
          } else {
            this.term <- paste0("I(", these.subterms[1], "^2)")
          }
          
          interactions[i] <- this.term
        }
      }
      
      if (!is.null(interactions)){
        model.string <- paste( paste(main.terms, collapse = " + "),
                               paste(interactions, collapse = " + "),
                               sep = " + ")
      } else {
        model.string <- paste(main.terms, collapse = " + ")
      }
      
      model.string <- paste0("co ~ ", model.string)
      
      output[[b]] <- list(formula = model.string,
                          crit = ramp.fit$cri.list[ramp.fit$cri.loc])
    }
  }
  
  for (b in 1:length(output)){
    form.mat[b,a] <- output[[b]]$formula
    crit.mat[b,a] <- output[[b]]$crit
  }
  
  rm(output)
  gc()
  
  time.diff <- Sys.time() - start.time
  print(paste0("previous ebic.val took this long: ", time.diff))
  
}


# Find index of best model for each EBIC gamma value
min.ind.vec <- apply(crit.mat, 2, which.min)

# List to hold best model for each EBIC gamma value
best.form.list <- vector(mode = "list", length = length(ebic.vals))

# Recreate lm data matrix
lm.data.matrix <- as.data.frame(cbind(y, X))
names(lm.data.matrix)[1] <- "co"

# Recreate best model for each EBIC gamma value
for (i in 1:length(ebic.vals)){
  this.best.form <- formula(form.mat[min.ind.vec[i], i])
  best.form.list[[i]] <- lm(this.best.form, data = lm.data.matrix)
}

# Create vector to booleans indicating where unique models are located
new.form <- vector(length = length(ebic.vals))
new.form[1] <- T

for (i in 2:length(ebic.vals)){
  this.formula <- formula(best.form.list[[i]])
  previous.formula <- formula(best.form.list[[i-1]])
  new.form[i] <- this.formula != previous.formula
}

# Vector holding only the unique models
unique.forms <- best.form.list[new.form]

write.csv(unique.forms, paste0("best_models.csv"))
write.csv(form.mat,     paste0("formula_matrix.csv"))
write.csv(crit.mat,     paste0("EBIC_gamma_matrix.csv"))
