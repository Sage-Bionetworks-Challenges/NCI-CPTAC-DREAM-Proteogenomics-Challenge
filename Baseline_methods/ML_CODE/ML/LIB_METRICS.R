library("Hmisc")

# obs and pred need to be numeric for all functions!

# function to calculate the concordance index
cIDX <- function(pred, obs) {
  # calculate the concordance index 
  rcorr.cens(pred, obs)[1]
}

# function to calculate the root mean square error
# observed (independent variable, x)
# predicted (dependent variable, y)
RMSE <- function(obs, pred) {
  sqrt(mean((obs - pred)^2, na.rm = TRUE)) 
}

# coefficient of determination / R squared
R2 <- function(obs, pred) {
  # total variation in the observed data
  SSyy = sum((obs-mean(obs))^2)
  # sqared error of prediction
  SSE = sum((obs-pred)^2)
  1 - SSE / SSyy
}

# The proportion of variance explained (PVE)
PVE <- function(obs, pred) {
  # covariance sum of squares
  SSxy <- sum((obs-mean(obs)) * (pred-mean(pred)))
  # sum of squares for variable X
  SSx  <- sum((obs-mean(obs))^2)
  # sum of squares for variable Y
  SSy  <- sum((pred-mean(pred))^2)
  
  numerator <- SSxy^2
  denominator <- SSx*SSy
  if (denominator > 0)
    return(numerator/ denominator)
  else
    return(0)
}