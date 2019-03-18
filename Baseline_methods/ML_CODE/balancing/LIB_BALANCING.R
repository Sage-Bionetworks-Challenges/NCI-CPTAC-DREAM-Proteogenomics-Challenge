# ------------------------------------------------------------------------------------------
# DESC: Balance a vector of response data, by splitting the data in n bins of equal size. 
#       This function resamples the data so that each bin has afterwards the size of the 
#       largest bin. Empty bins remain empty. In case the modulo of largest bin and bin to 
#       balance is larger than zero, randomly samples will be taken from the bin to balance
#       for equalling size with the largest bin. 
# IN:   res   =>  Vector containing the drug response data. This vector needs to have the 
#                 the COSMIC ids as name.
#       nBins =>  Number of splits (default = 10).
# OUT:  Returns a balanced response vector. 
# ------------------------------------------------------------------------------------------
balance <- function(res, nBins=10) {
  
  if (length(res) < 10) 
    stop("The response vector neesds at least 10 elements!")
  
  if (!is.numeric(res)) 
    stop("Response need to be a numeric vector!")
  
  if (is.null(names(res)))
    stop("Response vector needs names, which are the ids of each entry!")
  
  if (length(unique(names(res))) != length(names(res)))
    stop("The names of the response vector being not unique IDs!")
  
  membership <- as.numeric(cut(res, nBins))
  counts<-table(membership)
  
  # how often each bin will be resampled
  nResample_new <- floor(max(counts) / counts)
  
  # how many random samples taken from each bin
  nDraw_new <- max(counts) %% counts
  
  bRes <- c()
  for (member in names(counts)) {
    binRes <- res[membership == member]
    bRes <- c(bRes, rep(binRes, nResample_new[member]))
    # only if a random number of sample has to be drawn to fill up the bin
    if (nDraw_new[member]>0) {
      bRes <- c(bRes, sample(binRes, nDraw_new[member]))
    }
  }
  
  # hist is ignoring "breaks" parameter, therefore alternative visualization. 
  # plot(cut(res, nBins), las=2)
  # plot(cut(bRes, nBins), las=2)
  
  return(sample(bRes))
}