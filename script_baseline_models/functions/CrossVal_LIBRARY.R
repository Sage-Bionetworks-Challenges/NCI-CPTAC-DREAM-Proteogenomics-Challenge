#-----------------------------------------------------------------------------------------------------------
#DESC:    This function performs a "nFold" crossvalidation on a given input called "entity". 
#         It is testing the created list "cv" and throws ERROR?s, if those are not passed(see
#         ERROR-section below).
#IN:      entity        => vector, which should be crossvalidated
#         nFold         => numeric value in charge of which tpye of "nFold" crossvalidation will 
#                          be applied to "entity"
#         seed          => this parameter decides if a seed will be set or not 
#OUT:     Returns the list "cv", which contains as many lists as the size of "nFold". 
#         Each element in "cv" is called "fold_XXX".
#         Each of these "fold_XXX" contains three sets: "test", "xTrain" and "train"   
#ERROR:   This function throws an ERROR if:
#           - "entity" is not unique
#           - in the sets of every "fold_XXX" not only the values of "entity" exist
#           - there is a difference between the length of the union of the sets for every "fold_XX"
#             and the length of "entity"
#           - in the sets of one fold the values are the same
#-----------------------------------------------------------------------------------------------------------

crossvalidation<-function(entity, nFold, seed) {

# set seed
if(seed){
  set.seed(754)
}

#creating a shuffled vector from the values of "entity"
vec <- entity[sample(length(entity))]

#defining the size of each set
size <- floor(length(vec)/nFold)
binSize <- rep(size, nFold)
modulo <- length(vec) %% size
if (modulo > 0) {
  binSize[1:modulo] <- binSize[1:modulo] +1
}

#setting the start indicies and the end indicies for the sets
binEndIdx <- cumsum(binSize)
binStardIdx <- c(1, binEndIdx[1:(nFold-1)]+1)

#iterating over all the values of "vec" and putting them into a set
#creating the list "cv" with the sublists and sets in them
cv <- list()
for (i in 1:nFold){
  
  testIdx <- binStardIdx[i]:binEndIdx[i]
  test <- vec[testIdx]
  
  xTrainIdx <- binStardIdx[(i %% nFold) + 1]:binEndIdx[(i %% nFold) + 1]
  xTrain <- vec[xTrainIdx]
  
  train <- vec[setdiff(1:length(vec), union(testIdx, xTrainIdx))]
  
  cv[[paste("fold_", i, sep="")]] <- list(test=test, xTrain=xTrain, train=train)
}


#TEST cases for cross-validation

#TEST if entity is unique otherwise throws ERROR
if(length(unique(entity))!=length(entity)) stop("'entity' is NOT unique",call. =FALSE)
  
for (i in 1:nFold){
  
  #TEST that there is no difference between the length of the union of the sets in one fold and the length of "entity"
  if(length(union(c(cv[[i]]$test, cv[[i]]$xTrain, cv[[i]]$train), entity)) != length(entity)){ 
    stop("difference between the length of union of the sets in fold_",print(i),"and the length of 'entity'",call. =FALSE)
  }
  #TEST if the elements of train, xTrain and test set all exist in entity(so no added or removed)  
  if(length(c(cv[[i]]$test, cv[[i]]$xTrain, cv[[i]]$train))!= length(entity)){
    stop("in the sets of every 'fold_XXX' don't only exist  the values of 'entity'",call. =FALSE)
  }
  
  #TEST that there are no equal elements in xTrain and in test, within one fold
  if (length(intersect(cv[[i]]$xTrain, cv[[i]]$test)) > 0){
    stop("equal elements in the sets 'xTrain' and 'test' of the fold_",print(i),call. =FALSE)
  }
  
  #TEST that there are no equal elements in xTrain and in train, within one fold
  if (length(intersect(cv[[i]]$xTrain, cv[[i]]$train)) > 0){
    stop("equal elements in the sets 'xTrain' and 'train' of the fold_",print(i),call. =FALSE)
  }
  
  #TEST that there are no equal elements in test and in train, within one fold
  if (length(intersect(cv[[i]]$test, cv[[i]]$train)) > 0){
    stop("equal elements in the sets 'test' and 'train' of the fold_",print(i),call. =FALSE)
  }
 }
 return(cv) 
}

############################################# STRATIFIED CROSS VALIDATION #############################################
#######################################################################################################################
library(caret)

stratified_crossvalidation<-function(entity, nFold, seed) {
  
  # set seed
  if(seed){
    set.seed(754)
  }
  
  folds <- createFolds(factor(entity), k = nFold, list = FALSE)
   
  vec <- c()
  binStardIdx <- c()
  binEndIdx <- c()
  for( i in 1:nFold) {
    vec <- c(vec , entity[folds == i])
    if(i ==1) {
    binStardIdx <- i
    binEndIdx <- length(entity[folds == i]) } else {
      binStardIdx <- c(binStardIdx , binEndIdx[i-1]+1 )
      binEndIdx <- c(binEndIdx , binEndIdx[i-1]+length(entity[folds == i]) )
   }
  }

  #iterating over all the values of "vec" and putting them into a set
  #creating the list "cv" with the sublists and sets in them
  cv <- list()
  for (i in 1:nFold){
    
    testIdx <- binStardIdx[i]:binEndIdx[i]
    test <- vec[testIdx]
    
    xTrainIdx <- binStardIdx[(i %% nFold) + 1]:binEndIdx[(i %% nFold) + 1]
    xTrain <- vec[xTrainIdx]
    
    train <- vec[setdiff(1:length(vec), union(testIdx, xTrainIdx))]
    
    cv[[paste("fold_", i, sep="")]] <- list(test=test, xTrain=xTrain, train=train)
  }
  
  
  #TEST cases for cross-validation
  
  #TEST if entity is unique otherwise throws ERROR
  if(length(unique(entity))!=length(entity)) stop("'entity' is NOT unique",call. =FALSE)
  
  for (i in 1:nFold){
    
    #TEST that there is no difference between the length of the union of the sets in one fold and the length of "entity"
    if(length(union(c(cv[[i]]$test, cv[[i]]$xTrain, cv[[i]]$train), entity)) != length(entity)){ 
      stop("difference between the length of union of the sets in fold_",print(i),"and the length of 'entity'",call. =FALSE)
    }
    #TEST if the elements of train, xTrain and test set all exist in entity(so no added or removed)  
    if(length(c(cv[[i]]$test, cv[[i]]$xTrain, cv[[i]]$train))!= length(entity)){
      stop("in the sets of every 'fold_XXX' don't only exist  the values of 'entity'",call. =FALSE)
    }
    
    #TEST that there are no equal elements in xTrain and in test, within one fold
    if (length(intersect(cv[[i]]$xTrain, cv[[i]]$test)) > 0){
      stop("equal elements in the sets 'xTrain' and 'test' of the fold_",print(i),call. =FALSE)
    }
    
    #TEST that there are no equal elements in xTrain and in train, within one fold
    if (length(intersect(cv[[i]]$xTrain, cv[[i]]$train)) > 0){
      stop("equal elements in the sets 'xTrain' and 'train' of the fold_",print(i),call. =FALSE)
    }
    
    #TEST that there are no equal elements in test and in train, within one fold
    if (length(intersect(cv[[i]]$test, cv[[i]]$train)) > 0){
      stop("equal elements in the sets 'test' and 'train' of the fold_",print(i),call. =FALSE)
    }
  }
  return(cv) 
}
