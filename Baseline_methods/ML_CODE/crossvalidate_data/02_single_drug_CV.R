#-----------------------------------------------------------------------------------------------------------
#DESC:    This function performs a SingleDrug "nFold" crossvalidation.
#IN:      entities  => vector representing the names of Celllines.
#                      The crossvalidation is done on this vector.
#         nFold     => numeric value in charge of which tpye of "nFold" crossvalidation will 
#                      be applied to "entities".
#OUT:     Returns a list "cv", which contains as many lists as the size of "nFold". 
#         Each element in "cv" is called "fold_XXX".
#         Each of these "fold_XXX" contains three sets: "test", "xTrain" and "train". 
#ERROR:   This function throws an ERROR if: 
#             - "entities" isn?t of the class "numeric"
#             - an ERROR occurs in the function "crossvalidation"
#-----------------------------------------------------------------------------------------------------------

singlecross <- function(entities,nFold){
  
  #TEST that entities is numeric
  if(class(entities)!="numeric") stop("'entities' needs to be numeric",call. =FALSE)
  
  #crossvalidating entities and testing the results
  cv<-crossvalidation(entities,nFold)
  
  return (cv)
}