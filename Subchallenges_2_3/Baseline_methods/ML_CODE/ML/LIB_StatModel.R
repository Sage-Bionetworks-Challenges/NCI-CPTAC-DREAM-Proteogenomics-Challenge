# #######################################################################################################################################
# RandomForest Function
# #######################################################################################################################################

#-----------------------------------------------------------------------------------------------------------------------------------------
#DESC:  The function "trainRF" trains a "randomForest" model with the data from "train_feat". 
#       This function will return a list containing the best model and a sublist with the c-indxes for the train-, xTrain-set.
#
#IN:    train_feat      ==> train values obtained from crossvalidating the observated feature data.
#                           the model will be trained with this values
#       train_obs       ==> observed drug_response for the values of train_feat
#       xTrain_feat     ==> xTrain values obtained from crossvalidating the observed feature data.
#                           necessary to figure out when to stop the training of the model with train_feat
#       xTrain_obs      ==> observed drug_response for the values of xTrain_feat
#       nTreeOffset     ==> the number of trees the first generated model will have
#       nTreeStep       ==> number of trees to grow(or to be added to an existing ensemble of trees)
#                           default value = 1
#       maxTrees        ==> maximum value for the value of nTrees
#                           if nTrees bigger than maxTrees the training will be stoppped
#                           default value = 1000
#
#OUT:   returns a list that contains the best model created and a sublist with the c-indexes for the train-,xTrain-set and the trees
#       grown by the best performing model
#-----------------------------------------------------------------------------------------------------------------------------------------

trainRF <-function(train_feat, train_obs, xTrain_feat, xTrain_obs, nTreeOffset=10, nTreeStep=10, maxTrees=1000) {
  
  # check if xTrain_feat has more than two elements and is of type "numeric"
  if(nrow(xTrain_feat)<2 || class(xTrain_feat)=="numeric"){
    stop("xTrain_feat contains less elements than 2. So it is not possible to calculate a Cindex value, that's why the
         code will be exited")
  }
  #######################################################################
  # check that there are no NA values in any committed parameter
  #######################################################################
  if(any(is.na(xTrain_feat))){
    stop("xTrain_feat contains NA values which are causing problems! This should never happen...")
  }
  
  if(any(is.na(xTrain_obs))){
    stop("xTrain_obs contains NA values which are causing problems! This should never happen...")
  }
  
  if(any(is.na(train_feat))){
    stop("train_feat contains NA values which are causing problems! This should never happen...")
  }
  
  if(any(is.na(train_obs))){
    stop("train_obs contains NA values which are causing problems! This should never happen...")
  }
  
  # check if train_obs is unique!!! if it is unique throw Exception
  if(length(unique(train_obs))==1){
    stop("train_obs is unique!!! This should never happen...")
  }
  
  # check if there's variance in train_feat, if not throw Exception
  if(length(unique((as.vector(train_feat)))) == 1){
    stop("No variance in train_feat, all values are unique!")
  }
  
  # HINT:
  # concordance index cannot be calculated in case the xTrain observations are all the same. In this 
  # particular case, no rank can be build, and therefore, no cross-training is possible. As a solution,
  # the number of trees will be set to the max tree number (default=1000). 
  if(length(unique(xTrain_obs))==1) {
    return(list(model=randomForest(train_feat, y=train_obs, ntree=nTreeOffset), 
                performance=NA,
                nTree_best_performer=nTreeOffset,
                nTrees=NA))
  } 
  else {
    # set the parameter in charge of stopping the training to false, so that the training can start
    stopTrain <- FALSE
    # create a model, train it with the values of train_feat and calculate the c-index
    rf_old <- randomForest(train_feat, y=train_obs, ntree=nTreeOffset)
    cIdx_train <- cIDX(predict(rf_old, train_feat), train_obs)[1]
    
    # feed the model with the values from xtrain and calculate the c-index
    cIdx_xTrain <- cIDX(predict(rf_old, xTrain_feat), xTrain_obs)[1]
    # set cIdx_xTrain_old  to the value of cidx_xTrain 
    cIdx_xTrain_old <- cIdx_xTrain
    
    nTree <- nTreeOffset
    
    # train as long as stopTrain is false and the numberr of trees in the model is smaller than maxTrees
    while (!stopTrain  & rf_old$ntree < maxTrees) {
      # add nTreeStep trees to the existing model
      rf_new <- grow(rf_old, nTreeStep)
      # feed the model with the values from xtrain and calculate the c-index, define result as cIdx_xTrain_new
      cIdx_xTrain_new <-cIDX(predict(rf_new, xTrain_feat), xTrain_obs)[1]
      # only if there was an improved of the c-index between the old one and the new one
      if (cIdx_xTrain_old > cIdx_xTrain_new) {
        # stop training
        stopTrain <- TRUE
      } else {
        # set actual model to the previous one(rf_old) and the actual c-index of it to the previous one as well
        rf_old <- rf_new
        cIdx_xTrain_old <- cIdx_xTrain_new
      }
      # save all calculated c-indexes for the xTrain-set and the train-set in one list each
      cIdx_xTrain <- c(cIdx_xTrain, cIdx_xTrain_new)
      cIdx_train <- c(cIdx_train, cIDX(predict(rf_new, train_feat), train_obs)[1])
      
      nTree <- c(nTree, rf_new$ntree)
    }
  }
  
  return(list(model=rf_old, 
              performance=list(train_cIdx=cIdx_train, xTrain_cIdx=cIdx_xTrain),
              nTree_best_performer=rf_old$ntree,
              nTrees=nTree))
}

# #######################################################################################################################################
# ElasticNet Function
# #######################################################################################################################################

#-----------------------------------------------------------------------------------------------------------------------------------------
#DESC:  The function "trainEN" trains a "elasticNet" model with the data from "train_feat". 
#       This function will return a list containing the best model, the lambda calculated for the best model, 
#       the cIndex for the train and xTrain set and the weights of the features.
#
#IN:    train_feat      ==> train values obtained from crossvalidating the observated feature data.
#                           the model will be trained with this values
#       train_obs       ==> observed drug_response for the values of train_feat
#       xTrain_feat     ==> xTrain values obtained from crossvalidating the observed feature data.
#                           neccessary to find the best value for lambda
#       xTrain_obs      ==> observed drug_response for the values of xTrain_feat
#
# OUT:  model           ==> Elastic net model
#       lambda          ==> Chosen lambda for further predictions (parameter which has been optimize)
#       xTrain_cIDX     ==> Concordance indices of the cross-train data
#       train_cIDX      ==> Concordance indices of the train data
#       featWeights     ==> weigths for each feature of the chosen model
#-----------------------------------------------------------------------------------------------------------------------------------------

trainEN <- function(train_feat, train_obs, xTrain_feat, xTrain_obs) {
  
  # check if xTrain_feat has more than two elements and is of type "numeric"
  if(nrow(xTrain_feat)<2 || class(xTrain_feat)=="numeric"){
    stop("xTrain_feat contains less elements than 2. So it is not possible to calculate a Cindex value, that's why the
         code will be exited")
  }
  #######################################################################
  # check that there are no NA values in any committed parameter
  #######################################################################
  if(any(is.na(xTrain_feat))){
    stop("xTrain_feat contains NA values which are causing problems! This should never happen...")
  }
  
  if(any(is.na(xTrain_obs))){
    stop("xTrain_obs contains NA values which are causing problems! This should never happen...")
  }
  
  if(any(is.na(train_feat))){
    stop("train_feat contains NA values which are causing problems! This should never happen...")
  }
  
  if(any(is.na(train_obs))){
    stop("train_obs contains NA values which are causing problems! This should never happen...")
  }
  
  # check if train_obs is unique!!! if it is unique throw Exception
  if(length(unique(train_obs))==1){
    stop("No variance in train_obs!!! This should never happen...")
  }
  
  # check if there's variance in train_feat, if not throw Exception
  if(length(unique((as.vector(train_feat)))) == 1){
    stop("No variance in train_feat, all values are unique!")
  }
  
  # HINT:
  # concordance index cannot be calculated in case the xTrain observations are all the same. In this 
  # particular case, no rank can be build, and therefore, no cross-training is possible. As a solution,
  # lambda is set to the first quantile of fit$lambda. 
  # alpha is used with default settings
  if(length(unique(xTrain_obs))==1) {
    
    fit <- glmnet(train_feat, train_obs)
    lambda <- fit$lambda[length(fit$lambda) - length(fit$lambda) / 4]                
    
    return(list(model=fit,
                lambda=lambda,
                xTrain_cIDX=NA,
                train_cIDX=NA,
                featWeights=fit$beta[,lambda==fit$lambda]))
  } else{
    # finding optimal EN mixing parameter (alpha) 
    # and optimal number of input features (lambda)
    xTrain_cIDX <- NA
    for (a in seq(1, 0, -0.1)) {
      # train an elasticNet model with the train_feat and train_obs
      current_fit <- glmnet(train_feat, train_obs, alpha=a)
      
      # predict with the xTrain_feat
      current_xTrain_pred <- predict(current_fit, type="response",newx=xTrain_feat)
      
      # find best lambda by calculating the xTrain_cidx
      current_xTrain_cIDX <- c()
      for (i in 1:ncol(current_fit$beta)) {
        current_xTrain_cIDX <- c(current_xTrain_cIDX, cIDX(current_xTrain_pred[,i], xTrain_obs))
      }
      
      # choose current model with this alpha if best fit
      if (length(xTrain_cIDX) == 1 || (max(xTrain_cIDX) < max(current_xTrain_cIDX))) {
        xTrain_cIDX <- current_xTrain_cIDX
        fit <- current_fit
        lambda <- current_fit$lambda[which(current_xTrain_cIDX == max(current_xTrain_cIDX))[1]]
      }
    }
    
    # train set performance across differnt lambda's
    train_pred <- predict(fit, type="response",newx=train_feat)
    train_cIDX <- c()
    for (i in 1:ncol(fit$beta)) {
      train_cIDX <- c(train_cIDX, cIDX(train_pred[,i], train_obs))
    }
  }
  
  return(list(model=fit,
              lambda=lambda,
              xTrain_cIDX=xTrain_cIDX,
              train_cIDX=train_cIDX,
              featWeights=fit$beta[,lambda==fit$lambda]))
}