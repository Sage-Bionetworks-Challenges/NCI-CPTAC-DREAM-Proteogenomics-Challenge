
##################################################################### Standard Elastic Net ########################################################################
###################################################################################################################################################################
library(glmnet)
trainEN <- function(train_feat, train_obs, xTrain_feat, xTrain_obs, alpha) {
  fit <- glmnet(train_feat, train_obs, alpha = alpha) ## betas and lambdas are computed from the training set
  
  xTrain_pred <- predict(fit,type="response",newx=xTrain_feat) ## rows are samples, columns are each step of determination of beta coefficients 
  train_pred <- predict(fit,type="response",newx=train_feat)
  
  xTrain_cIDX <- c()
 
  for (i in 1:ncol(fit$beta)) {
    xTrain_cIDX <- c(xTrain_cIDX, cIDX(xTrain_pred[,i], xTrain_obs))  ## Harrell c-index or concordance C
    ## fit$beta is a matrix of the beta weights across the different steps for determining the betas
    ## we save all the concordance index for every step of beta for training set and cross validation set
  }
  
  lambda <- fit$lambda[tail(which(xTrain_cIDX == max(xTrain_cIDX)), n=1)]
  ## choose the lambda for the highest c index in the cross validation set
  
  return(list(model=fit,
              lambda=lambda,
              xTrain_cIDX=xTrain_cIDX,
              featWeights=fit$beta[,lambda==fit$lambda]))
}

predictEN <- function(model, lambda, feat) {
  return(predict(model,type="response",newx=feat)[,model$lambda == lambda])
}

# i <- 1 ; fold <- 10  ;  iteration <- 1  ;  alpha <- 1 ; index <- i

run_EN <- function(index,RES, features, fold, iteration, alpha) {
  
  RES <- as.matrix(RES)
  features <- as.matrix(features)
  
  res <- RES[index, ]
  res <- res[!is.na(res)]
  if(length(res) < 30) { fold<-3 }
  
  pearson_table <- c()
  RMSE_table <- c()
  for (i in 1:iteration) {     ##  start iteration of k fold cross validation
    crossVal <- crossvalidation(names(res), nFold=fold, seed=F)
    observation <- c()
    prediction <- c()
    
    for (i in 1:fold) {
      
      train_obs <- res[as.character(crossVal[[i]]$train)]
      train_feat <- as.matrix(features[names(train_obs), ])
      rownames(train_feat) <- names(train_obs)
      train_feat[is.na(train_feat)] <- 0  ## replace NA by 0
      
      # prepare cross-trainset
      xTrain_obs <- res[as.character(crossVal[[i]]$xTrain)]
      xTrain_feat <- as.matrix(features[names(xTrain_obs), ])
      if (ncol(xTrain_feat) == 1 ) {
        xTrain_feat <- t(xTrain_feat)
      }
      rownames(xTrain_feat) <- names(xTrain_obs)
      xTrain_feat[is.na(xTrain_feat)] <- 0  ## replace NA by 0
      
      # train an elastic net model
      trainOut <- trainEN(train_feat, train_obs, xTrain_feat, xTrain_obs, alpha)
      str(trainOut)                       
      
      # extract test set
      test_obs <- res[as.character(crossVal[[i]]$test)]
      test_feat <- as.matrix(features[names(test_obs), ])
      rownames(test_feat) <- names(test_obs)
      test_feat[is.na(test_feat)] <- 0  ## replace NA by 0
      
      # predict with test set
      test_pred <- predictEN(trainOut$model, trainOut$lambda, test_feat) 
      
      # compute pearson correlation and RMSE
      library(Hmisc) ; library(hydroGOF)
      x <- matrix(c(test_obs, test_pred),nrow=length(test_obs),ncol=2) 
      x <- na.omit(x)
      if( var(x[,1]) ==0) { x[1,1] <- x[1,1] + 1E-100} 
      if( var(x[,2]) ==0) { x[1,2] <- x[1,2] + 1E-100}
      
      # different metrics
      if(length(x[,1]) < 3) {
        pcorr <- 0 # pearson correlation
        pval <- 1 # p value
        result <- matrix(c(pcorr, pval, length(test_obs)),nrow=1,ncol=3) 
        colnames(result) <- c("pcorr", "pval", "Nobs")
        pearson_table <- rbind (pearson_table, result)
      } else {
      
      temp <- cor.test(x[,1], x[,2])
      pcorr <- temp$estimate # pearson correlation
      pval <- temp$p.value # p value
      result <- matrix(c(pcorr, pval, length(test_obs)),nrow=1,ncol=3) 
      colnames(result) <- c("pcorr", "pval", "Nobs")
      pearson_table <- rbind (pearson_table, result)
      
      RMSE_table <- c(RMSE_table , rmse(x[,1], x[,2], na.rm=TRUE) )
      }
    }
  } ## end of iteration 
  
  result_table <- list(pearson_table, RMSE_table) ; names(result_table) <- c("pcorr","rmse")
  return(result_table)
}

# i <- 1 ; fold <- 3  ;  iteration <- 1  ;  alpha <- 1 ; index <- i

create_data_setting4 <- function(RES_subset) {
  features_col_subset <- features_col [ colnames(RES_subset) , ]
  features_row_subset <- features_row [ rownames(RES_subset) , ]
  
  res <- as.numeric(RES_subset)
  res_name <- c() 
  features <- c()
  for(col in 1:length(colnames(RES_subset))) {
    for(row in 1:length(rownames(RES_subset))) {
      res_name <- c(res_name, paste0(rownames(RES_subset)[row],"_",colnames(RES_subset)[col]))
      features <- rbind(features, as.numeric(c(features_col_subset[col, ],features_row_subset[row, ])) )
    }
  }
  names(res) <- res_name
  rownames(features) <- res_name
  colnames(features) <- c(colnames(features_col), colnames(features_row))
  return(list(res,features))
}

create_Macau_combo <- function(features) {
  feature_combo_names <- expand.grid(colnames(features_col),colnames(features_row))
  feature_combo_names <- paste0(feature_combo_names$Var1,"_",feature_combo_names$Var2)
  
  Macau_combo <- matrix(nrow=length(rownames(features)),ncol = length(feature_combo_names))
  for(sample in 1:length(rownames(features))) {
    features_col_value <- features[sample , 1:length(colnames(features_col)) ]
    features_row_value <- features[sample , length(colnames(features_col))+1:length(colnames(features_row))]
    Macau_combo[sample, ] <- mapply(`*`,features_col_value,rep(features_row_value,each=length(features_col_value))) 
  }

  colnames(Macau_combo) <- feature_combo_names
  rownames(Macau_combo) <- rownames(features)
  return(Macau_combo)
}

run_EN_setting4 <- function(RES,features_col,features_row, fold, iteration, alpha) {

  pearson_table <- c()
  RMSE_table <- c()
  for (i in 1:iteration) {     ##  start iteration of k fold cross validation
    crossVal_row <- crossvalidation(rownames(RES), nFold=fold, seed=F)
    crossVal_col <- crossvalidation(colnames(RES), nFold=fold, seed=F)
   
    observation <- c()
    prediction <- c()
    
    for (i in 1:1) {
      
      # prepare trainset
      mat <- RES[as.character(crossVal_row[[i]]$train),as.character(crossVal_col[[i]]$train)]
      L <- create_data_setting4(RES_subset = mat)
      train_obs <- L[[1]] ; train_feat <- L[[2]]
      train_obs <- train_obs[complete.cases(train_obs)]
      train_feat <- train_feat[ names(train_obs) , ]
      
      # prepare cross-trainset
      mat <- RES[as.character(crossVal_row[[i]]$xTrain),as.character(crossVal_col[[i]]$xTrain)]
      L <- create_data_setting4(RES_subset = mat)
      xTrain_obs <- L[[1]] ; xTrain_feat <- L[[2]]
      xTrain_obs <- xTrain_obs[complete.cases(xTrain_obs)]
      xTrain_feat <- xTrain_feat[ names(xTrain_obs) , ]
      
      # train an elastic net model
      trainOut <- trainEN(train_feat, train_obs, xTrain_feat, xTrain_obs, alpha)
      str(trainOut)                       
      
      # extract test set
      mat <- RES[as.character(crossVal_row[[i]]$test),as.character(crossVal_col[[i]]$test)]
      L <- create_data_setting4(RES_subset = mat)
      test_obs <- L[[1]] ; test_feat <- L[[2]]
      test_obs <- test_obs[complete.cases(test_obs)]
      test_feat <- test_feat[ names(test_obs) , ]
      
      # predict with test set
      test_pred <- predictEN(trainOut$model, trainOut$lambda, test_feat) 
      
      # compute pearson correlation and RMSE
      library(Hmisc) ; library(hydroGOF)
      x <- matrix(c(test_obs, test_pred),nrow=length(test_obs),ncol=2) 
      x <- na.omit(x)
      if( var(x[,1]) ==0) { x[1,1] <- x[1,1] + 1E-100} 
      if( var(x[,2]) ==0) { x[1,2] <- x[1,2] + 1E-100}
      
      # different metrics
      if(length(x[,1]) < 3) {
        pcorr <- 0 # pearson correlation
        pval <- 1 # p value
        result <- matrix(c(pcorr, pval, length(test_obs)),nrow=1,ncol=3) 
        colnames(result) <- c("pcorr", "pval", "Nobs")
        pearson_table <- rbind (pearson_table, result)
      } else {
        
        temp <- cor.test(x[,1], x[,2])
        pcorr <- temp$estimate # pearson correlation
        pval <- temp$p.value # p value
        result <- matrix(c(pcorr, pval, length(test_obs)),nrow=1,ncol=3) 
        colnames(result) <- c("pcorr", "pval", "Nobs")
        pearson_table <- rbind (pearson_table, result)
        
        RMSE_table <- c(RMSE_table , rmse(x[,1], x[,2], na.rm=TRUE) )
      }
    }
  } ## end of iteration 
  
  result_table <- list(pearson_table, RMSE_table) ; names(result_table) <- c("pcorr","rmse")
  return(result_table)
}


# i <- 1 ; fold <- 3  ;  iteration <- 1  ;  alpha <- 1 ; index <- rownames(RES)[i]

run_EN_TRANSFER <- function(index, fold, iteration, alpha) {
  res <- RES[index, ]
  res <- res[!is.na(res)]
  
  result_table <- c()
  for (i in 1:iteration) {     ##  start iteration of k fold cross validation. i <- 1
    crossVal <- crossvalidation(names(res), nFold=fold, seed=F)
    observation <- c()
    prediction <- c()
    
    for (i in 1:fold) {
      
      train_obs <- res[as.character(crossVal[[i]]$train)]
      train_feat <- as.matrix(features[names(train_obs), ])
      rownames(train_feat) <- names(train_obs)
      train_feat[is.na(train_feat)] <- 0  ## replace NA by 0
      
      # prepare cross-trainset
      xTrain_obs <- res[as.character(crossVal[[i]]$xTrain)]
      xTrain_feat <- as.matrix(features[names(xTrain_obs), ])
      if (ncol(xTrain_feat) == 1 ) {
        xTrain_feat <- t(xTrain_feat)
      }
      rownames(xTrain_feat) <- names(xTrain_obs)
      xTrain_feat[is.na(xTrain_feat)] <- 0  ## replace NA by 0
      
      # train an elastic net model
      trainOut <- trainEN(train_feat, train_obs, xTrain_feat, xTrain_obs, alpha)
      str(trainOut)     
      
      # transfer set
      test_obs <- RES_test[drug, ]
      test_feat <- features_test
      test_feat[is.na(test_feat)] <- 0  ## replace NA by 0
      
      # predict with test set
      test_pred <- predictEN(trainOut$model, trainOut$lambda, test_feat) 
      
      # compute pearson correlation
      library(Hmisc)
      x <- matrix(c(test_obs, test_pred),nrow=length(test_obs),ncol=2) 
      x <- na.omit(x)
      if( var(x[,1]) ==0) { x[1,1] <- x[1,1] + 0.0001}
      if( var(x[,2]) ==0) { x[1,2] <- x[1,2] + 0.0001}
      
      if(length(x[,1]) < 4) {
        pcorr <- NA # pearson correlation
        pval <- NA # p value
        result <- matrix(c(pcorr, pval, length(test_obs)),nrow=1,ncol=3) 
        colnames(result) <- c("pcorr", "pval", "Nobs")
        result_table <- rbind (result_table, result)
      } else {
        
        temp <- cor.test(x[,1], x[,2])
        pcorr <- temp$estimate # pearson correlation
        pval <- temp$p.value # p value
        result <- matrix(c(pcorr, pval, length(test_obs)),nrow=1,ncol=3) 
        colnames(result) <- c("pcorr", "pval", "Nobs")
        result_table <- rbind (result_table, result)
      }
    }
  } ## end of iteration of k fold cross validation
  return(result_table)
}

SAVE_MODEL <- function(index, fold, iteration, alpha ) {
  res <- RES[index, ]
  res <- res[!is.na(res)]
  if( (length(res)/fold) < 4) { fold<-3 }
 
  average_weight <- c()
  average_lambda <- c()
  for (i in 1:iteration) {     ##  start iteration of k fold cross validation. i <- 1
    crossVal <- stratified_crossvalidation(names(res), nFold=fold, seed=F)
    observation <- c()
    prediction <- c()
    
    for (i in 1:fold) {
      
      train_obs <- res[as.character(crossVal[[i]]$train)]
      train_feat <- as.matrix(features[names(train_obs), ])
      rownames(train_feat) <- names(train_obs)
      train_feat[is.na(train_feat)] <- 0  ## replace NA by 0
      
      # prepare cross-trainset
      xTrain_obs <- res[as.character(crossVal[[i]]$xTrain)]
      xTrain_feat <- as.matrix(features[names(xTrain_obs), ])
      if (ncol(xTrain_feat) == 1 ) {
        xTrain_feat <- t(xTrain_feat)
      }
      rownames(xTrain_feat) <- names(xTrain_obs)
      xTrain_feat[is.na(xTrain_feat)] <- 0  ## replace NA by 0
      
      # train an elastic net model
      trainOut <- trainEN(train_feat, train_obs, xTrain_feat, xTrain_obs, alpha)
      str(trainOut)     
      
      # save weight on all available data: FINAL MODEL
      common <- intersect(names(res), rownames(features)) ; res <- res[common] ; features <- features[common, ]

      
      # save parameters
      lambda <- trainOut$lambda
      average_lambda <- cbind(average_lambda ,  lambda)
      
    } ## end of k fold cross validation
  } ## end of iteration 
  average_lambda <- mean(average_lambda)

  fit <- glmnet(features, res, alpha = alpha, lambda=average_lambda)
#   featWeights <- fit$beta ; weight <- as(featWeights, "matrix")
#   average_weight <- cbind(average_weight ,  weight[ ,1] )
#   return(average_weight)
#   # predict with test set
#   test_pred <- predictEN(fit,fit$lambda,features) 
  return(fit)
}

##################################################################### Logistic regression #####################################################################
###############################################################################################################################################################

library(glmnet)
library(ROCR)
trainEN_logit <- function(train_feat, train_obs, xTrain_feat, xTrain_obs, alpha) {
 
  fit <- glmnet(train_feat, train_obs, alpha = alpha, family = "binomial") ## betas and lambdas are computed from the training set
  
  xTrain_pred <- predict(fit,type="response",newx=xTrain_feat) ## rows are samples, columns are each step of determination of beta coefficients 
  train_pred <- predict(fit,type="response",newx=train_feat)
  
  xTrain_AUC <- c()
  for (i in 1:ncol(fit$beta)) {
    
    pr <- prediction( xTrain_pred[ ,i] , xTrain_obs )
    prf <- performance(pr, measure = "tpr", x.measure = "fpr")
    #  plot(prf)
    auc <- performance(pr, measure = "auc")
    auc <- auc@y.values[[1]]
    
    xTrain_AUC <- c(xTrain_AUC , auc )
    ## fit$beta is a matrix of the beta weights across the different steps for determining the betas
    ## we save all the concordance index for every step of beta for training set and cross validation set
  }
  
  lambda <- fit$lambda[tail(which(xTrain_AUC == max(xTrain_AUC)), n=1)]
  ## choose the lambda for the highest AUC in the cross validation set
  
  return(list(model=fit,
              lambda=lambda,
              xTrain_AUC=xTrain_AUC,
              featWeights=fit$beta[,lambda==fit$lambda]))
}

predictEN_logit <- function(model, lambda, feat) {
  return(predict(model,type="response",newx=feat)[,model$lambda == lambda])
}

#  i <- 6 ; fold <- 3 ; iteration <- 1  ;  alpha <- 0  ;  index <- rownames(RES)[i]

run_EN_logit <- function(RES,index,features, fold, iteration, alpha) {
  res <- RES[index, ]
  res <- res[!is.na(res)]
  
  result_table <- c()
  for (i in 1:iteration) {     ##  start iteration of k fold cross validation
    crossVal <- stratified_crossvalidation(names(res), nFold=fold, seed=F)
    
    for (i in 1:fold) {
      
      train_obs <- res[as.character(crossVal[[i]]$train)]
      train_feat <- as.matrix(features[names(train_obs), ])
      rownames(train_feat) <- names(train_obs)
      train_feat[is.na(train_feat)] <- 0  ## replace NA by 0
      
      # prepare cross-trainset
      xTrain_obs <- res[as.character(crossVal[[i]]$xTrain)]
      xTrain_feat <- as.matrix(features[names(xTrain_obs), ])
      if (ncol(xTrain_feat) == 1 ) {
        xTrain_feat <- t(xTrain_feat)
      }
      rownames(xTrain_feat) <- names(xTrain_obs)
      xTrain_feat[is.na(xTrain_feat)] <- 0  ## replace NA by 0
      
      # train an elastic net model
      trainOut <- trainEN_logit (train_feat, train_obs, xTrain_feat, xTrain_obs, alpha)
      str(trainOut)           
      
      # extract test set
      test_obs <- res[as.character(crossVal[[i]]$test)]
      test_feat <- as.matrix(features[names(test_obs), ])
      rownames(test_feat) <- names(test_obs)
      test_feat[is.na(test_feat)] <- 0  ## replace NA by 0
      
      # predict with test set
      test_pred <- predictEN (trainOut$model, trainOut$lambda, test_feat) 
      
      # compute AUC
      pr <- prediction(test_pred, test_obs )
      prf <- performance(pr, measure = "tpr", x.measure = "fpr")
      auc <- performance(pr, measure = "auc")
      auc <- auc@y.values[[1]]
      result_table <- c(result_table , auc)
    } ## end of iteration of k fold cross validation
  } 
  result_mean <- mean(result_table, na.rm=T)
  result_sd <- sd(result_table, na.rm=T)
  L <- list(result_mean,result_sd) ; names(L) <- c("mean","sd")
  return(L)
}



##################################################################### features ranking ########################################################################
###############################################################################################################################################################

EN_feature_selection <- function (index, RES,  features, fold, iteration, alpha) {
  
  RES <- as.matrix(RES)
  features <- as.matrix(features)
  
  res <- RES[index, ]
  res <- res[!is.na(res)]
  target <- colnames(features)
  
  if( (length(res)/fold) < 4) { fold<-3 }
  
  features_rank <- c() ; weight_matrix <- c()
  for (i in 1:iteration) {   
    
    crossVal <- stratified_crossvalidation(names(res), nFold=fold, seed=F)
    
    for (i in 1:fold) {
      
      features <- features[,sample(ncol(features))]
      
      train_obs <- res[as.character(crossVal[[i]]$train)]
      train_feat <- as.matrix(features[names(train_obs), ])
      rownames(train_feat) <- names(train_obs)
      train_feat[is.na(train_feat)] <- 0  ## replace NA by 0
      
      # prepare cross-trainset
      xTrain_obs <- res[as.character(crossVal[[i]]$xTrain)]
      xTrain_feat <- as.matrix(features[names(xTrain_obs), ])
      if (ncol(xTrain_feat) == 1 ) {
        xTrain_feat <- t(xTrain_feat)
      }
      rownames(xTrain_feat) <- names(xTrain_obs)
      xTrain_feat[is.na(xTrain_feat)] <- 0  ## replace NA by 0
      
      # add test set to training set, since we don't need the test set
      test_obs <- res[as.character(crossVal[[i]]$test)]
      test_feat <- as.matrix(features[names(test_obs), ])
      train_feat <- rbind(train_feat, test_feat)
      train_obs <- c(train_obs, test_obs)
      
      # train an elastic net model
      trainOut <- trainEN(train_feat, train_obs, xTrain_feat, xTrain_obs, alpha)
      str(trainOut)                       
      
      # get the rank of the features
      B <- trainOut$model$beta
      #   B <- as(B, "matrix")
      B[B==0] <- NA
      weight <- rowMeans(B, na.rm = TRUE) 
      weight <- weight[match(target, names(weight))]
      weight_matrix <- cbind(weight_matrix, weight)
      
      varSeqNumber <- sort(rowSums(trainOut$model$beta == 0), index.return=TRUE)$ix
      varSeqName <- colnames(features)[varSeqNumber]
      varSeqName <- matrix(varSeqName) ; varSeqName <- cbind(varSeqName, 1:length(varSeqName[ ,1])) ; rownames(varSeqName) <- varSeqName[ ,1] ; varSeqName <- varSeqName[ ,-1]
      varSeqName <- varSeqName[match(target, names(varSeqName))]
      
      features_rank <- cbind(features_rank, varSeqName)
      
    }
  }

  storage.mode(features_rank) <- "numeric"
  fea_list_for_1 <- rowMeans(features_rank) ; names(fea_list_for_1) <- target 
  fea_list_for_1 <- fea_list_for_1 [ order(fea_list_for_1, decreasing = FALSE) ]  ; ## fea_list_for_1 [1:10]
  weight <- rowMeans(weight_matrix, na.rm = TRUE)  ; ## weight[1:10]
  weight <- weight[match(names(fea_list_for_1), names(weight))]
  fea_TABLE_for_1 <- cbind(fea_list_for_1, weight) ; colnames(fea_TABLE_for_1) <- c("average_apparition_rank", "average_weight" )
  return(fea_TABLE_for_1)
}


EN_feature_selection_weight <- function (RES,index,features, fold, iteration, alpha) {
  res <- RES[index, ]
  res <- res[!is.na(res)]
  target <- colnames(features)
  
  if( (length(res)/fold) < 4) { fold<-3 }
  
  features_weight <- c()
  for (i in 1:iteration) {   
    
    crossVal <- crossvalidation(names(res), nFold=fold, seed=F)
    
    for (i in 1:fold) {
      
      features <- features[,sample(ncol(features))]
      
      train_obs <- res[as.character(crossVal[[i]]$train)]
      train_feat <- as.matrix(features[names(train_obs), ])
      rownames(train_feat) <- names(train_obs)
      train_feat[is.na(train_feat)] <- 0  ## replace NA by 0
      
      # prepare cross-trainset
      xTrain_obs <- res[as.character(crossVal[[i]]$xTrain)]
      xTrain_feat <- as.matrix(features[names(xTrain_obs), ])
      if (ncol(xTrain_feat) == 1 ) {
        xTrain_feat <- t(xTrain_feat)
      }
      rownames(xTrain_feat) <- names(xTrain_obs)
      xTrain_feat[is.na(xTrain_feat)] <- 0  ## replace NA by 0
      
      # train an elastic net model
      trainOut <- trainEN(train_feat, train_obs, xTrain_feat, xTrain_obs, alpha)
      str(trainOut)                       
      
      # get the rank of the features
      B <- trainOut$model$beta
      # B_m <- as(B, "matrix")
      total_weight <- rowMeans(B)
      
      total_weight <- total_weight[match(target, names(total_weight))]
      
      features_weight <- cbind ( features_weight, total_weight)
      
    }
  }
  features_weight <- rowMeans(features_weight)
  features_weight <- features_weight[ order(-abs(features_weight)) ]
  return (features_weight)
}

EN_feature_selection_weight_LOGIT <- function (RES,index,features, fold, iteration, alpha) {
  res <- RES[index, ]
  res <- res[!is.na(res)]
  target <- colnames(features)
  
  if( (length(res)/fold) < 4) { fold<-3 }
  
  features_weight <- c()
  for (i in 1:iteration) {   
    
    crossVal <- crossvalidation(names(res), nFold=fold, seed=F)
    
    for (i in 1:fold) {
      
      features <- features[,sample(ncol(features))]
      
      train_obs <- res[as.character(crossVal[[i]]$train)]
      train_feat <- as.matrix(features[names(train_obs), ])
      rownames(train_feat) <- names(train_obs)
      train_feat[is.na(train_feat)] <- 0  ## replace NA by 0
      
      # prepare cross-trainset
      xTrain_obs <- res[as.character(crossVal[[i]]$xTrain)]
      xTrain_feat <- as.matrix(features[names(xTrain_obs), ])
      if (ncol(xTrain_feat) == 1 ) {
        xTrain_feat <- t(xTrain_feat)
      }
      rownames(xTrain_feat) <- names(xTrain_obs)
      xTrain_feat[is.na(xTrain_feat)] <- 0  ## replace NA by 0
      
      # train an elastic net model
      trainOut <- trainEN_logit(train_feat, train_obs, xTrain_feat, xTrain_obs, alpha)
      str(trainOut)                       
      
      # get the rank of the features
      B <- trainOut$model$beta
      # B_m <- as(B, "matrix")
      total_weight <- rowMeans(B)
      
      total_weight <- total_weight[match(target, names(total_weight))]
      
      features_weight <- cbind ( features_weight, total_weight)
      
    }
  }
  features_weight <- rowMeans(features_weight)
  features_weight <- features_weight[ order(-abs(features_weight)) ]
  return (features_weight)
}

EN_feature_selection_weight_LOGIT_APPEARANCE <- function (RES,index,features, fold, iteration, alpha) {
  res <- RES[index, ]
  res <- res[!is.na(res)]
  target <- colnames(features)
  
  if( (length(res)/fold) < 4) { fold<-3 }
  
  features_weight <- c()
  for (i in 1:iteration) {   
    
    crossVal <- crossvalidation(names(res), nFold=fold, seed=F)
    
    for (i in 1:fold) {
      
      features <- features[,sample(ncol(features))]
      
      train_obs <- res[as.character(crossVal[[i]]$train)]
      train_feat <- as.matrix(features[names(train_obs), ])
      rownames(train_feat) <- names(train_obs)
      train_feat[is.na(train_feat)] <- 0  ## replace NA by 0
      
      # prepare cross-trainset
      xTrain_obs <- res[as.character(crossVal[[i]]$xTrain)]
      xTrain_feat <- as.matrix(features[names(xTrain_obs), ])
      if (ncol(xTrain_feat) == 1 ) {
        xTrain_feat <- t(xTrain_feat)
      }
      rownames(xTrain_feat) <- names(xTrain_obs)
      xTrain_feat[is.na(xTrain_feat)] <- 0  ## replace NA by 0
      
      # train an elastic net model
      trainOut <- trainEN_logit(train_feat, train_obs, xTrain_feat, xTrain_obs, alpha)
      str(trainOut)                       
      
      # get the rank of the features
      B <- trainOut$model$beta
      B_m <- as(B, "matrix")
      total_weight <- B_m
      
      total_weight <- total_weight[match(target, rownames(total_weight)), ]
      
      features_weight <- cbind ( features_weight, total_weight)
      
    }
  }
  features_weight <- t(features_weight)
  return (features_weight)
}


run_EN_setting4_APPEARANCE <- function(RES,features_col,features_row, fold, iteration, alpha) {
  
  RES <- data.matrix(RES)
  features_col <- data.matrix(features_col)
  features_row <- data.matrix(features_row)
  
  features_weight <- c()
  for (i in 1:iteration) {     ##  start iteration of k fold cross validation
    crossVal_row <- crossvalidation(rownames(RES), nFold=fold, seed=F)
    crossVal_col <- crossvalidation(colnames(RES), nFold=fold, seed=F)
    
    observation <- c()
    prediction <- c()
    
    for (i in 1:1) {
      
      # prepare trainset
      mat <- RES[as.character(crossVal_row[[i]]$train),as.character(crossVal_col[[i]]$train)]
      L <- create_data_setting4(RES_subset = mat)
      train_obs <- L[[1]] ; train_feat <- L[[2]]
      train_obs <- train_obs[complete.cases(train_obs)]
      train_feat <- train_feat[ names(train_obs) , ]
      
      # prepare cross-trainset
      mat <- RES[as.character(crossVal_row[[i]]$xTrain),as.character(crossVal_col[[i]]$xTrain)]
      L <- create_data_setting4(RES_subset = mat)
      xTrain_obs <- L[[1]] ; xTrain_feat <- L[[2]]
      xTrain_obs <- xTrain_obs[complete.cases(xTrain_obs)]
      xTrain_feat <- xTrain_feat[ names(xTrain_obs) , ]
      
      # train an elastic net model
      trainOut <- trainEN(train_feat, train_obs, xTrain_feat, xTrain_obs, alpha)
      str(trainOut)    
      
      # get the rank of the features
      B <- trainOut$model$beta
      B_m <- as(B, "matrix")
      total_weight <- B_m
      
      features_weight <- cbind ( features_weight, total_weight)
    }
  } ## end of iteration 
  features_weight <- t(features_weight)
  return(features_weight)
}


