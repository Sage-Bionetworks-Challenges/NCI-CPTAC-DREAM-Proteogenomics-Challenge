# fold <- 10  ;  i <- 1  ;  clinvar = 'DAYS_TO_LAST_FOLLOWUP'  ;  alpha <- 0.5
trainSurvival <- function(CAT_train, stime_train, status_train, CAT_xTrain, stime_xTrain, status_xTrain, alpha) {
  library("glmnet")
  library('survcomp')
  fit <- glmnet( x=CAT_train , Surv(stime_train, status_train), family = "cox", alpha = alpha ) ## betas and lambdas are computed from the training set
  xTrain_pred = predict ( fit , newx = CAT_xTrain  ) ## rows are samples, columns are each step of determination of beta coefficients 
  
  #  data = list(x = CAT_train, time = stime_train, status = status_train)
  #  debug(betterPathCalc)  debug(SGL)
  #  fit_SGL <- SGL(data, index , type = "cox", alpha=alpha, nlam = 20)
  
  xTrain_cIDX <- c()
  for (i in 1:ncol(fit$beta)) {
    
    # Determine concordance
    xTrain_pred = predict ( fit , newx = CAT_xTrain, s = fit$lambda[i] )
    cindex_validation = concordance.index(xTrain_pred, surv.time = stime_xTrain, surv.event=status_xTrain, method = "noether") # package survcomp
    xTrain_cIDX <- c(xTrain_cIDX, cindex_validation$c.index) ; 
  }
  
  ## choose the lambda for the highest c index in the cross validation set
  xTrain_cIDX[is.na(xTrain_cIDX)] <- 0
  lambda <- fit$lambda[tail(which(xTrain_cIDX == max(xTrain_cIDX)), n=1)]
  
  return(list(model=fit, lambda=lambda) ) 
}


RUN_survival <- function(CAT, CLINICS, fold, iteration, alpha, event=event,clinvar=clinvar) {
  result_AUC <- c()
  for (i in 1:iteration) {     ##  start iteration of k fold cross validation
    crossVal <- stratified_crossvalidation(colnames(CAT), nFold=fold, seed=F)
    
    for (i in 1:fold) {
      
      # prepare trainset
      samples2study_train <- as.character(crossVal[[i]]$train) 
      CAT_train <- CAT [ , samples2study_train  ] 
      status_train = (CLINICS[ samples2study_train , event] != 0) + 0
      stime_train = as.numeric(CLINICS[ samples2study_train , clinvar ])
      
      # prepare cross-trainset
      samples2study_xTrain <- as.character(crossVal[[i]]$xTrain) 
      CAT_xTrain <- CAT [ , samples2study_xTrain  ] 
      status_xTrain = (CLINICS[ samples2study_xTrain , event] != 0) + 0
      stime_xTrain = as.numeric(CLINICS[ samples2study_xTrain , clinvar ])
      
      # train a survival model
      CAT_train <- as.matrix(CAT_train) ;  CAT_train <- t(CAT_train) ; 
      CAT_xTrain <- as.matrix(CAT_xTrain) ;  CAT_xTrain <- t(CAT_xTrain) ; 
      CAT_train[is.na(CAT_train)] <- 0  ## replace NA by 0
      CAT_xTrain[is.na(CAT_xTrain)] <- 0  ## replace NA by 0
      
      trainOut <- trainSurvival (CAT_train, stime_train, status_train, CAT_xTrain, stime_xTrain, status_xTrain, alpha)
      fit <- trainOut$model
      lambda <- trainOut$lambda
      
      # extract test set
      samples2study_test <- as.character(crossVal[[i]]$test)
      CAT_test <- CAT [ , samples2study_test] ; CAT_test <- t(CAT_test) ; 
      status_test = (CLINICS[ samples2study_test , event] != 0) + 0
      stime_test = as.numeric(CLINICS[ samples2study_test , clinvar ])
      
      # Create survival estimates on test data
      test_pred = predict(fit , newx = CAT_test , s = lambda)
      
      library(pROC)
      AUC <- auc(status_test,test_pred[ ,1]   ) 
      
#       # Determine AUC
#       tt <- unique(sort(stime_test[status_test == 1]))
#       AUC_store <- c()
#       for(i in 1:length(tt)) {
#         e <- tdrocc(test_pred, surv.time = stime_test, surv.event=status_test, time=tt[i] )
#         AUC_store <- c(AUC_store, e$AUC)
#       }
      result_AUC <- c( result_AUC , AUC ) 
    }
    
  } ## end of iteration of k fold cross validation
  return(result_AUC)
}


run_MOFA <- function(mat_list,numFactors,MOFA_iteration,dataDir) {
  library(MOFAtools)
  
  dataDir <- paste0(dataDir,"/",sample(1:100000, 1) ) 
  dir.create(file.path(dataDir) , recursive = T)
  
  MOFAobject <- createMOFAobject(mat_list)
  DirOptions <- list(
    "dataDir" = dataDir, # Folder to store the input matrices as .txt files, it can be a simple temporary folder
    "outFile" = paste0(dataDir,"/result.hdf5") # Output file of the model (use hdf5 extension)
  )
  ModelOptions <- getDefaultModelOpts(MOFAobject)
  ModelOptions$numFactors <- numFactors
  ModelOptions$likelihood <- likelihood
  TrainOptions <- getDefaultTrainOpts()
  TrainOptions$maxiter <- MOFA_iteration
  TrainOptions$DropFactorThreshold <- 0
  TrainOptions
  MOFAobject <- prepareMOFA(MOFAobject, 
                            DirOptions = DirOptions,
                            ModelOptions = ModelOptions,
                            TrainOptions = TrainOptions
                           )
  
  MOFAobject <- runMOFA(MOFAobject, DirOptions, mofaPath = "/home/my871390/.local/bin/mofa" )  ## "/home/my871390/.local/bin/mofa" "/Users/miyang/anaconda/bin/mofa"
  Factor_matrix <- MOFAobject@Expectations$Z$E[ ,-1]
  colnames(Factor_matrix) <- paste0("Factor_",1:length(colnames(Factor_matrix)))
  
  W_1 <- MOFAobject@Expectations$SW$predicted_protein$E[ ,-1] ; rownames(W_1) <- paste0("W_1_",rownames(W_1))
  W_2 <- MOFAobject@Expectations$SW$mRNA$E[ ,-1] ; rownames(W_2) <- paste0("W_2_",rownames(W_2))
  W_3 <- MOFAobject@Expectations$SW$mutation$E[ ,-1] ; rownames(W_3) <- paste0("W_3_",rownames(W_3))
  W <- rbind(W_1,W_2,W_3) ; colnames(W) <- colnames(Factor_matrix)
  
  unlink(dataDir, recursive = T)
  return(list(Factor_matrix, W)) 
}

matrix.sort <- function(matrix) {
  if (nrow(matrix) != ncol(matrix)) stop("Not diagonal")
  if(is.null(rownames(matrix))) rownames(matrix) <- 1:nrow(matrix)
  row.max <- apply(matrix,1,which.max)
  if(all(table(row.max) != 1)) stop("Ties cannot be resolved")
  matrix[names(sort(row.max)),]
}

RUN_survival_MOFA_1 <- function(CAT, CLINICS, fold, iteration,numFactors, MOFA_iteration,MOFA_dir, alpha, event=event,clinvar=clinvar) {
  result_AUC <- c()
  for (i in 1:iteration) {     ##  start iteration of k fold cross validation
    crossVal <- stratified_crossvalidation(colnames( CAT[[1]] ), nFold=fold, seed=F)
    
    for (i in 1:fold) {
      
      # prepare trainset
      samples2study_train <- as.character(crossVal[[i]]$train) 
      CAT_train <- lapply(CAT, `[`, samples2study_train)
      status_train = (CLINICS[ samples2study_train , event] != 0) + 0
      stime_train = as.numeric(CLINICS[ samples2study_train , clinvar ])
      
      # prepare cross-trainset
      samples2study_xTrain <- as.character(crossVal[[i]]$xTrain) 
      CAT_xTrain <- lapply(CAT, `[`, samples2study_xTrain)
      status_xTrain = (CLINICS[ samples2study_xTrain , event] != 0) + 0
      stime_xTrain = as.numeric(CLINICS[ samples2study_xTrain , clinvar ])
      
      # extract test set
      samples2study_test <- as.character(crossVal[[i]]$test)
      CAT_test <- lapply(CAT, `[`, samples2study_test)
      status_test = (CLINICS[ samples2study_test , event] != 0) + 0
      stime_test = as.numeric(CLINICS[ samples2study_test , clinvar ])
      
      # APPLY MOFA
      L  <- run_MOFA(CAT_train,numFactors,MOFA_iteration,dataDir=MOFA_dir)
      CAT_train <- L[[1]] ; W_train <- L[[2]]
      L  <- run_MOFA(CAT_xTrain,numFactors,MOFA_iteration,dataDir=MOFA_dir)
      CAT_xTrain <- L[[1]] ; W_xTrain <- L[[2]]
      L  <- run_MOFA(CAT_test,numFactors,MOFA_iteration,dataDir=MOFA_dir)
      CAT_test <- L[[1]] ; W_test <- L[[2]]
      
      common <- Reduce(intersect, list(rownames(W_train),rownames(W_xTrain),rownames(W_test)))
      W_train <- W_train[common, ]
      W_xTrain <- W_xTrain[common, ]
      W_test <- W_test[common, ]
      
      # Find and reorder the Factors with respect to the training set (as reference)
      # THIS IS NOT CORRECT AS IT DOESN'T TAKE THE SIGN INTO ACCOUNT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      correlation <- cor(W_train,W_xTrain) 
      correlation_sorted <- t(matrix.sort(t(correlation))) 
      CAT_xTrain <- CAT_xTrain[ ,colnames(correlation_sorted) ] 
      colnames(CAT_xTrain) <- colnames(CAT_train)
      
      correlation <- cor(W_train,W_test)
      correlation_sorted <- t(matrix.sort(t(correlation)))
      CAT_test <- CAT_test[ ,colnames(correlation_sorted) ] 
      colnames(CAT_test) <- colnames(CAT_train)
      
      unloadNamespace("MOFAtools")
      trainOut <- trainSurvival (CAT_train, stime_train, status_train, CAT_xTrain, stime_xTrain, status_xTrain, alpha)
      fit <- trainOut$model
      lambda <- trainOut$lambda
      
      # Create survival estimates on test data
      test_pred = predict(fit , newx = CAT_test , s = lambda)
      
      library(pROC)
      AUC <- auc(status_test,test_pred[ ,1]   ) 
      result_AUC <- c( result_AUC , AUC ) 
    }
    
  } ## end of iteration of k fold cross validation
  return(result_AUC)
}

RUN_survival_MOFA_2 <- function(CAT, CLINICS, fold, iteration,numFactors, MOFA_iteration,MOFA_dir, alpha, event=event,clinvar=clinvar) {
  result_AUC <- c()
  for (i in 1:iteration) {     ##  start iteration of k fold cross validation
    crossVal <- stratified_crossvalidation(colnames( CAT[[1]] ), nFold=fold, seed=F)
    
    for (i in 1:fold) {
      
      # prepare trainset
      samples2study_train <- as.character(crossVal[[i]]$train) 
      CAT_train <- lapply(CAT, `[`, samples2study_train)
      status_train = (CLINICS[ samples2study_train , event] != 0) + 0
      stime_train = as.numeric(CLINICS[ samples2study_train , clinvar ])
      
      # prepare cross-trainset
      samples2study_xTrain <- as.character(crossVal[[i]]$xTrain) 
      CAT_xTrain <- lapply(CAT, `[`, samples2study_xTrain)
      status_xTrain = (CLINICS[ samples2study_xTrain , event] != 0) + 0
      stime_xTrain = as.numeric(CLINICS[ samples2study_xTrain , clinvar ])
      
      # extract test set
      samples2study_test <- as.character(crossVal[[i]]$test)
      CAT_test <- lapply(CAT, `[`, samples2study_test)
      status_test = (CLINICS[ samples2study_test , event] != 0) + 0
      stime_test = as.numeric(CLINICS[ samples2study_test , clinvar ])
      
      # APPLY MOFA
      L  <- run_MOFA(CAT_train,numFactors,MOFA_iteration,dataDir=MOFA_dir)
      CAT_train <- L[[1]] ; W_train <- L[[2]]

      
      unloadNamespace("MOFAtools")
      trainOut <- trainSurvival (CAT_train, stime_train, status_train, CAT_xTrain, stime_xTrain, status_xTrain, alpha)
      fit <- trainOut$model
      lambda <- trainOut$lambda
      
      # Create survival estimates on test data
      test_pred = predict(fit , newx = CAT_test , s = lambda)
      
      library(pROC)
      AUC <- auc(status_test,test_pred[ ,1]   ) 
      result_AUC <- c( result_AUC , AUC ) 
    }
    
  } ## end of iteration of k fold cross validation
  return(result_AUC)
}

cox_survival_weight <- function(CAT, CLINICS, fold, iteration, alpha,no_event=0, clinvar='DAYS_TO_LAST_FOLLOWUP') {

  target <- rownames(CAT)
  
  features_weight <- c()
  for (i in 1:iteration) {     ##  start iteration of k fold cross validation
    crossVal <- crossvalidation(colnames(CAT), nFold=fold, seed=F)
    
    for (i in 1:fold) {
      
      CAT <- CAT[ sample(nrow(CAT)) , ]
      
      samples2study_train <- as.character(crossVal[[i]]$train) 
      CAT_train <- CAT [ , samples2study_train  ] 
      status_train = (CLINICS[ samples2study_train , ]$EVENT != no_event) + 0
      stime_train = as.numeric(CLINICS[ samples2study_train , clinvar ])
      
      # prepare cross-trainset
      samples2study_xTrain <- as.character(crossVal[[i]]$xTrain) 
      CAT_xTrain <- CAT [ , samples2study_xTrain  ] 
      status_xTrain = (CLINICS[ samples2study_xTrain , ]$EVENT != no_event) + 0
      stime_xTrain = as.numeric(CLINICS[ samples2study_xTrain , clinvar ])
      
      # train a survival model
      CAT_train <- as.matrix(CAT_train) ;  CAT_train <- t(CAT_train) ; 
      CAT_xTrain <- as.matrix(CAT_xTrain) ;  CAT_xTrain <- t(CAT_xTrain) ; 
      
      trainOut <- trainSurvival (CAT_train, stime_train, status_train, CAT_xTrain, stime_xTrain, status_xTrain, alpha)
     
      # get the rank of the features
      B <- trainOut$model$beta
      # B_m <- as(B, "matrix")
      total_weight <- rowSums(B)
      
      total_weight <- total_weight[match(target, names(total_weight))]
      
      features_weight <- cbind ( features_weight, total_weight)
      
    }
    
  } ## end of iteration of k fold cross validation
  features_weight <- rowMeans(features_weight)
  features_weight <- features_weight[ order(-abs(features_weight)) ]
  return (features_weight)
}


survival_analysis_continuous = function(CAT, CLINICS, samples2study, main = '', clinvar = 'DAYS_TO_LAST_FOLLOWUP' ){
  status = (CLINICS[ samples2study , ]$VITAL_STATUS != 'Alive') + 0
  stime = as.numeric(CLINICS[ samples2study , clinvar ])
  age = as.numeric(CLINICS[ samples2study , ]$AGE)
  gender = CLINICS[ samples2study , ]$GENDER
  RESULTS = list()
  n.dead.cases = table(status)['1']
  n.dead.cases[ is.na(n.dead.cases) ] = 1
  
  if( n.dead.cases > 3 ) {
    for (predictor in rownames(CAT)){
      X = CAT[predictor, samples2study]
      X <- as.numeric(X)  
      
      mysurfit = survfit(Surv(stime, status) ~ X)  #  summary(mysurfit)
      if( length(unique(gender)) > 1 )
        mysurfit2 = coxph(Surv(stime, status) ~ age+gender+X)
      if( length(unique(gender)) == 1 )   
        mysurfit2 = coxph(Surv(stime, status) ~ age+X)   #  summary(mysurfit2)
      X_pval = summary(mysurfit2)$coefficients['X', 'Pr(>|z|)']
      X_coef = unname(coef(mysurfit2)['X'])
      RESULTS[[predictor]] = c(predictor= predictor, X_pval = X_pval, 'X_coef' = X_coef ) 
      mysurvdiff = survdiff(Surv(stime, status) ~ X) # log-rank test
      X_pval[ is.na(X_pval) ] = 10
      #  summary( mysurfit2 ) ; basehaz(mysurfit2) predictSurvProb
      #  put here for 3 levels discretization
    }
    RESULTS = lapply(RESULTS, function(x) data.frame(t(data.frame(x, stringsAsFactors = F)), stringsAsFactors = F) )
    mRESULTS = melt(RESULTS, id.vars = colnames(RESULTS[[1]]) )
    names(mRESULTS)[ names(mRESULTS) == 'L1' ] = 'predictor'
    mRESULTS$X_adjpval = p.adjust(mRESULTS$X_pval, method = 'fdr')
    return(mRESULTS)
  }
}


guan_rank <- function(x) {
total_rank <- c()
for(i in 1:length(x[,1])) {
  if(x$EVENT[i] == 1) {
    tB_greater_than_tA <- length(which(x$DAYS_TO_LAST_FOLLOWUP > x$DAYS_TO_LAST_FOLLOWUP[i]))  
    
    equal <- which( x$DAYS_TO_LAST_FOLLOWUP == x$DAYS_TO_LAST_FOLLOWUP[i] ) ; equal <- equal[ !equal == i]
    tB_same_tA_sB1 <- length(which(x$EVENT[equal] == 1))/2
    
    tB_lower_than_tA <- intersect( which(x$DAYS_TO_LAST_FOLLOWUP <= x$DAYS_TO_LAST_FOLLOWUP[i]) , which(x$EVENT == 0) ) 
    tB_lower_than_tA <- x[tB_lower_than_tA, "DAYS_TO_LAST_FOLLOWUP" ] 
    
    if(length(tB_lower_than_tA) != 0) {
      ratio <- c()
      for(k in 1:length(tB_lower_than_tA)) {
        survival_A <- length(which(x$DAYS_TO_LAST_FOLLOWUP >= x$DAYS_TO_LAST_FOLLOWUP[i] )  ) 
        survival_B <- length(which(x$DAYS_TO_LAST_FOLLOWUP >= tB_lower_than_tA[k] )   ) 
        ratio <- c(ratio, survival_A/survival_B)
      }
      tB_lower_than_tA <- sum(ratio)
    } else { tB_lower_than_tA <- 0 }
    
    rank <- tB_greater_than_tA + tB_same_tA_sB1 + tB_lower_than_tA
    total_rank <- c(total_rank, rank)
  } else {
    
    tB_greater_than_tA_sB0 <- intersect( which(x$DAYS_TO_LAST_FOLLOWUP >= x$DAYS_TO_LAST_FOLLOWUP[i]) , which(x$EVENT == 0) ) 
    tB_greater_than_tA_sB0 <- x[tB_greater_than_tA_sB0, "DAYS_TO_LAST_FOLLOWUP" ] 
    ratio <- c()
    for(k in 1:length(tB_greater_than_tA_sB0)) {
      survival_A <- length( which(x$DAYS_TO_LAST_FOLLOWUP >= x$DAYS_TO_LAST_FOLLOWUP[i] ) ) 
      survival_B <- length( which(x$DAYS_TO_LAST_FOLLOWUP >= tB_greater_than_tA_sB0[k] ) ) 
      ratio <- c(ratio, 1 - (survival_B/survival_A)*0.5 ) 
    }
    tB_greater_than_tA_sB0 <- sum(ratio)
    
    tB_greater_than_tA_sB1 <- intersect( which(x$DAYS_TO_LAST_FOLLOWUP >= x$DAYS_TO_LAST_FOLLOWUP[i]) , which(x$EVENT == 0) ) 
    tB_greater_than_tA_sB1 <- x[tB_greater_than_tA_sB1, "DAYS_TO_LAST_FOLLOWUP" ] 
    ratio <- c()
    for(k in 1:length(tB_greater_than_tA_sB1)) {
      survival_A <- length( which(x$DAYS_TO_LAST_FOLLOWUP >= x$DAYS_TO_LAST_FOLLOWUP[i] ) ) 
      survival_B <- length( which(x$DAYS_TO_LAST_FOLLOWUP >= tB_greater_than_tA_sB1[k] ) ) 
      ratio <- c(ratio, 1 - (survival_B/survival_A) ) 
    }
    tB_greater_than_tA_sB1 <- sum(ratio)
    
    tB_lower_than_tA_sB0 <- intersect( which(x$DAYS_TO_LAST_FOLLOWUP <= x$DAYS_TO_LAST_FOLLOWUP[i]) , which(x$EVENT == 0) ) 
    tB_lower_than_tA_sB0 <- x[tB_lower_than_tA_sB0, "DAYS_TO_LAST_FOLLOWUP" ] 
    ratio <- c()
    for(k in 1:length(tB_lower_than_tA_sB0)) {
      survival_A <- length( which(x$DAYS_TO_LAST_FOLLOWUP >= x$DAYS_TO_LAST_FOLLOWUP[i] )  ) 
      survival_B <- length( which(x$DAYS_TO_LAST_FOLLOWUP >= tB_lower_than_tA_sB0[k] ) ) 
      ratio <- c(ratio, (survival_A/survival_B)*0.5 ) 
    }
    tB_lower_than_tA_sB0 <- sum(ratio)
    
    rank <- tB_greater_than_tA_sB0 + tB_greater_than_tA_sB1 + tB_lower_than_tA_sB0
    total_rank <- c(total_rank, rank)
  }
}
return(total_rank)
}

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

run_EN_guan_rank <- function(CAT, CLINICS, fold, iteration, alpha, no_event=0, clinvar='DAYS_TO_LAST_FOLLOWUP') {
  result_AUC <- c()
  for (i in 1:iteration) {     ##  start iteration of k fold cross validation
    crossVal <- stratified_crossvalidation(colnames(CAT), nFold=fold, seed=F)
    
    for (i in 1:fold) {
      
      samples2study_train <- as.character(crossVal[[i]]$train) 
      CAT_train <- CAT [ , samples2study_train  ] 
      status_train = (CLINICS[ samples2study_train , ]$EVENT != no_event) + 0
      stime_train = as.numeric(CLINICS[ samples2study_train , clinvar ])
      CLINICS_train <- CLINICS[samples2study_train, ]
      
      # prepare cross-trainset
      samples2study_xTrain <- as.character(crossVal[[i]]$xTrain) 
      CAT_xTrain <- CAT [ , samples2study_xTrain  ] 
      status_xTrain = (CLINICS[ samples2study_xTrain , ]$EVENT != no_event) + 0
      stime_xTrain = as.numeric(CLINICS[ samples2study_xTrain , clinvar ])
      CLINICS_xTrain <- CLINICS[samples2study_xTrain, ]
      
      # train a survival model
      CAT_train <- as.matrix(CAT_train) ;  CAT_train <- t(CAT_train) ; 
      CAT_xTrain <- as.matrix(CAT_xTrain) ;  CAT_xTrain <- t(CAT_xTrain) ; 
      CAT_train[is.na(CAT_train)] <- 0  ## replace NA by 0
      CAT_xTrain[is.na(CAT_xTrain)] <- 0  ## replace NA by 0
      
#     train_rank <- guan_rank(CLINICS_train) ; train_rank <- normalization_01(train_rank) ; names(train_rank) <- rownames(CLINICS_train) #  
#     xTrain_rank <- guan_rank(CLINICS_xTrain); xTrain_rank <- normalization_01(xTrain_rank) ; names(xTrain_rank) <- rownames(CLINICS_xTrain) # 
      
      # adjust rank of training
      train_rank <- guan_rank(CLINICS_train) ; names(train_rank) <- rownames(CLINICS_train)
      ratio_event_train <- table(status_train)[2]/length(status_train) ; limit <- quantile(train_rank, 1-ratio_event_train)
      train_rank_high_risk <- normalization_01(train_rank[train_rank>=limit])/2 + 0.5
      train_rank_low_risk  <- normalization_01(train_rank[train_rank<limit])/2 
      train_rank_adjusted <- c(train_rank_low_risk,train_rank_high_risk)
      train_rank_adjusted <- train_rank_adjusted[names(train_rank)]
      
      # adjust rank of cross training
      xTrain_rank <- guan_rank(CLINICS_xTrain) ; names(xTrain_rank) <- rownames(CLINICS_xTrain)
      ratio_event_xTrain <- table(status_xTrain)[2]/length(status_xTrain) ; limit <- quantile(xTrain_rank, 1-ratio_event_xTrain)
      xTrain_rank_high_risk <- normalization_01(xTrain_rank[xTrain_rank>=limit])/2 + 0.5
      xTrain_rank_low_risk  <- normalization_01(xTrain_rank[xTrain_rank<limit])/2 
      xTrain_rank_adjusted <- c(xTrain_rank_low_risk,xTrain_rank_high_risk)
      xTrain_rank_adjusted <- xTrain_rank_adjusted[names(xTrain_rank)]
      
      # train an elastic net model
      trainOut <- trainEN(CAT_train,train_rank_adjusted, CAT_xTrain,xTrain_rank_adjusted, alpha)
      str(trainOut) 
      
      # extract test set
      samples2study_test <- as.character(crossVal[[i]]$test)
      CAT_test <- CAT [ , samples2study_test] ; CAT_test <- t(CAT_test) ; 
      CLINICS_test <- CLINICS[samples2study_test, ]
      status_test = (CLINICS[ samples2study_test , ]$EVENT != no_event) + 0
      stime_test = as.numeric(CLINICS[ samples2study_test , clinvar ])
      
      # predict with test set
      test_pred <- predictEN(trainOut$model, trainOut$lambda, CAT_test) 
      test_pred <- as.matrix(test_pred)
    
      # Determine AUC
      tdroc <- tdrocc(test_pred, surv.time = stime_test, surv.event=status_test, time=1)
      AUC <- tdroc$AUC
      result_AUC <- c( result_AUC ,   AUC  ) 
      
      plot(x=1-tdroc$spec, y=tdroc$sens, type="l", xlab="1 - specificity",ylab="sensitivity", xlim=c(0, 1), ylim=c(0, 1))
      lines(x=c(0,1), y=c(0,1), lty=3, col="red")
    }
  } ## end of iteration of k fold cross validation
  return(result_AUC)
}

survival_analysis_continuous = function(CAT, CLINICS){
  library(reshape)
  library(survminer)
  status = CLINICS$EVENT
  stime = CLINICS$DAYS_TO_LAST_FOLLOWUP
  age = CLINICS$AGE
  gender = CLINICS$GENDER
  RESULTS = list()
  n.dead.cases = table(status)['1']
  n.dead.cases[ is.na(n.dead.cases) ] = 1

  if( n.dead.cases > 3 ) {
    for (predictor in colnames(CAT)){  ##  predictor="HSH2D"
      X = CAT[  , predictor  ]
      X <- as.numeric(X)  
      
      # mysurfit2 = coxph(Surv(stime, status) ~ X)  #  summary(mysurfit)
      if( length(unique(gender)) > 1 )
        mysurfit2 = coxph(Surv(stime, status) ~ age+gender+X)
      if( length(unique(gender)) == 1 )   
        mysurfit2 = coxph(Surv(stime, status) ~ age+X)   #  summary(mysurfit2)
      X_pval = summary(mysurfit2)$coefficients['X', 'Pr(>|z|)']  #  Cox p-value tests if strata is significantly associated with survival.
      X_coef = unname(coef(mysurfit2)['X'])
      RESULTS[[predictor]] = c(predictor= predictor, X_pval = X_pval, 'X_coef' = X_coef ) 
      # mysurvdiff = survdiff(Surv(stime, status) ~ X) # 
      X_pval[ is.na(X_pval) ] = 10
      #  summary( mysurfit2 ) ; basehaz(mysurfit2) predictSurvProb
      #  put here for 3 levels discretization
    }
    RESULTS = lapply(RESULTS, function(x) data.frame(t(data.frame(x, stringsAsFactors = F)), stringsAsFactors = F) )
    mRESULTS = melt(RESULTS, id.vars = colnames(RESULTS[[1]]) )
    names(mRESULTS)[ names(mRESULTS) == 'L1' ] = 'predictor'
    mRESULTS$X_adjpval = p.adjust(mRESULTS$X_pval, method = 'fdr')
  }
  mRESULTS <- mRESULTS[ ,c(1,2,5,3)]
#   fit <- coxph(Surv(stime, status) ~ age+gender+X)
#   plot(survfit(fit)) 
  return(mRESULTS)
}

survival_KM_plot <- function(data,title,size=25) {
  library(survminer)
  fit <- survfit( coxph ( Surv(DAYS_TO_LAST_FOLLOWUP, EVENT) ~ strata(Biomarker)+AGE+GENDER, data = data) )
  # fit <- survfit( Surv(DAYS_TO_LAST_FOLLOWUP, EVENT) ~ Biomarker, data = data )
  
  size <- size
  ggsurvplot(
    fit,                  # survfit object with calculated statistics.
    data = data,          # data used to fit survival curves.
    risk.table = F,       # show risk table.
    pval = F,             # show p-value of log-rank test.
    pval.size = size-15,
    conf.int = F,         # show confidence intervals for 
    # point estimates of survival curves.
    # xlim = c(0,5000),   # present narrower X axis, but not affect
    # survival estimates.
    xlab = "Time in days",     # customize X axis label.
    break.time.by = 1500,      # break X axis in time intervals by 500.
    ggtheme = theme_light(),   # customize plot and risk table with a theme.
    risk.table.y.text.col = T, # colour risk table text annotations.
    risk.table.y.text = F,     # show bars instead of names in text annotations
    title = title,
    # in legend of risk table
    font.title=size, font.subtitle=size, font.caption=size, font.x=size, font.y=size, font.tickslab=size, font.legend=size
  )
}


# library("survminer") ##  http://www.sthda.com/english/rpkgs/survminer/
# require("survival")
# fit <- survfit(Surv(time, status) ~ sex, data = lung)
# 
# ggsurvplot(fit, data = lung, size=1, censor.size = 6 )  
# 
# text_size <- 20
# ggsurvplot(
#   fit, 
#   data = lung, 
#   size = 1,                 # change line size
#   censor.size = 6,
#   palette = 
#     c("#CC0000", "#0066CC"),# custom color palettes
#   conf.int = TRUE,          # Add confidence interval
#   pval = TRUE,              # Add p-value
#   risk.table = TRUE,        # Add risk table
#   risk.table.col = "strata",# Risk table color by groups
#   legend.labs = 
#     c("Male", "Female"),    # Change legend labels
#   risk.table.height = 0.25, # Useful to change when you have multiple groups
#   risk.table.font = 7 , # text size of the content of table
#   ggtheme = theme(legend.position="none",plot.title=element_text(size=text_size+4,hjust=0.5),
#                   axis.text=element_text(colour="black",size=text_size),
#                   axis.title.x=element_text(colour="black",size=text_size),axis.title.y=element_text(colour="black",size=text_size) )   
# )
# 
# 
# ggsurv <- ggsurvplot(
#   fit,                     # survfit object with calculated statistics.
#   data = lung,             # data used to fit survival curves.
#   risk.table = TRUE,       # show risk table.
#   pval = TRUE,             # show p-value of log-rank test.
#   conf.int = TRUE,         # show confidence intervals for 
#   # point estimates of survival curves.
#   palette = c("#CC0000", "#0066CC"),
#   xlim = c(0,500),         # present narrower X axis, but not affect
#   # survival estimates.
#   xlab = "Time in days",   # customize X axis label.
#   break.time.by = 100,     # break X axis in time intervals by 500.
#   risk.table.y.text.col = T,# colour risk table text annotations.
#   risk.table.font = 7 , # text size of the content of table
#   risk.table.height = 0.25, # the height of the risk table
#   risk.table.y.text = T,# show bars instead of names in text annotations
#   # in legend of risk table.
#   ncensor.plot = TRUE,      # plot the number of censored subjects at time t
#   ncensor.plot.height = 0.25,
#   conf.int.style = "step",  # customize style of confidence intervals
#   surv.median.line = "hv",  # add the median survival pointer.
#   legend.labs = c("Male", "Female"),    # change legend labels.
#   ggtheme = theme(legend.position="none",plot.title=element_text(size=text_size+4,hjust=0.5),
#                   axis.text=element_text(colour="black",size=text_size),
#                   axis.title.x=element_text(colour="black",size=text_size),axis.title.y=element_text(colour="black",size=text_size) )   
# )
# ggsurv









