setwd("~/Documents/RWTH Aachen/ML CODE/GDSC_web_v5.0")
library(glmnet)

drug <- "Nutlin-3a (-)"

# ###############################################################################################
# Elastic net functions
# ###############################################################################################

# train elastic net
# IN: copy from other training method.. 
# OUT: model       <= Elastic net model
#      lambda      <= Chosen lambda for further predictions (parameter which has been optimize)
#      xTrain_cIDX <= Concordance indices of the cross-train data
#      train_cIDX  <= Concordance indices of the train data
#      featWeights <= weigths for each feature of the chosen model
trainEN <- function(train_feat, train_obs, xTrain_feat, xTrain_obs) {
  fit <- glmnet(train_feat, train_obs)
  
  xTrain_pred <- predict(fit,type="response",newx=xTrain_feat)
  train_pred <- predict(fit,type="response",newx=train_feat)
  
  xTrain_cIDX <- c()
  train_cIDX <- c()
  for (i in 1:ncol(fit$beta)) {
    xTrain_cIDX <- c(xTrain_cIDX, cIDX(xTrain_pred[,i], xTrain_obs))
    train_cIDX <- c(train_cIDX, cIDX(train_pred[,i], train_obs))
  }
  lambda <- fit$lambda[tail(which(xTrain_cIDX == max(xTrain_cIDX)), n=1)]
  
  return(list(model=fit,
              lambda=lambda,
              xTrain_cIDX=xTrain_cIDX,
              train_cIDX=train_cIDX,
              featWeights=fit$beta[,lambda==fit$lambda]))
}

# train elastic net
# IN:  model  <= Elastic net model
#      lambda <= Feature 
# OUT: prediction <= Prediction of feature with given model and chosen lambda.
predictEN <- function(model, lambda, feat) {
  return(predict(model,type="response",newx=feat)[,model$lambda == lambda])
}

# ###############################################################################################
# Load neccessary matrices
# ###############################################################################################

OMAUC_RELEASED <- t(read.csv("drugRes_OMAUC.csv", row.names=1, check.names=F))
DRUG_MASTER_LIST <- read.csv("drug_info.csv")

features <- read.csv("bem.csv", row.names=1, check.names=F)

source("~/Documents/RWTH Aachen/ML CODE/crossvalidate_data/CrossVal_LIBRARY.R")
source("~/Documents/RWTH Aachen/ML CODE/ML/LIB_METRICS.R")

# #############################################################################
# Crossvalidation for all Drugs
# ############################################################################

drugID <- DRUG_MASTER_LIST$DRUG_ID[DRUG_MASTER_LIST$DRUG_NAME==drug]

# only analyse data where we have response
res <- OMAUC_RELEASED[as.character(drugID), ]
res <- res[!is.na(res)]

# only cell lines from NCI-60  really ?!
str(res)
str(features)
res <- res[intersect(names(res), rownames(features))] ## those names are cell line COSMIC ID

# cross-validate
crossVal <- crossvalidation(as.numeric(names(res)), nFold=10, seed=123)

# #############################################################################
# Training a ElasticNet model
# #############################################################################

# only do first fold
train_obs <- res[as.character(crossVal[[1]]$train)]
train_feat <- as.matrix(features[names(train_obs), ])
rownames(train_feat) <- names(train_obs)

# prepare cross-trainset
xTrain_obs <- res[as.character(crossVal[[1]]$xTrain)]
xTrain_feat <- as.matrix(features[names(xTrain_obs), ])
rownames(xTrain_feat) <- names(xTrain_obs)

# train a randomForest model
trainOut <- trainEN(train_feat, train_obs, xTrain_feat, xTrain_obs)
str(trainOut$model)
# #############################################################################
# plots after training
# #############################################################################

# plot elastic net with all lambda
nFeat <- 5 # how many feature to plot as parameter
plot(trainOut$model, 
     xvar="lambda", 
     label=TRUE, 
     col=1:8, 
     main=drug)
str(trainOut$model$beta)
varSeqNumber <- sort(rowSums(trainOut$model$beta == 0), index.return=TRUE)$ix
legend("topright",
       legend=paste("(", varSeqNumber[1:nFeat], ") ", rownames(trainOut$model$beta)[varSeqNumber[1:nFeat]], sep=""), 
       lty=0,cex=.6)

# plot cross-training
plot(log(trainOut$model$lambda), trainOut$train_cIDX, type="l", 
     main="cross-training", xlab="log(lambda)", ylab="c-index")
points(log(trainOut$model$lambda), trainOut$train_cIDX)
lines(log(trainOut$model$lambda), trainOut$xTrain_cIDX, col="red")
points(log(trainOut$model$lambda), trainOut$xTrain_cIDX, col="red")
abline(v=log(trainOut$lambda), col="purple", lwd=3)
legend("bottomleft", 
       legend=c("training", "cross-training", "chosen lambda"), 
       lty=1, lwd=c(1,1,3), pch=c(1,1,NA),
       col=c("black", "red", "purple"))

# plot feature weigths plot (TODO, make more pretty) ;-)
#barplot(sort(trainOut$featWeights),las=2, cex.names=0.5, cex.main=1, cex.lab=3, main=drug)

# feature importances plot
#barplot(sort(abs(trainOut$featWeights)),las=2, cex.names=0.5, main=drug)

# #######################################################################
# Predict with trained model
# #######################################################################

# extract test set
test_obs <- res[as.character(crossVal[[1]]$test)]
test_feat <- as.matrix(features[names(test_obs), ])
rownames(test_feat) <- names(test_obs)

# predict with test set
test_pred <- predictEN(trainOut$model, trainOut$lambda, test_feat)

# prediction versus observation plot
plot(test_pred, test_obs)

obs <- unname(test_obs, force = FALSE)
pred <- unname(test_pred, force = FALSE)
cor(obs, pred)

