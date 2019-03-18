# ###############################################################################################
# Get values commited from the bash script
# ###############################################################################################

start_time <- Sys.time()

# gets the arguments given to the script
args <- commandArgs(trailingOnly=TRUE)

# define parameters according to the values commited from the bash script 
balancing <- args[1]
cellFeatureSet <- args[2]
drugFeatureSet <- args[3]
tissue_spec <- args [4]
output_path <- args[5]
src_path <- args[6]
stat_model <- args[7]
drug_mat <- args[8]
cv_fold <- args[9]
cv_fold <- as.numeric(cv_fold)
cv_idep <- args[10]
cv_method <- args[11]
cv_foldNR <- args[12]
cv_foldNR <- as.numeric(cv_foldNR)

# ###############################################################################################
# Load libraries
# ###############################################################################################

library(randomForest)
library(Hmisc)
library(glmnet)

# #############################################################################
# Source functions
# #############################################################################
setwd(src_path)

source("script/R/ML/loadDrugFeatSet_multi.R")
source("script/R/ML/loadFeatSet_multi.R")
source("script/R/ML/loadDrug_res.R")
source("script/R/ML/createFolders_shrinkToTissue.R")
source("script/R/ML/LIB_StatModel.R")
source("script/R/ML/LIB_METRICS.R")
source("script/R/ML/LIB_BALANCING.R")
source("script/R/crossvalidate_data/03_multi_drug_CV.R")
source("script/R/crossvalidate_data/CrossVal_LIBRARY.R")

# ###############################################################################################
# Load drug response matrix
# ###############################################################################################

drug_res <- loadDrug_res(drug_mat)

# ###############################################################################################
# Load drug feature sets
# ###############################################################################################

# define combined drug feature matrix
drug_features <- loadDrugData(drugFeatureSet)

# shrink drug_res to those the drugId's, which are included in drug_features and drug_res
drug_res <- drug_res[intersect(rownames(drug_features), rownames(drug_res)), ]

# ###############################################################################################
# Load cell feature sets
# ###############################################################################################

# define combined cell feature matrix with all commited feature sets
cell_features <- loadCellData(cellFeatureSet)

# shrink drug_res to the cellId's, which are included in the cell features
drug_res <- drug_res[, rownames(cell_features)]

# ###############################################################################################
# Decide which kind of tissue specific analysis should be done
# and create general folder structure
# ###############################################################################################

result <- getFolders_shrinkToTissue(tissue_spec, drug_res, cell_features, output_path, drug_mat, 
                                    paste(cellFeatureSet, drugFeatureSet, sep=":"), balancing)


#print(paste("Loading completed time needed: ", as.numeric(difftime(Sys.time(), start_time), units="mins"), "mins" ))

# define drug_res(might have been shrinked for tissue specific analysis)
drug_res <- result$drug_res

# define feature matrix(might have been shrinked for tissue specific analysis)
cell_features <- result$features

# set the feature_folder
feature_folder <- result$feat_folder

# crossproduct of drug_res
cross_Ids <- as.matrix(expand.grid(rownames(drug_res), colnames(drug_res)))
colnames(cross_Ids) <- c("drug_id", "cell_id")

# create empty train matrix
trainMat <- c()

#create empty observation vector
obsVec <- c()

start_time <- Sys.time()

#obsVec:
obsVec <- expand.grid(drug_res)
obsVec <- as.numeric(as.matrix(obsVec))
names(obsVec) <- paste(cross_Ids[,"drug_id"], "_", cross_Ids[, "cell_id"], sep="") 

#trainMat:
trainMat <- cbind(drug_features[cross_Ids[, "drug_id"], ], cell_features[cross_Ids[, "cell_id"], ])
rownames(trainMat) <-  paste(cross_Ids[,"drug_id"], "_", cross_Ids[, "cell_id"], sep="") 

#print(paste("cross_Ids loop completed time needed: ", as.numeric(difftime(Sys.time(), start_time), units="mins"), "mins" ))


#start_time <- Sys.time()

###############################################################################################
# Crossvalidation
###############################################################################################

# check which entity should be independent
if(cv_idep == "drug"){
  
  cv <- multicross(as.numeric(rownames(drug_res)), as.numeric(colnames(drug_res)), cv_fold, cv_method, c("drug_id", "cell_id"))
} else {
  
  cv <- multicross(as.numeric(colnames(drug_res)), as.numeric(rownames(drug_res)), cv_fold, cv_method, c("cell_id", "drug_id"))
}

#print(paste("CV completed time needed: ", as.numeric(difftime(Sys.time(), start_time), units="mins"), "mins" ))

#################################################################################################
# Train section
#################################################################################################

# generate a directory for the statistical method that will be used later on
path_to_statModel_folder <- file.path(feature_folder, stat_model)
dir.create(path_to_statModel_folder, showWarnings=FALSE)

#start_time <- Sys.time()
# train for each cv
  
  # prepare train set
  if(cv_idep == "drug") {
    
    train_obs <- obsVec[as.character(paste(cv[[cv_foldNR]]$train[, 1], cv[[cv_foldNR]]$train[, 2], sep = "_"))]
  } else {
    
    train_obs <- obsVec[as.character(paste(cv[[cv_foldNR]]$train[, 2], cv[[cv_foldNR]]$train[, 1], sep = "_"))]
  }
  
  # remove NA values and entries in train_obs
  train_obs <- train_obs[complete.cases(train_obs)]
  
  # check if trainset should be balanced or not
  if(balancing){
    train_obs <- balance(train_obs)
  }
  
  train_feat <- trainMat[names(train_obs), ]
  rownames(train_feat) <- names(train_obs)
  
  
  # prepare cross-trainset
  if(cv_idep == "drug") {
    
    xTrain_obs <- obsVec[as.character(paste(cv[[cv_foldNR]]$xTrain[, 1], cv[[cv_foldNR]]$xTrain[, 2], sep = "_"))]
  } else {
    
    xTrain_obs <- obsVec[as.character(paste(cv[[cv_foldNR]]$xTrain[, 2], cv[[cv_foldNR]]$xTrain[, 1], sep = "_"))]
  }
  
  # remove NA values and entries
  xTrain_obs <- xTrain_obs[complete.cases(xTrain_obs)]
  
  xTrain_feat <- trainMat[names(xTrain_obs), ]
  rownames(xTrain_feat) <- names(xTrain_obs)
  
  
  # prepare testset
  if(cv_idep == "drug") {
    
    test_obs <- obsVec[as.character(paste(cv[[cv_foldNR]]$test[, 1], cv[[cv_foldNR]]$test[, 2], sep = "_"))]
  } else {
    
    test_obs <- obsVec[as.character(paste(cv[[cv_foldNR]]$test[, 2], cv[[cv_foldNR]]$test[, 1], sep = "_"))]
  }
  
  # remove NA values and entries
  test_obs <- test_obs[complete.cases(test_obs)]
  
  test_feat <- trainMat[names(test_obs), ]
  rownames(test_feat) <- names(test_obs)
  
  #print(paste("TrainSets completed time needed: ", as.numeric(difftime(Sys.time(), start_time), units="mins"), "mins" ))
  
  # #######################################################################
  # Train a randomForest model
  # #######################################################################
  if(stat_model=="RF"){
    
    # train a randomForest model
    trainOut <- trainRF(train_feat, train_obs, xTrain_feat, xTrain_obs)
    
    # predict with test set
    test_pred <- predict(trainOut$model, test_feat)
  }
  
  # #######################################################################
  # Train an elasticNet model
  # #######################################################################
  #start_time <- Sys.time()
  if(stat_model=="EN"){
    
    # train an elasticNet model
    trainOut  <- trainEN(train_feat, train_obs, xTrain_feat, xTrain_obs)
    
    # predict with test set
    test_pred <- predict(trainOut$model,type="response",newx=test_feat)[,trainOut$model$lambda == trainOut$lambda]
  }
  #print(paste("Training completed time needed: ", as.numeric(difftime(Sys.time(), start_time), units="mins"), "mins" ))
  
  
  # save trainOut
  if(cv_method == "random") {
    
    save_path_trainOut <- file.path(path_to_statModel_folder, paste("trainOut_", stat_model, "_random_", cv_foldNR, ".RData", sep=""))
  } else {
    
    save_path_trainOut <- file.path(path_to_statModel_folder, paste("trainOut_", stat_model, "_leave_out_", cv_idep, "_", cv_foldNR, ".RData", sep=""))
  }
  save(trainOut, file=save_path_trainOut)

#save predMat to the folder of the current feature 
if(cv_method == "random") {

  save_path_predVec <- file.path(path_to_statModel_folder, paste("predVec_", stat_model, "_random_", cv_foldNR, ".RData", sep=""))
} else {
  
  save_path_predVec <- file.path(path_to_statModel_folder, paste("predVec_", stat_model,"_leave_out_", cv_idep, "_", cv_foldNR, ".RData", sep=""))
}
save(test_pred, file=save_path_predVec)

#write and save a text file containing the execution time of this script
#end_time <- Sys.time()

#elapsed_time <- end_time-start_time

#write(as.numeric(difftime(end_time, start_time), units="mins"), file = paste(path_to_statModel_folder, "/timeFile.txt", sep=""))
