#-----------------------------------------------------------------------------------------------------------------------------------------
# This script is the main script needed in order to run a single drug analysis with a crossvalidation for the input data.
# It reads the input from the bash script it's been called from and loads and manipulates the necessary data according to those arguments.
# This script performs a crossvalidation for the input data and trains a statistical modell with a algorithm according to the 6th
# committed argument. A matrix containing all the predicted values will be stored for each feature and the 10 best models for the 
# statistical algorithm used.
# List of committed paramaters:
#       - 1st arg: a boolean value. Only if this value is TRUE the trainSet will be balanced
#       - 2nd arg: a numeric value that decides which kind of crossvalidation will be done
#       - 3rd arg: a string containing the features, which should be used for the training. If more than one feature set
#                  is committed, they need to be seperated by "_"(underline).
#       - 4th arg: a string containing a tissue type. If the committed string is not "all_tissue" a tissue specific analysis
#                  with the committed tissue type will be done, otherwise a non specific.
#       - 5th arg: a string containing the path, where all the created data will be stored
#       - 6th arg: a string deciding which statistical model will be used
#       - 7th arg: a string deciding which matrix for the drug reponse will be used
#       - 8th arg: a boolean value, if it is TRUE a seed will be set for the crossvalidation
#-----------------------------------------------------------------------------------------------------------------------------------------
setwd("/nfs/nobackup2/saezgrp/homes/micha")

# ###############################################################################################
# Load libraries
# ###############################################################################################

library(randomForest)
library(Hmisc)
library(glmnet)
# #############################################################################
# Source functions
# #############################################################################

source("script/R/ML/loadFeatureSet.R")
source("loadDrug_res.R")
source("script/R/ML/createFolders_shrinkToTissue")
source("script/R/ML/CrossVal_LIBRARY.R")
source("script/R/ML/StatModel_LIBRARY.R")
source("script/R/ML/LIB_METRICS.R")
source("script/R/ML/LIB_BALANCING.R")
# ###############################################################################################
# Get values commited from the bash script
# ###############################################################################################

# gets the arguments given to the script
args <- commandArgs(trailingOnly=TRUE)

# define parameters according to the values commited from the bash script 
balancing <- args[1]
nFold <- args[2]
nFold <-  as.numeric(nFold)
featureSet <- args[3]
tissue_spec <- args[4]
general_path <- args[5]
stat_model <- args[6]
drug_mat <- args[7]
seed <-  args[8]

# ###############################################################################################
# Load neccessary drug matrices and check input for stat_model
# ###############################################################################################

drug_res <- loadDrug_res(drug_mat, stat_model)
#################################################################################################
# Define the features matrix according to the input for featureSet
#################################################################################################

features <- loadData(featureSet)
# #############################################################################
# Decide which kind of tissue specific analysis should be done
# and create general folder structure
# #############################################################################

result <- getFolders_shrinkToTissue(tissue_spec, drug_res, features, general_path, drug_mat, featureSet, balancing)

# define drug_res(might have been shrinked for tissue specific analysis)
drug_res <- result$drug_res

# define feature matrix(might have been shrinked for tissue specific analysis)
features <- result$features

# set the feature_folder
feature_folder <- result$feat_folder
# #############################################################################
# Crossvalidation for all Drugs
# #############################################################################

# declaring a list, which is going to be filled with 
# the crossvaledated cosmic id's of the cell-lines
crossVal <- list()

for(i in 1:nrow(drug_res)){
  
  # remove cell lines without drug response
  res <-drug_res[i,]
  cells <- as.character(names(res[!is.na(res)]))
  
  # using only cells with drug response and complete feature
  cell_withFeat <- rownames(features)[!is.na(rowSums(features))]
  cells <- intersect(cells, cell_withFeat)
  
  # DON'T DO CROSSVALIDATION WITH res<2*nFold and enter NA in crossval
  if (length(cells) < 2*nFold){
    crossVal[[as.character(rownames(drug_res)[i])]] <- NA
  }else{
    
    #doing crossvalidation
    crossVal[[as.character(rownames(drug_res)[i])]] <- crossvalidation(as.numeric(cells), nFold, seed)
    
    # if there is no variance in train nor in xtrain for any fold, set to NA in crossval
    for (fold in 1:nFold){
      if(length(unique(res[as.character(crossVal[[as.character(rownames(drug_res)[i])]][[fold]]$train)]))==1){
        crossVal[[as.character(rownames(drug_res)[i])]] <- NA
        break
      } else if (length(unique(res[as.character(crossVal[[as.character(rownames(drug_res)[i])]][[fold]]$xTrain)]))==1){
        crossVal[[as.character(rownames(drug_res)[i])]] <- NA
        break
      }
    }
  }
}

# check if only NA values are inside crossVal, throw Exception, if that is the case
if(!any(!is.na(crossVal))){
  stop(paste("There is not one drug, which has a given res>2*nFold for ", tissue_spec, sep=""))
}

# #############################################################################
# Training a Model
# #############################################################################

# create a new matrix with exactly the same size, colnames and rownames as drug_res matrix but filled with "NA"
predMat <- matrix(NA, nrow=nrow(drug_res), ncol=ncol(drug_res))
colnames(predMat) <- colnames(drug_res)
rownames(predMat) <- rownames(drug_res)
predMat <- t(predMat)

# train for each drug
for(d_Idx in 1:nrow(drug_res)){
  
  # only train a statistical model, if a crossvalidation has been done for the drug with the id d_Idx
  if(is.list(crossVal[[d_Idx]])){
    
    # getting the drug response only for that one drug
    res <- drug_res[d_Idx,]
    res <- res[!is.na(res)]
      
    # generate a directory for each drug
    path_to_drug_folder <- file.path(feature_folder, paste("drug", rownames(drug_res)[d_Idx], sep=""))
    dir.create(path_to_drug_folder, showWarnings=FALSE)
  
    # generate a directory for the statistical method that will be used later on
    path_to_statModel_folder <- file.path(path_to_drug_folder, stat_model)
    dir.create(path_to_statModel_folder, showWarnings=FALSE)
    
    for(n in 1:nFold){
    
      train_obs <- res[as.character(crossVal[[d_Idx]][[n]]$train)]
      # check if trainset should be balanced or not
      if(balancing){
        train_obs <- balance(train_obs)
      }
      train_feat <- features[names(train_obs), ]
      rownames(train_feat) <- names(train_obs)
      
      # check if there's variance in train_feat 
      if(length(unique((as.vector(train_feat)))) == 1){
        
        # delete drug folder
        unlink(path_to_drug_folder, recursive=TRUE)
        
        # save NA values for this drug in predMat
        predMat[, rownames(drug_res)[d_Idx]] <- NA
        
        # quit the loop and end training for this drug
        break
      }
  
      # prepare cross-trainset
      xTrain_obs <- res[as.character(crossVal[[d_Idx]][[n]]$xTrain)]
      xTrain_feat <- features[names(xTrain_obs), ]
      rownames(xTrain_feat) <- names(xTrain_obs)
      
      # prepare testset
      test_obs <- res[as.character(crossVal[[d_Idx]][[n]]$test)]
      test_feat <- features[names(test_obs), ]
      rownames(test_feat) <- names(test_obs)
      
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
      if(stat_model=="EN"){
        
        # train an elasticNet model
        trainOut  <- trainEN(train_feat, train_obs, xTrain_feat, xTrain_obs)
        
        # predict with test set
        test_pred <- predict(trainOut$model,type="response",newx=test_feat)[,trainOut$model$lambda == trainOut$lambda]
      }
      
      # save trainOut
      save_path_trainOut <- file.path(path_to_statModel_folder, paste("trainOut_", n, ".RData", sep=""))
      save(trainOut,file=save_path_trainOut)
      
      # save predicted values in to the matrix predMat
      predMat[names(test_pred), rownames(drug_res)[d_Idx]] <- test_pred
    }
  }else{
    # in case there are less than nFold training samples
    predMat[, rownames(drug_res)[d_Idx]] <- NA
  }
}
#save predMat to the folder of the current feature 
save_path_predMat <- file.path(feature_folder, paste("predictions_", stat_model, ".RData", sep="" ))
save(predMat, file=save_path_predMat)
