#-----------------------------------------------------------------------------------------------------------------------------------------
# This script is the main script needed in order to run a single drug analysis with bootstrap. This script runs the analysis only
# for one drug, which one is used is decided by the commited argument number 7.
# It takes input from the bash file it has been called from, trains a statistical model with the algorithm according to argument
# number 5. After that for each feature type all drugs are stored in a folder containing a folder for the statistical model, where
# a matrix with the feature importance and a matrix with the performance values for the predicted values.
# List of committed paramaters:
#       - 1st arg: a boolean value. Only if this value is TRUE the trainSet will be balanced
#       - 2nd arg: a string containing the features, which should be used for the training. If more than one feature set
#                  is committed, they are seperated by "_"(underline).
#       - 3th arg: a string containing a tissue type. If the committed string is not "all_tissue" a tissue specific analysis
#                  with the committed tissue type will be done, otherwise a non specific.
#       - 4th arg: a string containing the path, where all the created data will be stored
#       - 5th arg: source path
#       - 6th arg: a string deciding which statistical model will be used
#       - 7th arg: a string deciding which matrix for the drug reponse will be used
#       - 8th arg: a numeric value that decides, which row of the drug_res matrix will be used
#-----------------------------------------------------------------------------------------------------------------------------------------

# ###############################################################################################
# Get values commited from the bash script
# ###############################################################################################

# gets the arguments given to the script
args <- commandArgs(trailingOnly=TRUE)

# define parameters according to the values commited from the bash script 
balancing <- args[1]
featureSet <- args[2]
tissue_spec <- args[3]
outputPath <- args[4]
src_path <- args[5]
stat_model <- args[6]
drug_mat <- args[7]
rownumber <- as.numeric(args[8])

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

source("~/Documents/RWTH Aachen/ML CODE/ML/loadFeatureSet.R")
source("~/Documents/RWTH Aachen/ML CODE/ML/loadDrug_res.R")
source("~/Documents/RWTH Aachen/ML CODE/ML/createFolders_shrinkToTissue.R")
source("~/Documents/RWTH Aachen/ML CODE/ML/LIB_StatModel.R")
source("~/Documents/RWTH Aachen/ML CODE/ML/LIB_METRICS.R")
source("~/Documents/RWTH Aachen/ML CODE/balancing/LIB_BALANCING.R")

# ###############################################################################################
# Load neccessary drug matrices and check input for stat_model
# ###############################################################################################

drug_res <- loadDrug_res(drug_mat)
#################################################################################################
# Define the features matrix according to the input for featureSet
#################################################################################################

features <- loadData(featureSet)
# #############################################################################
# Decide which kind of tissue specific analysis should be done
# and create general folder structure
# #############################################################################

result <- getFolders_shrinkToTissue(tissue_spec, drug_res, features, outputPath, drug_mat, featureSet, balancing)

# define drug_res(might have been shrinked for tissue specific analysis)
drug_res <- result$drug_res

# define feature matrix(might have been shrinked for tissue specific analysis)
features <- result$features

# set the feature_folder
feature_folder <- result$feat_folder
# #############################################################################
# Training section 
# ############################################################################

# how many models should be build
nSteps <- 10

# define the drug_id corresponding to the commited rownumber
drug_id <-  rownames(drug_res)[rownumber]

# remove cell lines without drug response
drug_res <- drug_res[drug_id,]
cells <- as.character(names(drug_res[!is.na(drug_res)]))

# using only cells with drug response and complete feature
cell_withFeat <- rownames(features)[!is.na(rowSums(features))]
cells <- intersect(cells, cell_withFeat)

# DON'T DO TRAINING WITH res<2*nFold and stop the script for this drug
if (length(cells) < 10){
  stop(paste("less than 10 cell lines have a response for drug_id:", drug_id, sep=""))
} else{

  # define size of one bin for boot strapping (test set size)
  binSize <- ceil(length(cells) / 10)
  
  # for evaluating performance, save output
  predMat <- matrix(NA, ncol=binSize, nrow=nSteps)
  obsMat <- matrix(NA, ncol=binSize, nrow=nSteps)
  cellIDmat <- matrix(NA, ncol=binSize, nrow=nSteps)
  
  # create a matrix "impMat" with the featurenames of the features matrix
  impMat <- matrix(NA, ncol=ncol(features), nrow=nSteps)
  colnames(impMat) <- colnames(features)
  
  # train 1000 models if possible
  for (i in 1:(nSteps+1)){

    # train as long as i is smaller than 1000
    if(i<=nSteps){
      
      # sample COSMIC ID'S
      ids <- sample(cells)
         
      # devide ids into three samples
      test <- ids[1:binSize]
      xTrain <- ids[(binSize+1):(2*binSize)]
      train <- ids[(2*binSize+1):length(ids)]
      
      #####################################################################################################################
      # TEST section for the bootstrapping
      #####################################################################################################################
      
      #TEST if all COSMIC ID's have been used
      if(sum(length(test),length(train),length(xTrain))!=length(cells)) 
          stop(paste("not all COSMIC ID's are either used in train, xTrain or test in step", i, sep=""),call. =FALSE)
      
      #TEST that there are no equal COSMIC ID's in train and xtrain
      if(length(intersect(train,xTrain))>0) 
        stop(paste("train and xTrain contain the same values in step:", i, sep=""))
      
      #TEST that there are no equal COSMIC ID's in train and test
      if(length(intersect(train,test))>0) 
        stop(paste("train and test contain the same values in step:", i, sep=""))
      
      #TEST that there are no equal COSMIC ID's in xTrain and test
      if(length(intersect(xTrain,test))>0) 
        stop(paste("xTrain and test contain the same values in step:", i, sep=""))
      
      #TEST if there's variance in the train set
      if(length(unique(drug_res[as.character(train)]))==1){
        # go to the next model for this drug with new sampled cell lines
        next
      }
      
      #TEST if there's variance in the xTrain set
      if(length(unique(drug_res[as.character(xTrain)]))==1){
        # go to the next model for this drug with new sampled cell lines
        next
      }
      
      ##############################################################################################################
      # Train statistical model
      ##############################################################################################################
      
      # prepare trainset
      train_obs <- drug_res[as.character(train)]
      # check if trainset should be balanced or not
      if(balancing){
        train_obs <- balance(train_obs)
      }
      train_feat <- features[names(train_obs), ]
      rownames(train_feat) <- names(train_obs)
      
      # check if there's variance in train_feat 
      if(length(unique((as.vector(train_feat)))) == 1){
        # go to the next model for this drug with new sampled cell lines
        next
      }
      
      # prepare cross-trainset
      xTrain_obs <- drug_res[as.character(xTrain)]
      xTrain_feat <- features[names(xTrain_obs), ]
      rownames(xTrain_feat) <- names(xTrain_obs)
      
      # prepare testset
      test_obs <- drug_res[as.character(test)]
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
        
        # define best model importance
        best_model_importance <- trainOut$model$importance
      }
      
      # #######################################################################
      # Train an elasticNet model
      # #######################################################################
      if(stat_model=="EN"){
        # train an elasticNet model
        trainOut  <- trainEN(train_feat, train_obs, xTrain_feat, xTrain_obs)
        
        # predict with test set
        test_pred <- predict(trainOut$model,type="response",newx=test_feat)[,trainOut$model$lambda == trainOut$lambda]
        
        # define best model importance
        best_model_importance <- as.matrix(trainOut$featWeights) * apply(features[cells,], 2, sd)
      }
      
      # #######################################################################
      # store output matrix
      # #######################################################################
      
      predMat[i, ] <- test_pred
      obsMat[i, ] <- test_obs
      cellIDmat[i, ] <- names(test_obs)
    
      # #######################################################################
      # Fill feat importance matrix
      # #######################################################################
 
      impMat[i, rownames(best_model_importance)] <- best_model_importance
          
    } else{  
      # #######################################################################
      # Create folder structure for each drug
      # #######################################################################
      
      # generate a directory for each drug
      path_to_drug_folder <- file.path(feature_folder, paste("drug", drug_id, sep=""))
      dir.create(path_to_drug_folder, showWarnings=FALSE)
      
      # generate a directory for the statistical method 
      path_to_statModel_folder <- file.path(path_to_drug_folder, stat_model)
      dir.create(path_to_statModel_folder, showWarnings=FALSE)
      
      # save performance matrix
      save_path_prediction <- file.path(path_to_statModel_folder, "pred.ro")
      save(predMat, obsMat, cellIDmat, file=save_path_prediction)
      
      # calculate the mean and the standard deviation of the impMat
      vMean <- apply(impMat, 2, mean)
      vStd <- apply(impMat, 2, sd)
      
      # save the standard deviation and the mean in a matrix
      effectSize <- cbind(vMean, vStd)
      # save featImp mat 
      save_path_imp <- file.path(path_to_statModel_folder, "effectSize.ro")
      save(effectSize, file=save_path_imp)
    }
  }
}