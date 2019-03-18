library("Hmisc", lib.loc="/usr/lib/R/site-library")
library(hydroGOF)

#-------------------------------------------------------------------------------------------------------------------------------
#DESC: This function creates a plot of pred vs observations. First of all the function loads the prediction matrix of the 
#      commited feature set and also creates the observations matrix. Then all cindexes for one drug are calculated.
#IN:   drug_id         ==> the drug_id of the drug for which the plot will be generated
#      search_path     ==> the path to the drug_folder where the predictions are stored
#OUT: returns a plot 
#-------------------------------------------------------------------------------------------------------------------------------
create_pred_plot <- function(drug_id, search_path){
  
  # defining the drug_name fitting to the current drug_id
  drug_name <- as.character(DRUG_MASTER_LIST$DRUG_NAME[DRUG_MASTER_LIST$DRUG_ID==drug_id])
  
  # go into the folder of one feature
  path_to_feature_folder <- file.path(search_path)
  # load predMat
  load(paste(path_to_feature_folder, "/predictions.RData", sep=""))
  
  # load observation matrix
  load("data/pipeline/drug/DRUG_RESPONSE.ro")
  DRUG_RESPONSE_RELEASED <- t(OMAUC_RELEASED)
  
  # define predictions for the drug with the committed drug_id
  pred <- predMat[, colnames(predMat)==drug_id]
  
  # define observations for the drug with the committed drug_id
  obs <- DRUG_RESPONSE_RELEASED[, colnames(DRUG_RESPONSE_RELEASED)==drug_id]
  
  # calculate performance on test set and round result on three digits
  cIdx_test <- formatC(rcorr.cens(obs, pred)[1])
  
  # calculate Spearman-correlation coefficient and round result on three digits
  value_spearman <- formatC(cor(obs, pred, method="spearman", use ="na.or.complete"))
  
  # calculate RMSE and round result on three digits
  value_rmse <- formatC(rmse(obs, pred, na.rm=TRUE)) 
  
  # create the plot 
  plot(obs, pred, 
       xlab="Obs",
       ylab="Pred",
       main=drug_name, 
       sub=paste("cIdx=",cIdx_test, ", spearman=",value_spearman, ", RMSE=",value_rmse, sep=""))
}