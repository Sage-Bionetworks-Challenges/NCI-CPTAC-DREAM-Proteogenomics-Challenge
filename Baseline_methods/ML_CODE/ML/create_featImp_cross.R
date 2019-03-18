#-------------------------------------------------------------------------------------------------------------------------------
#DESC:  This function creates a barplot with feature importance for analysis, where a crossvalidation has been performed. 
#       It's collecting the featureImportance from all ten models of a 
#       certain drug and calculates the average of it. After this the barplot will be created and saved at the current working
#       directory with the name similiar to the value of "plot_name".
#IN:    drug_id         ==> the drug_id of the drug for which the plot will be generated
#       search_path     ==> a string containing the path to the folder, which contains all drug folders
#       plot_name       ==> the name of the stored plot
#       stat_model      ==> the name of the statistical model, that has been used to built the results
#       save_path       ==> the path, where the plot will be saved
#       tissue_spec     ==> a string defining the tissue_type, that has been used for the model
#       balance         ==> a string containing, whether the train_set has been balanced or not
#OUT:   returns a plot named after the value of "plot_name"
#-------------------------------------------------------------------------------------------------------------------------------
create_feature_imp <- function(drug_id, search_path, plot_name, stat_model, save_path, tissue_spec, balance){

  # defining the drug_name fitting to the current drug_id
  drug_name <- as.character(DRUG_MASTER_LIST$DRUG_NAME[DRUG_MASTER_LIST$DRUG_ID==drug_id])
  
  # go into the folder of one feature
  path_to_feature_folder <- file.path(search_path)

  # go into the folder of one drug
  path_to_drug_folder <- file.path(path_to_feature_folder, paste("drug", drug_id, sep=""))

  # go into the folder for models
  path_to_model <- file.path(path_to_drug_folder, stat_model)
  
  # declaring the number of how many models are there for one drug
  nFold <- length(list.files(path=path_to_model))

  # loop over every model stored
  for(n in 1:nFold){
    
    if(stat_model=="RF"){
      
      # load a model 
      load(paste(path_to_model, "/trainOut_", n, ".RData", sep=""))
      
      # define best model importance
      best_model_importance <- trainOut$model$importance
    }
    
    if(stat_model=="EN"){
      
      # load a model 
      load(paste(path_to_model, "/trainOut_", n, ".RData", sep=""))
      
      # define best model importance
      best_model_importance <- as.matrix(trainOut$featWeights)
    }
    
    
    # declare only for the first run of the loop a matrix "impMat"
    if(n==1){
      # create a matrix "impMat" with ncol=nrow(best_model_importance) and nrow=nFold
      impMat <- matrix(ncol=nrow(best_model_importance), nrow=nFold)
      # declare the colnames of "impMat" as the rownames of "best_model_importance"
      colnames(impMat) <- rownames(best_model_importance)
    }
    
    # fill the "impMat" with the values of the "best_model_importance"
    impMat[n, rownames(best_model_importance)] <- best_model_importance
  }
  
  # create a matrix "mean_impMat" with the same size as "best_model_importance"
  mean_impMat <- matrix(nrow=nrow(best_model_importance), ncol=2)
  # declare the rownames of "mean_impMat" as the colnames of "impMat"
  rownames(mean_impMat) <- colnames(impMat)
  # declaring the colnames of "mean_impMat"
  colnames(mean_impMat) <- c("mean", "standard_deviation")
  
  # loop over all columns of "impMat"
  for(fimp_Idx in 1:ncol(impMat)){
    # fill the column "mean" of the matrix "mean_impMat" with the means of one column of "impMat"
    mean_impMat[fimp_Idx, 1] <- mean(impMat[, fimp_Idx])
    # fill the column "standard_deviation" of the matrix "mean_impMat" with the standard deviation of one column of "impMat"
    mean_impMat[fimp_Idx, 2] <- sd(impMat[, fimp_Idx])
  }
  
  # create the folder where everything will be saved
  general_path <- file.path("plots_of_interest")
  dir.create(general_path, showWarnings = FALSE)
  
  analysis_path <- file.path(general_path, save_path)
  dir.create(analysis_path, showWarnings = FALSE)
  
  balance_path <-  file.path(analysis_path, balance)
  dir.create(balance_path, showWarnings = FALSE)
  
  tissue_path <-  file.path(balance_path, tissue_spec)
  dir.create(tissue_path, showWarnings = FALSE)
  
  # define save_path for the "averaged_feature_importance" plot
  save_path_importance <- file.path(tissue_path, plot_name)
  
  # save "importance_plot"
  #svg(file=save_path_importance, width=0.18*ncol(impMat)+0.75)
  pdf(file=save_path_importance, width=0.18*ncol(impMat)+0.75)
  par(mar=c(21, 4, 4, 2)+0.1)
  barplot(sort(mean_impMat[,1]),
          las=2,
          main=drug_name)
  dev.off()
}