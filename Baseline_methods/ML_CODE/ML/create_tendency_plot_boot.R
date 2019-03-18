#---------------------------------------------------------------------------------------------------------------------------------------
# DESC: This function creates a tendency plot for the feature importance of a certain drug. With this plot you can compare the 
#       the performance of two different feature types and how those features behave as related to the feature TISSUE. 
#       The performance of one feature will be plotted against the one of the other. Round points represent that values, where at least
#       one of the two features outperforms TISSUE, the other points are displayed as rectangles.
#       You can decide, which parameter should be used to estimate the performance. The Wilcoxon pvalue will be added to the plot as
#       well as a boxlpot for each feature type at the side. 
# IN:   - stat_model: a string deciding with which type of statistical model the results have been made
#       - calc_type:  the string deciding the method to calculate performance of the features predicted values
#       - tissue:     the type of tissue specific analysis, that lead to the results
#       - balance:    a string containing, whether the train_set has been balanced or not
#       - search:     a string containing the name of the all sorounding folder for that analysis
#       - feature:    a vector containing the two features, which should be compared and important as a third value the type of TISSUE
#                     matrix to use for comparison
#       - drug_mat:   a string containing the name of the drug response matrix used to built up those results
#---------------------------------------------------------------------------------------------------------------------------------------
tendency <-  function(stat_model, calc_type, tissue, balance, search, feature, drug_mat){
  
  #############################################################################
  # Collect necessary data
  #############################################################################
  
  if(drug_mat=="fda"){
    
    load("data/pipeline/drug/fda_drug_response.ro")
    obs <- t(fda)
  }
  
  if(drug_mat=="moa"){
    
    load("data/pipeline/drug/moa_drug_response.ro")
    obs <- t(moa)
  }
  
  if(drug_mat=="sanger"){
    
    load("data/pipeline/drug/DRUG_RESPONSE.ro")
    load("data/pipeline/drug/DRUG_DESCRIPTION.ro")
    load("data/pipeline/drug/DRUG_MASTER_LIST.ro")
    obs <- t(OMAUC_RELEASED)
  }
  
  # build a path to the folders
  search_path <- file.path(paste(search, "/", drug_mat, "/", tissue, "/", balance, "/", sep=""))
  
  # create pred matrix filled with NA
  pred <- matrix(NA, nrow=length(feature), ncol=ncol(obs))
  rownames(pred) <- c(feature)
  colnames(pred) <- colnames(obs)
  
  for (f_Idx in 1:length(feature)){
    
    # go into the folder of one feature
    path_to_feature_folder <- file.path(paste(search_path, feature[f_Idx], sep=""))
    
    # read all drug folders
    drug_folders <- basename(list.dirs(path_to_feature_folder, recursive=FALSE, full.names=FALSE))
    
    # loop over all drug folders
    for (d_Idx in 1:length(drug_folders)){
      
      # define path to the actual drug
      drug_path <- file.path(path_to_feature_folder, drug_folders[d_Idx])
      
      load(paste(drug_path, "/", stat_model, "/performance_", stat_model, ".RData", sep=""))
      
      # define drug_id corresponding to the name of the drug folder
      drug_id <- strsplit(drug_folders[d_Idx], "drug")[[1]][2]
      
      # define the feature name for the current drug
      feature_name <- feature[f_Idx]
      
      # remove rows with NA values from performance matrix
      performance <- na.omit(performance)
      
      # define the c_index for the current drug 
      if(calc_type=="cindex"){
        
        res <- mean(performance[,1])
      }
      
      # define the pearson cor for the current drug 
      if(calc_type=="PEARSON"){
        
        res <- mean(performance[,4])
      }
      
      # define the spearman cor for the current drug 
      if(calc_type=="SPEARMAN"){
        
        res <- mean(performance[,3])
      }
      
      # define the rmse for the current drug 
      if(calc_type=="RMSE"){
        
        res <- mean(performance[,2])
      }
      
      # add res to the pred matrix
      pred[feature_name, drug_id] <- res
    }
  }
  
  # remove those columns only containing NA values
  pred <- pred[,colSums(is.na(pred)) != nrow(pred)]
  
  # remove those rows only containing NA values
  pred <- pred[rowSums(is.na(pred)) != ncol(pred),]
  
  # remove rows that contain at least one NA value
  pred <- na.omit(pred)
  
  # define the drug_name for the current drug_id
  if(drug_mat=="sanger"){
    
    drugNames <- as.character(DRUG_MASTER_LIST$DRUG_NAME[sapply(colnames(pred), function(x) which(x==DRUG_MASTER_LIST$DRUG_ID))])
  }
  if(drug_mat=="fda"){
    
    fda_trial <- read.csv("data/pipeline/drug/FDA_TRIAL_163.csv",check.names=FALSE,header=TRUE,row.names=1)
    drugNames <- as.character(fda_trial$Symbol[sapply(colnames(pred), function(x) which(x==rownames(fda_trial)))])
  }
  
  if(drug_mat=="moa"){
    
    moa_trial <- read.csv("data/pipeline/drug/MOA_344.csv",check.names=FALSE,header=TRUE,row.names=1)
    drugNames <- as.character(moa_trial$Symbol[sapply(colnames(pred), function(x) which(x==rownames(moa_trial)))])
  }
  
  # defining the column names new to the format: "drug_id:symbol"
  colnames(pred) <- paste(colnames(pred), drugNames, sep=":")
  
  #############################################################################
  # Control the Data
  #############################################################################
  
  # define a vector containing TRUE for the drugs where one or the other feature outperforms TISSUE 
  # and FALSE, if it is the other way around
  perform_mask <- pred[1,]>pred[3,] | pred[2,]>pred[3,]
  
  # vector for the colour encoding containing only "grey"
  colour <- rep("grey", length(perform_mask))
  
  # replace grey by black for the cases, where one or the other feature outperforms TISSUE
  colour[perform_mask] <-  "black"
  
  # vector for the layout of the points, first filled with "0"
  point_layout <- rep(0, length(perform_mask))
  
  # replace the "0" for the layout with "16" if perform_mask is TRUE
  point_layout[perform_mask] <- 16
  
  # shrink matrix to those columns, where one or the other feature outperforms TISSUE
  outperform_pred <- pred[,names(perform_mask[perform_mask==TRUE])]
  
  # calculate Wilcoxon pvalue
  wic_greater <- wilcox.test(outperform_pred[1,], outperform_pred[2,], paired=T, alternative="greater")
  wic_greater <- formatC(wic_greater$p.value, digits=2)
  
  wic_less <- wilcox.test(outperform_pred[1,], outperform_pred[2,], paired=T, alternative="less")
  wic_less <- formatC(wic_less$p.value, digits=2)
  #############################################################################
  # Visualize the Data
  #############################################################################
  
  # the name of the plot how it will be stored
  plot_name <- paste("tendency_plot_", feature[1], "_", feature[2], "_", drug_mat, "_", stat_model, "_", calc_type, ".pdf", sep="")
  
  # create the folder where everything will be saved
  general_path <- file.path("plots_of_interest")
  dir.create(general_path, showWarnings = FALSE)
  
  analysis_path <- file.path(general_path, search)
  dir.create(analysis_path, showWarnings = FALSE)
  
  drug_mat_path <- file.path(analysis_path, drug_mat)
  dir.create(drug_mat_path, showWarnings=FALSE)
  
  balance_path <-  file.path(drug_mat_path, balance)
  dir.create(balance_path, showWarnings = FALSE)
  
  tissue_path <-  file.path(balance_path, tissue)
  dir.create(tissue_path, showWarnings = FALSE)
  
  tendency_path <- file.path(tissue_path, "tendency")
  dir.create(tendency_path, showWarnings = FALSE)
  
  # define save_path for the "averaged_feature_importance" plot
  save_path_importance <- file.path(tendency_path, plot_name)
  
  # open a save function for the plot
  pdf(file=save_path_importance)
  
  # set position of the plot in the grid window
  par(fig=c(0,0.8,0,0.8), new=TRUE)
  
  # plot the two features against each other
  plot(pred[1,], pred[2,], col=colour, pch=point_layout, xlab=rownames(pred)[1], ylab=rownames(pred)[2])
  
  # add a legend to the plot window
  legend("topleft", legend=c(paste("p.value(x<y)= ", wic_less, sep=""), paste("p.value(x>y)= ", wic_greater, sep="")),
         cex=0.9, box.lwd=0, box.col="white")
  legend("bottomright", legend="linear_regression", cex=0.9, fill="red")
  
  # add line of origin to the plot
  abline(0,1)
  
  # add linear regression line to the plot
  abline(lm(outperform_pred[1,]~outperform_pred[2,]), col="red")
  
  par(fig=c(0,0.8,0.55,1), new=TRUE)
  boxplot(pred[1,], horizontal=TRUE, axes=FALSE)
  beeswarm(pred[1,], col="lightgrey", horizontal=TRUE, add=TRUE)
  par(fig=c(0.65,1,0,0.8),new=TRUE)
  boxplot(pred[2,], axes=FALSE)
  beeswarm(pred[2,], col="lightgrey", add=TRUE)
  if(calc_type=="cindex"){
    
    mtext(paste("Tendency_plot for C-Index", sep=""), side=3, outer=TRUE, line=-3) 
  }else{
    
    mtext(paste("Tendency_plot for ", calc_type, sep=""), side=3, outer=TRUE, line=-3) 
    
  }

  dev.off()
}
