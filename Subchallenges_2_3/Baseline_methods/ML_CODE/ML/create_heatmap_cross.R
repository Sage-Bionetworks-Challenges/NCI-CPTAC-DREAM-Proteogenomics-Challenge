#----------------------------------------------------------------------------------------------------------------------------------------
#DESC:  This function creates a heatmap for results been made with crossvalidated data. 
#       The values which will be plotted depend on the committed value of "calc_type". By default 
#       the C-index between the observations and predictions is calculated.
#IN:    stat_model  ==> a string in charge of, which statistical model was used for building the predictions    
#       calc_type   ==> the value of this string decides which type of value will be calculated between the observations and predictions
#                       by default the C-index between observations and predictions
#                       commit "RMSE"     --> the RMSE between observations and predictions will be calculated
#                              "RSQUARE"  --> first the NA entries of the observations and predictions matrix are removed via na.omit
#                                             then the RSQUARE is calculated between those two
#                              "PEARSON"  --> the pearson correlation between observations and predicitions is calcutated giving no
#                                             attendance to NA values
#                              "SPEARMAN" --> the spearman correlation between observations and predicitions is calcutated giving no
#                                             attendance to NA values
#       tissue      ==> a string deciding, which tissue_type was used for the training of a model
#       balance     ==> a string deciding, if the train_set has been balanced for training or not
#       search      ==> a string containing the path to the analysis folder
#       drug_mat    ==> a string containing the drug response matrix
#OUT: returns a heatmap
#----------------------------------------------------------------------------------------------------------------------------------------
plot_heatmap <- function(stat_model, calc_type="cindex", tissue, balance, search, drug_mat){
  
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
    
    # declaring an empty matrix which will be filled with cindexs
    performance <- c()
    
    # declaring an empty vecotr which will be filled with the feature setname
    colnames_performance <- c()
    
    # store the names of all folders, that exist in the directorx located at "search_path"
    folders <- basename(list.dirs(search_path, recursive=FALSE, full.names=FALSE))
    
#     folders <- c("TISSUE:NCI60",
# 				"TISSUE:NCI60_MS:NCI60",
# 				"AMPL:NCI60",
# 				"DEL:NCI60",
# 				"AMPL:NCI60_DEL:NCI60",
# 				"SV:NCI60",
# 				"SAAV:NCI60",
# 				"SV:NCI60_SAAV:NCI60",
# 				"FUSION:NCI60",
# 				"miRNA:NCI60",
# 				"FUSION:NCI60_miRNA:NCI60",
# 				"AMPL:NCI60_DEL:NCI60_FUSION:NCI60_miRNA:NCI60_SV:NCI60_SAAV:NCI60",
# 				"AMPL:NCI60_DEL:NCI60_FUSION:NCI60_miRNA:NCI60_SV:NCI60_SAAV:NCI60_TISSUE:NCI60_MS:NCI60",
# 				"SWATH",
# 				"SHOTGUN",
# 				"SWATH:highConf", 
# 				"SHOTGUN:highConf")
    

    for (f_Idx in 1:length(folders)){
      
      # go into the folder of one feature
      path_to_feature_folder <- file.path(paste(search_path, folders[f_Idx], sep=""))
    
      # load pred
      load(paste(path_to_feature_folder, "/predictions_", stat_model, ".RData", sep=""))

      pred <- predMat
      
      # shrinking the observation matrix to those cell lines, that are in predMat
      obs <- obs[rownames(predMat), ]
      # check if colnames of pred and obs match
      if(sum(colnames(pred) != colnames(obs))){
        stop("colnames of pred and DRUG_RESPONE_RELEASED dont match!")
      }
      
      # check if rownames of pred and obs match
      if(sum(rownames(pred) != rownames(obs))){
        stop("rownames of pred and DRUG_RESPONE_RELEASED dont match!")
      }
      
      if(calc_type=="cindex"){
        
        # define a vector which contains the cindexes for each prediction of one drug 
        perform <- sapply(1:ncol(obs), function(x)cIDX(obs[,x],pred[,x]))
      }
      
      else if(calc_type=="RMSE"){
        
        # define a vector which contains the rmses for each prediction of one drug
        perform <- sapply(1:ncol(obs), function(x)RMSE(obs[,x],pred[,x]))
      }
      
      else if(calc_type=="RSQUARE"){
        
        # define a vector which contains the rquares for each prediction of one drug 
        perform <- sapply(1:ncol(obs), function(x)R2(na.omit(obs[,x]), na.omit(pred[,x])))
      }
      
      else if(calc_type=="PEARSON"){
        
        # define a vector which contains the pearson correlations for each prediction of one drug 
        perform <- sapply(1:ncol(obs), function(x)cor(obs[,x], pred[,x], method="pearson", use ="na.or.complete"))
      }
      
      else if(calc_type=="SPEARMAN"){
        
        # define a vector which contains the pearson correlations for each prediction of one drug 
        perform <- sapply(1:ncol(obs), function(x)cor(obs[,x], pred[,x], method="spearman", use ="na.or.complete"))
      }
      
      # bind the perform to the matrix performance
      performance <- rbind(performance, perform)
      
      # add the actual name of the feature set to a vector 
      colnames_performance <- c(colnames_performance, folders[f_Idx])
    }
    
    # add rownames and colnames to performance
    rownames(performance) <- colnames_performance
    
    # define a vector with all drug names corresponding to all used drug id's
    if(drug_mat=="sanger"){
      
      drugNames <- as.character(DRUG_MASTER_LIST$DRUG_NAME[sapply(colnames(obs), function(x) which(x==DRUG_MASTER_LIST$DRUG_ID))])
    }
    if(drug_mat=="fda"){
      
      fda_trial <- read.csv("data/pipeline/drug/FDA_TRIAL_163.csv",check.names=FALSE,header=TRUE,row.names=1)
      drugNames <- as.character(fda_trial$Symbol[sapply(colnames(obs), function(x) which(x==rownames(fda_trial)))])
    }
    
    if(drug_mat=="moa"){
      
      moa_trial <- read.csv("data/pipeline/drug/MOA_344.csv",check.names=FALSE,header=TRUE,row.names=1)
      drugNames <- as.character(moa_trial$Symbol[sapply(colnames(obs), function(x) which(x==rownames(moa_trial)))])
    }
    colnames(performance) <- paste(colnames(obs),drugNames, sep=":")
    
    # remove NA values from the perfomance matrix
    performance <-performance[ , apply(performance, 2, function(x) !any(is.na(x)))]
    
    # create the folder where everything will be saved
    general_path <- file.path("plots_of_interest")
    dir.create(general_path, showWarnings = FALSE)
    
    analysis_path <- file.path(general_path, search)
    dir.create(analysis_path, showWarnings = FALSE)
    
    drug_mat_path <-  file.path(analysis_path, drug_mat)
    dir.create(drug_mat_path, showWarnings=FALSE)
    
    balance_path <-  file.path(drug_mat_path, balance)
    dir.create(balance_path, showWarnings = FALSE)
    
    tissue_path <-  file.path(balance_path, tissue)
    dir.create(tissue_path, showWarnings = FALSE)
    
    stat_model_path <- file.path(tissue_path, stat_model)
    dir.create(stat_model_path, showWarnings = FALSE)
    
    plot_path <- file.path(stat_model_path, paste(calc_type, ".pdf", sep=""))
    # ##############################################################################
    # Create a heatmap of performance 
    # ##############################################################################
    
    # define the title of the heatmap
    if(tissue!="all_tissue"){
      title_heatmap <- paste(calc_type, "(", tissue, ")", sep="")
    }else{
      title_heatmap <- calc_type
    }
   
    # save heatmap plots
    pdf(file=plot_path, width=55, height=10.625) 
    pheatmap(performance, main=title_heatmap, fontsize_row=15, fontsize_col=15)
    dev.off()
}
