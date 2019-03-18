#---------------------------------------------------------------------------------------------------------------------------------------
#DESC:  This function creates a heatmap for results been made with bootstrapped data. 
#       The values which will be plotted depend on the committed value of "calc_type". 
#IN:    stat_model  ==> a string in charge of, which statistical model was used for building the predictions    
#       calc_type   ==> the value of this string decides which type of value calculating will be done
#                       commit "cindex"   --> the c-index between observations and predictions will be calculated 
#                              "RMSE"     --> the RMSE between observations and predictions will be calculated
#                              "PEARSON"  --> the pearson correlation between observations and predicitions is calcutated giving no
#                                             attendance to NA values
#                              "SPEARMAN" --> the spearman correlation between observations and predicitions is calcutated giving no
#                                             attendance to NA values
#       tissue    ==> a string deciding, which tissue_type was used for the training of a model
#       balance   ==> a string deciding, if the train_set has been balanced for training or not
#       search    ==> a strin containing the path to the analysis folder
#       drug_mat  ==> a string containing the drug response matrix
#OUT: returns a heatmap
#---------------------------------------------------------------------------------------------------------------------------------------
plot_heatmap_boot <-  function(stat_model, calc_type, tissue, balance, search, drug_mat){
  
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
    
    if(drug_mat=="sanger_v16"){
      
      load("data/pipeline/drug/DRUG_RESPONSE_v16.ro")
      load("data/pipeline/drug/DRUG_DESCRIPTION.ro")
      load("data/pipeline/drug/DRUG_MASTER_LIST.ro")
      obs <- t(OMAUC_CLASSIFIED)
    }
    # build a path to the folders
    search_path <- file.path(paste(search, "/", drug_mat, "/", tissue, "/", balance, "/", sep=""))
    
    # store the names of all folders, that exist in the directory located at "search_path"
    feat_folders <- basename(list.dirs(search_path, recursive=FALSE, full.names=FALSE))
    
    # create pred matrix filled with NA
    pred <- matrix(NA, nrow=length(feat_folders), ncol=ncol(obs))
    rownames(pred) <- feat_folders
    colnames(pred) <- colnames(obs)
    
    # looping over all elements of "folders"
    for (f_Idx in 1:length(feat_folders)){
      
      # go into the folder of one feature
      path_to_feature_folder <- file.path(paste(search_path, feat_folders[f_Idx], sep=""))
      
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
        feature_name <- feat_folders[f_Idx]
        
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
    
    if(drug_mat=="sanger" | drug_mat=="sanger_v16"){
      
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
    
    
    # remove NA values from the perfomance matrix
#     pred <-pred[ , apply(pred, 2, function(x) !any(is.na(x)))]
    
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
    
    plot_path <- file.path(stat_model_path, paste(calc_type, "_", drug_mat,  ".pdf", sep=""))
    
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
    pheatmap(pred, main=title_heatmap, fontsize_row=15, fontsize_col=15)
    dev.off()
    
}