getFolders_shrinkToTissue <- function(tissue_spec, drug_res, features, general_path, drug_mat, featureSet, balancing){
  
  # #############################################################################
  # Decide which kind of tissue specific analysis should be done
  # #############################################################################
  
  # setting a boolean parameter "genFolders" to FALSE
  genFolders <- FALSE
  
  # check if tissue spec analysis or not
  if(tissue_spec!="all_tissue"){
    
    # load TISSUE matrix and remove missing cell lines
    load("data/pipeline/cell/TISSUE.ro")
    TISSUE <- TISSUE[match(rownames(features), rownames(TISSUE)),]
    
    # check if there are more than 2*nFold cell-lines for this tissue type
    if (sum(TISSUE[, tissue_spec]==1) < 2*10){
      stop(paste("There are too few cell lines from tissue type ", tissue_spec, 
                 "! Number of cell lines is smaller than 2*10", "!", sep=""))
    }
    else{
      # define a mask containing TRUE and FALSE values depending, if a cell line is responding to this tissue type or not
      tissue_mask <- rownames(TISSUE)[TISSUE[, tissue_spec]==1]
      
      # reduce drug_response to the cell lines, where tissue_mask has TRUE inside
      drug_res <- drug_res[, tissue_mask]
      
      # reduce features to the cell lines, where tissue_mask has TRUE inside
      features <- features[tissue_mask, ]
      
      # set genFolders to TRUE, so that folders will be created later on
      genFolders <- TRUE
    }
  } else {
    # set genFolders to TRUE, so that folders will be created later on
    genFolders <- TRUE
  }
  
  # ##############################################################################################
  # Generate folder structure
  # ##############################################################################################
  
  # only if genFolders equals TRUE the following folders will be created
  if (genFolders) {
    
    # create folder for the drug_response data used
    drug_res_folder <- file.path(general_path, drug_mat)
    dir.create(drug_res_folder, showWarnings = FALSE)
    
    # create tissue type specific folder
    spec_folder <- file.path(drug_res_folder, tissue_spec)
    dir.create(spec_folder, showWarnings = FALSE)
    
    # create balancing folder
    if(balancing==TRUE){
      
      balance_path <-  file.path(spec_folder, "balanced")
      dir.create(balance_path, showWarnings = FALSE)
      
    } else{
      
      balance_path <-  file.path(spec_folder, "unbalanced")
      dir.create(balance_path, showWarnings = FALSE)
    }
    
    # create feature folder
    feature_folder <- file.path(balance_path, featureSet)
    dir.create(feature_folder, showWarnings = FALSE)
  }
  
  return(list(features=features, drug_res=drug_res, feat_folder=feature_folder))
}