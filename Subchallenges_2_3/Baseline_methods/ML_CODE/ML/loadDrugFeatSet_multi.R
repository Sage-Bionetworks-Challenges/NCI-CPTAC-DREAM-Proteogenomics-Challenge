loadDrugData <-  function(drugFeatureSet){
  
  # split the commited string for the features by "_"
  split_featureSet <- strsplit(drugFeatureSet, "_")
  
  # define an empty vector features which will be filled with feature matrices
  feature <- c()
  
  # check if bioactivity 2.0 will be uses
  if ("bio2.0" %in% split_featureSet[[1]]) {
    
    # set flag to true --> drug_feature_matriy will contain only 137 drugs
    bioFlag <- TRUE 
  } else {
    
    bioFlag <- FALSE
  }
  
  
  # loop over all entries for split_featureSet
  for (l_Idx in 1 : length(split_featureSet[[1]])){

    #################################################
    # load bioactivity data
    #################################################
    
    if(split_featureSet[[1]][l_Idx] == "bio1.5"){ 
      
      drug_feat <- as.matrix(read.csv("data/pipeline/drug/drug_feat/bioactivity/prediction_1.5_stdev_with_gene_names.csv", 
                                      row.names=1))
      
      # shrink drug_feat to the drugs in bio_2.0
      if (bioFlag) {
        
        bio_2.0 <- as.matrix(read.csv("data/pipeline/drug/drug_feat/bioactivity/prediction_2.0_stdev_with_gene_names.csv", 
                                      row.names=1))
        
        drug_feat <- drug_feat[intersect(rownames(drug_feat),rownames(bio_2.0)),]
      }
      
      feature <- cbind(feature, drug_feat)
      

      
    }
      
    if(split_featureSet[[1]][l_Idx] == "bio2.0"){ 
      
      drug_feat <- as.matrix(read.csv("data/pipeline/drug/drug_feat/bioactivity/prediction_2.0_stdev_with_gene_names.csv", 
                                      row.names=1))
      
      feature <- cbind(feature, drug_feat)
    }
    
    #################################################
    # load chemical properties
    #################################################
    
    if(split_featureSet[[1]][l_Idx] == "chemProp"){ 
      
      drug_feat <- data.matrix(read.csv("data/pipeline/drug/drug_feat/chemProp/drug_properties.csv", 
                                      row.names=1))
      
      # shrink drug_feat to the drugs in bio_2.0
      if (bioFlag) {
        
        bio_2.0 <- as.matrix(read.csv("data/pipeline/drug/drug_feat/bioactivity/prediction_2.0_stdev_with_gene_names.csv", 
                                      row.names=1))
        
        drug_feat <- drug_feat[intersect(rownames(drug_feat),rownames(bio_2.0)),]
      }
      
      feature <- cbind(feature, drug_feat)
      
    }
    
    #################################################
    # load fingerprints data
    #################################################
    
    if(split_featureSet[[1]][l_Idx] == "FP4"){ 
      
      drug_feat <- data.matrix(read.csv("data/pipeline/drug/drug_feat/fingerprints/drug_ECFP_4.csv", 
                                      row.names=1))
      
      # shrink drug_feat to the drugs in bio_2.0
      if (bioFlag) {
        
        bio_2.0 <- as.matrix(read.csv("data/pipeline/drug/drug_feat/bioactivity/prediction_2.0_stdev_with_gene_names.csv", 
                                      row.names=1))
        
        drug_feat <- drug_feat[intersect(rownames(drug_feat),rownames(bio_2.0)),]
      }
      
      feature <- cbind(feature, drug_feat)
      
    }
    
    if(split_featureSet[[1]][l_Idx] == "FP6"){ 
      
      drug_feat <- data.matrix(read.csv("data/pipeline/drug/drug_feat/fingerprints/drug_FCFP_6.csv", 
                                      row.names=1))
      
      # shrink drug_feat to the drugs in bio_2.0
      if (bioFlag) {
        
        bio_2.0 <- as.matrix(read.csv("data/pipeline/drug/drug_feat/bioactivity/prediction_2.0_stdev_with_gene_names.csv", 
                                      row.names=1))
        
        drug_feat <- drug_feat[intersect(rownames(drug_feat),rownames(bio_2.0)),]
      }
      
      feature <- cbind(feature, drug_feat)
      
    }
    
    
    
  }
  
  return(feature)
}
  
  