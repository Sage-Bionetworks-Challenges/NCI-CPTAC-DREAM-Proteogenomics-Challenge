loadCellData <-  function(featureSet){
  
  # split the commited string for the features by "_"
  split_featureSet <- strsplit(featureSet, "_")
  
  # define an empty vector features which will be filled with feature matrices
  feature <- c()
  
  # check if GEX will be uses
  if ("GEX" %in% split_featureSet[[1]]) {
    
    # set flag to true --> cell_matrix will contain only 978 cells
    gexFlag <- TRUE 
  } else {
    
    gexFlag <- FALSE
  }
  
  
  # loop over all entries for split_featureSet
  for (l_Idx in 1 : length(split_featureSet[[1]])){
    
    #############################################################################################
    # load genFeat
    #############################################################################################
    if(split_featureSet[[1]][l_Idx] == "allGenetics"){ 
      load("data/pipeline/cell/AMPL.ro")
      load("data/pipeline/cell/DEL.ro")
      load("data/pipeline/cell/FUSION.ro")
      FUSION <- na.omit(FUSION)
      load("data/pipeline/cell/miRNA.ro")
      load("data/pipeline/cell/missenseMut.ro")
      load("data/pipeline/cell/truncMut.ro")
      feature <- cbind(feature, t(rbind(AMPL, DEL, FUSION, miRNA, missenseMut, truncMut)))
    }
    
    if(split_featureSet[[1]][l_Idx] == "tissue"){
      load("data/pipeline/cell/TISSUE.ro")
      
      # shrink cell_features to the cells in GEX
      if (gexFlag) {
        
        load("data/pipeline/cell/GEX_avg_2000.ro")  
        GEX <- final
        
        TISSUE <- TISSUE[intersect(rownames(TISSUE), rownames(GEX)),]
      }
      
      feature <- cbind(feature, TISSUE)
    }
    
    if(split_featureSet[[1]][l_Idx] == "missenseMut"){
      load("data/pipeline/cell/missenseMut.ro")
      feature <- cbind(feature, t(missenseMut))
    }
    
    if(split_featureSet[[1]][l_Idx] == "truncMut"){
      load("data/pipeline/cell/truncMut.ro")
      feature <- cbind(feature, t(truncMut))
    }
    
    if(split_featureSet[[1]][l_Idx] == "miRNA"){
      load("data/pipeline/cell/miRNA.ro")
      
      # shrink cell_features to the cells in GEX
      if (gexFlag) {
        
        load("data/pipeline/cell/GEX_avg_2000.ro")  
        GEX <- final
        
        miRNA <- miRNA[, intersect(colnames(miRNA), rownames(GEX))]
      }
      
      feature <- cbind(feature, t(miRNA))
    }
    
    if(split_featureSet[[1]][l_Idx] == "fusion"){
      load("data/pipeline/cell/FUSION.ro")
      FUSION <- na.omit(FUSION)
      
      # shrink cell_features to the cells in GEX
      if (gexFlag) {
        
        load("data/pipeline/cell/GEX_avg_2000.ro")  
        GEX <- final
        
        FUSION <- FUSION[, intersect(colnames(FUSION), rownames(GEX))]
      }
      
      feature <- cbind(feature, t(FUSION))
    }
    
    if(split_featureSet[[1]][l_Idx] == "amplification"){
      load("data/pipeline/cell/AMPL.ro")
      feature <- cbind(feature, t(AMPL))
    }
    
    if(split_featureSet[[1]][l_Idx] == "deletion"){
      load("data/pipeline/cell/DEL.ro")
      feature <- cbind(feature, t(DEL))
    }
    
    if(split_featureSet[[1]][l_Idx] == "CNV"){
      load("data/pipeline/cell/DEL.ro")
      load("data/pipeline/cell/AMPL.ro")
      feat <- rbind(AMPL, DEL)
      
      # shrink cell_features to the cells in GEX
      if (gexFlag) {
        
        load("data/pipeline/cell/GEX_avg_2000.ro")  
        GEX <- final
        
        feat <- feat[, intersect(colnames(feat), rownames(GEX))]
      }
      
      feature <- cbind(feature, t(feat))
    }
    
    if(split_featureSet[[1]][l_Idx] == "mutation"){
      load("data/pipeline/cell/truncMut.ro")
      load("data/pipeline/cell/missenseMut.ro")
      feat <- rbind(missenseMut, truncMut)
      
      # shrink cell_features to the cells in GEX
      if (gexFlag) {
        
        load("data/pipeline/cell/GEX_avg_2000.ro")  
        GEX <- final
        
        feat <- feat[, intersect(colnames(feat), rownames(GEX))]
      }
      
      feature <- cbind(feature, t(feat))
    }
    
    #############################################################################################
    # load methylome
    #############################################################################################
    if(split_featureSet[[1]][l_Idx] == "methylome"){
      load("data/pipeline/cell/METH_avg_2000.ro") 
      METH <- final
      colnames(METH) <- c(paste("Clust_Meth_" , 1:ncol(METH), sep=""))
      
      # shrink cell_features to the cells in GEX
      if (gexFlag) {
        
        load("data/pipeline/cell/GEX_avg_2000.ro")  
        GEX <- final
        
        METH <- METH[intersect(rownames(METH), rownames(GEX)),]
      }
      
      feature <- cbind(feature, METH)
    }
    
    #############################################################################################
    # load gene expression
    #############################################################################################
    if(split_featureSet[[1]][l_Idx] == "GEX"){
      load("data/pipeline/cell/GEX_avg_2000.ro")  
      GEX <- final
      colnames(GEX) <- c(paste("Clust_GEX_" , 1:ncol(GEX), sep=""))
      feature <- cbind(feature, GEX)
    }
  
  }
  
  feature <- na.omit(feature)
  return(feature)
}