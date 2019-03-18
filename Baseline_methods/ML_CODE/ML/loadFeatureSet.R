######################################################################################################################################
#DESC : This function loads different feature matrices, returns them and if more than one feature name is committed 
#	it combines the two feature matrices columnwise.
#IN   - feat: a string containing the names of feature matrices, sperated by "_"
#OUT  - feature: matrix, that has been built of all the matrices, according to the input string
######################################################################################################################################
loadData <-  function(feat){
  feature <- c()
  
  #############################################################################################
  # load genFeat
  #############################################################################################
  if(feat=="all_genetics"){ 
    load("data/pipeline/cell/AMPL.ro")
    load("data/pipeline/cell/DEL.ro")
    load("data/pipeline/cell/FUSION.ro")
    load("data/pipeline/cell/miRNA.ro")
    load("data/pipeline/cell/missenseMut.ro")
    load("data/pipeline/cell/truncMut.ro")
    feature <- t(rbind(AMPL, DEL, FUSION, miRNA, missenseMut, truncMut))
  }
  
  if(feat=="tissue"){
    load("data/pipeline/cell/TISSUE.ro")
    feature <- TISSUE
  }
   
  if(feat=="missense_mut"){
    load("data/pipeline/cell/missenseMut.ro")
    feature <- t(missenseMut)
  }
  
  if(feat=="trunc_mut"){
    load("data/pipeline/cell/truncMut.ro")
    feature <- t(truncMut)
  }
  
  if(feat=="miRNA"){
    load("data/pipeline/cell/miRNA.ro")
    feature <- t(miRNA)
  }
  
  if(feat=="fusion"){
    load("data/pipeline/cell/FUSION.ro")
    feature <- t(FUSION)
  }
  
  if(feat=="amplification"){
    load("data/pipeline/cell/AMPL.ro")
    feature <- t(AMPL)
  }
  
  if(feat=="deletion"){
    load("data/pipeline/cell/DEL.ro")
    feature <- t(DEL)
  }
  
  #############################################################################################
  # load methylome
  #############################################################################################
  if(feat=="methylome"){
    load("data/pipeline/cell/METHYLOME.ro")
    feature <- METHYLOME$DATA_ALL
  }
  
  #############################################################################################
  # load gene expression
  #############################################################################################
  if(feat=="GEX"){
    load("data/pipeline/cell/GEXRMA.ro")
    feature <- GEX$RMA_all
  }
      
  #############################################################################################
  # load sanger paper datasets
  #############################################################################################
  if (feat=="v17_BLCA_mut"){
    load("data/pipeline/v17/bem/BLCA_mut.ro")
    commonCell <- rownames(BLCA_mut)
    feature <- cbind(BLCA_mut[commonCell,])
  }
  if (feat=="v17_BLCA_methyl"){
    load("data/pipeline/v17/bem/BLCA_methyl.ro")
    commonCell <- rownames(BLCA_methyl)
    feature <- cbind(BLCA_methyl[commonCell,])
  }
  if (feat=="v17_BLCA_gex"){
    load("data/pipeline/v17/bem/BLCA_gex.ro")
    commonCell <- rownames(BLCA_gex)
    feature <- cbind(BLCA_gex[commonCell,])
  }
  if (feat=="v17_BLCA_cnv"){
    load("data/pipeline/v17/bem/BLCA_cnv.ro")
    commonCell <- rownames(BLCA_cnv)
    feature <- cbind(BLCA_cnv[commonCell,])
  }
  if (feat=="v17_BLCA_speed"){
    load("data/pipeline/v17/bem/BLCA_speed.ro")
    commonCell <- rownames(BLCA_speed)
    feature <- cbind(BLCA_speed[commonCell,])
  }
  if (feat=="v17_BLCA_mut.methyl"){
    load("data/pipeline/v17/bem/BLCA_mut.ro")
    load("data/pipeline/v17/bem/BLCA_methyl.ro")
    commonCell <- intersect(rownames(BLCA_mut), rownames(BLCA_methyl))
    feature <- cbind(BLCA_mut[commonCell,], BLCA_methyl[commonCell,])
  }
  if (feat=="v17_BLCA_mut.gex"){
    load("data/pipeline/v17/bem/BLCA_mut.ro")
    load("data/pipeline/v17/bem/BLCA_gex.ro")
    commonCell <- intersect(rownames(BLCA_mut), rownames(BLCA_gex))
    feature <- cbind(BLCA_mut[commonCell,], BLCA_gex[commonCell,])
  }
  if (feat=="v17_BLCA_mut.cnv"){
    load("data/pipeline/v17/bem/BLCA_mut.ro")
    load("data/pipeline/v17/bem/BLCA_cnv.ro")
    commonCell <- intersect(rownames(BLCA_mut), rownames(BLCA_cnv))
    feature <- cbind(BLCA_mut[commonCell,], BLCA_cnv[commonCell,])
  }
  if (feat=="v17_BLCA_mut.speed"){
    load("data/pipeline/v17/bem/BLCA_mut.ro")
    load("data/pipeline/v17/bem/BLCA_speed.ro")
    commonCell <- intersect(rownames(BLCA_mut), rownames(BLCA_speed))
    feature <- cbind(BLCA_mut[commonCell,], BLCA_speed[commonCell,])
  }
  if (feat=="v17_BLCA_methyl.gex"){
    load("data/pipeline/v17/bem/BLCA_methyl.ro")
    load("data/pipeline/v17/bem/BLCA_gex.ro")
    commonCell <- intersect(rownames(BLCA_methyl), rownames(BLCA_gex))
    feature <- cbind(BLCA_methyl[commonCell,], BLCA_gex[commonCell,])
  }
  if (feat=="v17_BLCA_methyl.cnv"){
    load("data/pipeline/v17/bem/BLCA_methyl.ro")
    load("data/pipeline/v17/bem/BLCA_cnv.ro")
    commonCell <- intersect(rownames(BLCA_methyl), rownames(BLCA_cnv))
    feature <- cbind(BLCA_methyl[commonCell,], BLCA_cnv[commonCell,])
  }
  if (feat=="v17_BLCA_methyl.speed"){
    load("data/pipeline/v17/bem/BLCA_methyl.ro")
    load("data/pipeline/v17/bem/BLCA_speed.ro")
    commonCell <- intersect(rownames(BLCA_methyl), rownames(BLCA_speed))
    feature <- cbind(BLCA_methyl[commonCell,], BLCA_speed[commonCell,])
  }
  if (feat=="v17_BLCA_gex.cnv"){
    load("data/pipeline/v17/bem/BLCA_gex.ro")
    load("data/pipeline/v17/bem/BLCA_cnv.ro")
    commonCell <- intersect(rownames(BLCA_gex), rownames(BLCA_cnv))
    feature <- cbind(BLCA_gex[commonCell,], BLCA_cnv[commonCell,])
  }
  if (feat=="v17_BLCA_gex.speed"){
    load("data/pipeline/v17/bem/BLCA_gex.ro")
    load("data/pipeline/v17/bem/BLCA_speed.ro")
    commonCell <- intersect(rownames(BLCA_gex), rownames(BLCA_speed))
    feature <- cbind(BLCA_gex[commonCell,], BLCA_speed[commonCell,])
  }
  if (feat=="v17_BLCA_cnv.speed"){
    load("data/pipeline/v17/bem/BLCA_cnv.ro")
    load("data/pipeline/v17/bem/BLCA_speed.ro")
    commonCell <- intersect(rownames(BLCA_cnv), rownames(BLCA_speed))
    feature <- cbind(BLCA_cnv[commonCell,], BLCA_speed[commonCell,])
  }
  if (feat=="v17_BLCA_mut.methyl.gex"){
    load("data/pipeline/v17/bem/BLCA_mut.ro")
    load("data/pipeline/v17/bem/BLCA_methyl.ro")
    load("data/pipeline/v17/bem/BLCA_gex.ro")
    commonCell <- intersect(rownames(BLCA_mut), intersect(rownames(BLCA_methyl), rownames(BLCA_gex)))
    feature <- cbind(BLCA_mut[commonCell,], BLCA_methyl[commonCell,], BLCA_gex[commonCell,])
  }
  if (feat=="v17_BLCA_mut.methyl.cnv"){
    load("data/pipeline/v17/bem/BLCA_mut.ro")
    load("data/pipeline/v17/bem/BLCA_methyl.ro")
    load("data/pipeline/v17/bem/BLCA_cnv.ro")
    commonCell <- intersect(rownames(BLCA_mut), intersect(rownames(BLCA_methyl), rownames(BLCA_cnv)))
    feature <- cbind(BLCA_mut[commonCell,], BLCA_methyl[commonCell,], BLCA_cnv[commonCell,])
  }
  if (feat=="v17_BLCA_mut.methyl.speed"){
    load("data/pipeline/v17/bem/BLCA_mut.ro")
    load("data/pipeline/v17/bem/BLCA_methyl.ro")
    load("data/pipeline/v17/bem/BLCA_speed.ro")
    commonCell <- intersect(rownames(BLCA_mut), intersect(rownames(BLCA_methyl), rownames(BLCA_speed)))
    feature <- cbind(BLCA_mut[commonCell,], BLCA_methyl[commonCell,], BLCA_speed[commonCell,])
  }
  if (feat=="v17_BLCA_mut.gex.cnv"){
    load("data/pipeline/v17/bem/BLCA_mut.ro")
    load("data/pipeline/v17/bem/BLCA_gex.ro")
    load("data/pipeline/v17/bem/BLCA_cnv.ro")
    commonCell <- intersect(rownames(BLCA_mut), intersect(rownames(BLCA_gex), rownames(BLCA_cnv)))
    feature <- cbind(BLCA_mut[commonCell,], BLCA_gex[commonCell,], BLCA_cnv[commonCell,])
  }
  if (feat=="v17_BLCA_mut.gex.speed"){
    load("data/pipeline/v17/bem/BLCA_mut.ro")
    load("data/pipeline/v17/bem/BLCA_gex.ro")
    load("data/pipeline/v17/bem/BLCA_speed.ro")
    commonCell <- intersect(rownames(BLCA_mut), intersect(rownames(BLCA_gex), rownames(BLCA_speed)))
    feature <- cbind(BLCA_mut[commonCell,], BLCA_gex[commonCell,], BLCA_speed[commonCell,])
  }
  if (feat=="v17_BLCA_mut.cnv.speed"){
    load("data/pipeline/v17/bem/BLCA_mut.ro")
    load("data/pipeline/v17/bem/BLCA_cnv.ro")
    load("data/pipeline/v17/bem/BLCA_speed.ro")
    commonCell <- intersect(rownames(BLCA_mut), intersect(rownames(BLCA_cnv), rownames(BLCA_speed)))
    feature <- cbind(BLCA_mut[commonCell,], BLCA_cnv[commonCell,], BLCA_speed[commonCell,])
  }
  if (feat=="v17_BLCA_methyl.gex.cnv"){
    load("data/pipeline/v17/bem/BLCA_methyl.ro")
    load("data/pipeline/v17/bem/BLCA_gex.ro")
    load("data/pipeline/v17/bem/BLCA_cnv.ro")
    commonCell <- intersect(rownames(BLCA_methyl), intersect(rownames(BLCA_gex), rownames(BLCA_cnv)))
    feature <- cbind(BLCA_methyl[commonCell,], BLCA_gex[commonCell,], BLCA_cnv[commonCell,])
  }
  if (feat=="v17_BLCA_methyl.gex.speed"){
    load("data/pipeline/v17/bem/BLCA_methyl.ro")
    load("data/pipeline/v17/bem/BLCA_gex.ro")
    load("data/pipeline/v17/bem/BLCA_speed.ro")
    commonCell <- intersect(rownames(BLCA_methyl), intersect(rownames(BLCA_gex), rownames(BLCA_speed)))
    feature <- cbind(BLCA_methyl[commonCell,], BLCA_gex[commonCell,], BLCA_speed[commonCell,])
  }
  if (feat=="v17_BLCA_methyl.cnv.speed"){
    load("data/pipeline/v17/bem/BLCA_methyl.ro")
    load("data/pipeline/v17/bem/BLCA_cnv.ro")
    load("data/pipeline/v17/bem/BLCA_speed.ro")
    commonCell <- intersect(rownames(BLCA_methyl), intersect(rownames(BLCA_cnv), rownames(BLCA_speed)))
    feature <- cbind(BLCA_methyl[commonCell,], BLCA_cnv[commonCell,], BLCA_speed[commonCell,])
  }
  if (feat=="v17_BLCA_gex.cnv.speed"){
    load("data/pipeline/v17/bem/BLCA_gex.ro")
    load("data/pipeline/v17/bem/BLCA_cnv.ro")
    load("data/pipeline/v17/bem/BLCA_speed.ro")
    commonCell <- intersect(rownames(BLCA_gex), intersect(rownames(BLCA_cnv), rownames(BLCA_speed)))
    feature <- cbind(BLCA_gex[commonCell,], BLCA_cnv[commonCell,], BLCA_speed[commonCell,])
  }
  if (feat=="v17_BLCA_mut.methyl.gex.cnv"){
    load("data/pipeline/v17/bem/BLCA_mut.ro")
    load("data/pipeline/v17/bem/BLCA_methyl.ro")
    load("data/pipeline/v17/bem/BLCA_gex.ro")
    load("data/pipeline/v17/bem/BLCA_cnv.ro")
    commonCell <- intersect(rownames(BLCA_mut), intersect(rownames(BLCA_methyl), intersect(rownames(BLCA_gex), rownames(BLCA_cnv))))
    feature <- cbind(BLCA_mut[commonCell,], BLCA_methyl[commonCell,], BLCA_gex[commonCell,], BLCA_cnv[commonCell,])
  }
  if (feat=="v17_BLCA_mut.methyl.gex.speed"){
    load("data/pipeline/v17/bem/BLCA_mut.ro")
    load("data/pipeline/v17/bem/BLCA_methyl.ro")
    load("data/pipeline/v17/bem/BLCA_gex.ro")
    load("data/pipeline/v17/bem/BLCA_speed.ro")
    commonCell <- intersect(rownames(BLCA_mut), intersect(rownames(BLCA_methyl), intersect(rownames(BLCA_gex), rownames(BLCA_speed))))
    feature <- cbind(BLCA_mut[commonCell,], BLCA_methyl[commonCell,], BLCA_gex[commonCell,], BLCA_speed[commonCell,])
  }
  if (feat=="v17_BLCA_mut.methyl.cnv.speed"){
    load("data/pipeline/v17/bem/BLCA_mut.ro")
    load("data/pipeline/v17/bem/BLCA_methyl.ro")
    load("data/pipeline/v17/bem/BLCA_cnv.ro")
    load("data/pipeline/v17/bem/BLCA_speed.ro")
    commonCell <- intersect(rownames(BLCA_mut), intersect(rownames(BLCA_methyl), intersect(rownames(BLCA_cnv), rownames(BLCA_speed))))
    feature <- cbind(BLCA_mut[commonCell,], BLCA_methyl[commonCell,], BLCA_cnv[commonCell,], BLCA_speed[commonCell,])
  }
  if (feat=="v17_BLCA_mut.gex.cnv.speed"){
    load("data/pipeline/v17/bem/BLCA_mut.ro")
    load("data/pipeline/v17/bem/BLCA_gex.ro")
    load("data/pipeline/v17/bem/BLCA_cnv.ro")
    load("data/pipeline/v17/bem/BLCA_speed.ro")
    commonCell <- intersect(rownames(BLCA_mut), intersect(rownames(BLCA_gex), intersect(rownames(BLCA_cnv), rownames(BLCA_speed))))
    feature <- cbind(BLCA_mut[commonCell,], BLCA_gex[commonCell,], BLCA_cnv[commonCell,], BLCA_speed[commonCell,])
  }
  if (feat=="v17_BLCA_methyl.gex.cnv.speed"){
    load("data/pipeline/v17/bem/BLCA_methyl.ro")
    load("data/pipeline/v17/bem/BLCA_gex.ro")
    load("data/pipeline/v17/bem/BLCA_cnv.ro")
    load("data/pipeline/v17/bem/BLCA_speed.ro")
    commonCell <- intersect(rownames(BLCA_methyl), intersect(rownames(BLCA_gex), intersect(rownames(BLCA_cnv), rownames(BLCA_speed))))
    feature <- cbind(BLCA_methyl[commonCell,], BLCA_gex[commonCell,], BLCA_cnv[commonCell,], BLCA_speed[commonCell,])
  }
  if (feat=="v17_BLCA_mut.methyl.gex.cnv.speed"){
    load("data/pipeline/v17/bem/BLCA_mut.ro")
    load("data/pipeline/v17/bem/BLCA_methyl.ro")
    load("data/pipeline/v17/bem/BLCA_gex.ro")
    load("data/pipeline/v17/bem/BLCA_cnv.ro")
    load("data/pipeline/v17/bem/BLCA_speed.ro")
    commonCell <- intersect(rownames(BLCA_mut), intersect(rownames(BLCA_methyl), intersect(rownames(BLCA_gex), intersect(rownames(BLCA_cnv), rownames(BLCA_speed)))))
    feature <- cbind(BLCA_mut[commonCell,], BLCA_methyl[commonCell,], BLCA_gex[commonCell,], BLCA_cnv[commonCell,], BLCA_speed[commonCell,])
  }
  if (feat=="v17_BRCA_mut"){
    load("data/pipeline/v17/bem/BRCA_mut.ro")
    commonCell <- rownames(BRCA_mut)
    feature <- cbind(BRCA_mut[commonCell,])
  }
  if (feat=="v17_BRCA_methyl"){
    load("data/pipeline/v17/bem/BRCA_methyl.ro")
    commonCell <- rownames(BRCA_methyl)
    feature <- cbind(BRCA_methyl[commonCell,])
  }
  if (feat=="v17_BRCA_gex"){
    load("data/pipeline/v17/bem/BRCA_gex.ro")
    commonCell <- rownames(BRCA_gex)
    feature <- cbind(BRCA_gex[commonCell,])
  }
  if (feat=="v17_BRCA_cnv"){
    load("data/pipeline/v17/bem/BRCA_cnv.ro")
    commonCell <- rownames(BRCA_cnv)
    feature <- cbind(BRCA_cnv[commonCell,])
  }
  if (feat=="v17_BRCA_speed"){
    load("data/pipeline/v17/bem/BRCA_speed.ro")
    commonCell <- rownames(BRCA_speed)
    feature <- cbind(BRCA_speed[commonCell,])
  }
  if (feat=="v17_BRCA_mut.methyl"){
    load("data/pipeline/v17/bem/BRCA_mut.ro")
    load("data/pipeline/v17/bem/BRCA_methyl.ro")
    commonCell <- intersect(rownames(BRCA_mut), rownames(BRCA_methyl))
    feature <- cbind(BRCA_mut[commonCell,], BRCA_methyl[commonCell,])
  }
  if (feat=="v17_BRCA_mut.gex"){
    load("data/pipeline/v17/bem/BRCA_mut.ro")
    load("data/pipeline/v17/bem/BRCA_gex.ro")
    commonCell <- intersect(rownames(BRCA_mut), rownames(BRCA_gex))
    feature <- cbind(BRCA_mut[commonCell,], BRCA_gex[commonCell,])
  }
  if (feat=="v17_BRCA_mut.cnv"){
    load("data/pipeline/v17/bem/BRCA_mut.ro")
    load("data/pipeline/v17/bem/BRCA_cnv.ro")
    commonCell <- intersect(rownames(BRCA_mut), rownames(BRCA_cnv))
    feature <- cbind(BRCA_mut[commonCell,], BRCA_cnv[commonCell,])
  }
  if (feat=="v17_BRCA_mut.speed"){
    load("data/pipeline/v17/bem/BRCA_mut.ro")
    load("data/pipeline/v17/bem/BRCA_speed.ro")
    commonCell <- intersect(rownames(BRCA_mut), rownames(BRCA_speed))
    feature <- cbind(BRCA_mut[commonCell,], BRCA_speed[commonCell,])
  }
  if (feat=="v17_BRCA_methyl.gex"){
    load("data/pipeline/v17/bem/BRCA_methyl.ro")
    load("data/pipeline/v17/bem/BRCA_gex.ro")
    commonCell <- intersect(rownames(BRCA_methyl), rownames(BRCA_gex))
    feature <- cbind(BRCA_methyl[commonCell,], BRCA_gex[commonCell,])
  }
  if (feat=="v17_BRCA_methyl.cnv"){
    load("data/pipeline/v17/bem/BRCA_methyl.ro")
    load("data/pipeline/v17/bem/BRCA_cnv.ro")
    commonCell <- intersect(rownames(BRCA_methyl), rownames(BRCA_cnv))
    feature <- cbind(BRCA_methyl[commonCell,], BRCA_cnv[commonCell,])
  }
  if (feat=="v17_BRCA_methyl.speed"){
    load("data/pipeline/v17/bem/BRCA_methyl.ro")
    load("data/pipeline/v17/bem/BRCA_speed.ro")
    commonCell <- intersect(rownames(BRCA_methyl), rownames(BRCA_speed))
    feature <- cbind(BRCA_methyl[commonCell,], BRCA_speed[commonCell,])
  }
  if (feat=="v17_BRCA_gex.cnv"){
    load("data/pipeline/v17/bem/BRCA_gex.ro")
    load("data/pipeline/v17/bem/BRCA_cnv.ro")
    commonCell <- intersect(rownames(BRCA_gex), rownames(BRCA_cnv))
    feature <- cbind(BRCA_gex[commonCell,], BRCA_cnv[commonCell,])
  }
  if (feat=="v17_BRCA_gex.speed"){
    load("data/pipeline/v17/bem/BRCA_gex.ro")
    load("data/pipeline/v17/bem/BRCA_speed.ro")
    commonCell <- intersect(rownames(BRCA_gex), rownames(BRCA_speed))
    feature <- cbind(BRCA_gex[commonCell,], BRCA_speed[commonCell,])
  }
  if (feat=="v17_BRCA_cnv.speed"){
    load("data/pipeline/v17/bem/BRCA_cnv.ro")
    load("data/pipeline/v17/bem/BRCA_speed.ro")
    commonCell <- intersect(rownames(BRCA_cnv), rownames(BRCA_speed))
    feature <- cbind(BRCA_cnv[commonCell,], BRCA_speed[commonCell,])
  }
  if (feat=="v17_BRCA_mut.methyl.gex"){
    load("data/pipeline/v17/bem/BRCA_mut.ro")
    load("data/pipeline/v17/bem/BRCA_methyl.ro")
    load("data/pipeline/v17/bem/BRCA_gex.ro")
    commonCell <- intersect(rownames(BRCA_mut), intersect(rownames(BRCA_methyl), rownames(BRCA_gex)))
    feature <- cbind(BRCA_mut[commonCell,], BRCA_methyl[commonCell,], BRCA_gex[commonCell,])
  }
  if (feat=="v17_BRCA_mut.methyl.cnv"){
    load("data/pipeline/v17/bem/BRCA_mut.ro")
    load("data/pipeline/v17/bem/BRCA_methyl.ro")
    load("data/pipeline/v17/bem/BRCA_cnv.ro")
    commonCell <- intersect(rownames(BRCA_mut), intersect(rownames(BRCA_methyl), rownames(BRCA_cnv)))
    feature <- cbind(BRCA_mut[commonCell,], BRCA_methyl[commonCell,], BRCA_cnv[commonCell,])
  }
  if (feat=="v17_BRCA_mut.methyl.speed"){
    load("data/pipeline/v17/bem/BRCA_mut.ro")
    load("data/pipeline/v17/bem/BRCA_methyl.ro")
    load("data/pipeline/v17/bem/BRCA_speed.ro")
    commonCell <- intersect(rownames(BRCA_mut), intersect(rownames(BRCA_methyl), rownames(BRCA_speed)))
    feature <- cbind(BRCA_mut[commonCell,], BRCA_methyl[commonCell,], BRCA_speed[commonCell,])
  }
  if (feat=="v17_BRCA_mut.gex.cnv"){
    load("data/pipeline/v17/bem/BRCA_mut.ro")
    load("data/pipeline/v17/bem/BRCA_gex.ro")
    load("data/pipeline/v17/bem/BRCA_cnv.ro")
    commonCell <- intersect(rownames(BRCA_mut), intersect(rownames(BRCA_gex), rownames(BRCA_cnv)))
    feature <- cbind(BRCA_mut[commonCell,], BRCA_gex[commonCell,], BRCA_cnv[commonCell,])
  }
  if (feat=="v17_BRCA_mut.gex.speed"){
    load("data/pipeline/v17/bem/BRCA_mut.ro")
    load("data/pipeline/v17/bem/BRCA_gex.ro")
    load("data/pipeline/v17/bem/BRCA_speed.ro")
    commonCell <- intersect(rownames(BRCA_mut), intersect(rownames(BRCA_gex), rownames(BRCA_speed)))
    feature <- cbind(BRCA_mut[commonCell,], BRCA_gex[commonCell,], BRCA_speed[commonCell,])
  }
  if (feat=="v17_BRCA_mut.cnv.speed"){
    load("data/pipeline/v17/bem/BRCA_mut.ro")
    load("data/pipeline/v17/bem/BRCA_cnv.ro")
    load("data/pipeline/v17/bem/BRCA_speed.ro")
    commonCell <- intersect(rownames(BRCA_mut), intersect(rownames(BRCA_cnv), rownames(BRCA_speed)))
    feature <- cbind(BRCA_mut[commonCell,], BRCA_cnv[commonCell,], BRCA_speed[commonCell,])
  }
  if (feat=="v17_BRCA_methyl.gex.cnv"){
    load("data/pipeline/v17/bem/BRCA_methyl.ro")
    load("data/pipeline/v17/bem/BRCA_gex.ro")
    load("data/pipeline/v17/bem/BRCA_cnv.ro")
    commonCell <- intersect(rownames(BRCA_methyl), intersect(rownames(BRCA_gex), rownames(BRCA_cnv)))
    feature <- cbind(BRCA_methyl[commonCell,], BRCA_gex[commonCell,], BRCA_cnv[commonCell,])
  }
  if (feat=="v17_BRCA_methyl.gex.speed"){
    load("data/pipeline/v17/bem/BRCA_methyl.ro")
    load("data/pipeline/v17/bem/BRCA_gex.ro")
    load("data/pipeline/v17/bem/BRCA_speed.ro")
    commonCell <- intersect(rownames(BRCA_methyl), intersect(rownames(BRCA_gex), rownames(BRCA_speed)))
    feature <- cbind(BRCA_methyl[commonCell,], BRCA_gex[commonCell,], BRCA_speed[commonCell,])
  }
  if (feat=="v17_BRCA_methyl.cnv.speed"){
    load("data/pipeline/v17/bem/BRCA_methyl.ro")
    load("data/pipeline/v17/bem/BRCA_cnv.ro")
    load("data/pipeline/v17/bem/BRCA_speed.ro")
    commonCell <- intersect(rownames(BRCA_methyl), intersect(rownames(BRCA_cnv), rownames(BRCA_speed)))
    feature <- cbind(BRCA_methyl[commonCell,], BRCA_cnv[commonCell,], BRCA_speed[commonCell,])
  }
  if (feat=="v17_BRCA_gex.cnv.speed"){
    load("data/pipeline/v17/bem/BRCA_gex.ro")
    load("data/pipeline/v17/bem/BRCA_cnv.ro")
    load("data/pipeline/v17/bem/BRCA_speed.ro")
    commonCell <- intersect(rownames(BRCA_gex), intersect(rownames(BRCA_cnv), rownames(BRCA_speed)))
    feature <- cbind(BRCA_gex[commonCell,], BRCA_cnv[commonCell,], BRCA_speed[commonCell,])
  }
  if (feat=="v17_BRCA_mut.methyl.gex.cnv"){
    load("data/pipeline/v17/bem/BRCA_mut.ro")
    load("data/pipeline/v17/bem/BRCA_methyl.ro")
    load("data/pipeline/v17/bem/BRCA_gex.ro")
    load("data/pipeline/v17/bem/BRCA_cnv.ro")
    commonCell <- intersect(rownames(BRCA_mut), intersect(rownames(BRCA_methyl), intersect(rownames(BRCA_gex), rownames(BRCA_cnv))))
    feature <- cbind(BRCA_mut[commonCell,], BRCA_methyl[commonCell,], BRCA_gex[commonCell,], BRCA_cnv[commonCell,])
  }
  if (feat=="v17_BRCA_mut.methyl.gex.speed"){
    load("data/pipeline/v17/bem/BRCA_mut.ro")
    load("data/pipeline/v17/bem/BRCA_methyl.ro")
    load("data/pipeline/v17/bem/BRCA_gex.ro")
    load("data/pipeline/v17/bem/BRCA_speed.ro")
    commonCell <- intersect(rownames(BRCA_mut), intersect(rownames(BRCA_methyl), intersect(rownames(BRCA_gex), rownames(BRCA_speed))))
    feature <- cbind(BRCA_mut[commonCell,], BRCA_methyl[commonCell,], BRCA_gex[commonCell,], BRCA_speed[commonCell,])
  }
  if (feat=="v17_BRCA_mut.methyl.cnv.speed"){
    load("data/pipeline/v17/bem/BRCA_mut.ro")
    load("data/pipeline/v17/bem/BRCA_methyl.ro")
    load("data/pipeline/v17/bem/BRCA_cnv.ro")
    load("data/pipeline/v17/bem/BRCA_speed.ro")
    commonCell <- intersect(rownames(BRCA_mut), intersect(rownames(BRCA_methyl), intersect(rownames(BRCA_cnv), rownames(BRCA_speed))))
    feature <- cbind(BRCA_mut[commonCell,], BRCA_methyl[commonCell,], BRCA_cnv[commonCell,], BRCA_speed[commonCell,])
  }
  if (feat=="v17_BRCA_mut.gex.cnv.speed"){
    load("data/pipeline/v17/bem/BRCA_mut.ro")
    load("data/pipeline/v17/bem/BRCA_gex.ro")
    load("data/pipeline/v17/bem/BRCA_cnv.ro")
    load("data/pipeline/v17/bem/BRCA_speed.ro")
    commonCell <- intersect(rownames(BRCA_mut), intersect(rownames(BRCA_gex), intersect(rownames(BRCA_cnv), rownames(BRCA_speed))))
    feature <- cbind(BRCA_mut[commonCell,], BRCA_gex[commonCell,], BRCA_cnv[commonCell,], BRCA_speed[commonCell,])
  }
  if (feat=="v17_BRCA_methyl.gex.cnv.speed"){
    load("data/pipeline/v17/bem/BRCA_methyl.ro")
    load("data/pipeline/v17/bem/BRCA_gex.ro")
    load("data/pipeline/v17/bem/BRCA_cnv.ro")
    load("data/pipeline/v17/bem/BRCA_speed.ro")
    commonCell <- intersect(rownames(BRCA_methyl), intersect(rownames(BRCA_gex), intersect(rownames(BRCA_cnv), rownames(BRCA_speed))))
    feature <- cbind(BRCA_methyl[commonCell,], BRCA_gex[commonCell,], BRCA_cnv[commonCell,], BRCA_speed[commonCell,])
  }
  if (feat=="v17_BRCA_mut.methyl.gex.cnv.speed"){
    load("data/pipeline/v17/bem/BRCA_mut.ro")
    load("data/pipeline/v17/bem/BRCA_methyl.ro")
    load("data/pipeline/v17/bem/BRCA_gex.ro")
    load("data/pipeline/v17/bem/BRCA_cnv.ro")
    load("data/pipeline/v17/bem/BRCA_speed.ro")
    commonCell <- intersect(rownames(BRCA_mut), intersect(rownames(BRCA_methyl), intersect(rownames(BRCA_gex), intersect(rownames(BRCA_cnv), rownames(BRCA_speed)))))
    feature <- cbind(BRCA_mut[commonCell,], BRCA_methyl[commonCell,], BRCA_gex[commonCell,], BRCA_cnv[commonCell,], BRCA_speed[commonCell,])
  }
  if (feat=="v17_COREAD_mut"){
    load("data/pipeline/v17/bem/COREAD_mut.ro")
    commonCell <- rownames(COREAD_mut)
    feature <- cbind(COREAD_mut[commonCell,])
  }
  if (feat=="v17_COREAD_methyl"){
    load("data/pipeline/v17/bem/COREAD_methyl.ro")
    commonCell <- rownames(COREAD_methyl)
    feature <- cbind(COREAD_methyl[commonCell,])
  }
  if (feat=="v17_COREAD_gex"){
    load("data/pipeline/v17/bem/COREAD_gex.ro")
    commonCell <- rownames(COREAD_gex)
    feature <- cbind(COREAD_gex[commonCell,])
  }
  if (feat=="v17_COREAD_cnv"){
    load("data/pipeline/v17/bem/COREAD_cnv.ro")
    commonCell <- rownames(COREAD_cnv)
    feature <- cbind(COREAD_cnv[commonCell,])
  }
  if (feat=="v17_COREAD_speed"){
    load("data/pipeline/v17/bem/COREAD_speed.ro")
    commonCell <- rownames(COREAD_speed)
    feature <- cbind(COREAD_speed[commonCell,])
  }
  if (feat=="v17_COREAD_mut.methyl"){
    load("data/pipeline/v17/bem/COREAD_mut.ro")
    load("data/pipeline/v17/bem/COREAD_methyl.ro")
    commonCell <- intersect(rownames(COREAD_mut), rownames(COREAD_methyl))
    feature <- cbind(COREAD_mut[commonCell,], COREAD_methyl[commonCell,])
  }
  if (feat=="v17_COREAD_mut.gex"){
    load("data/pipeline/v17/bem/COREAD_mut.ro")
    load("data/pipeline/v17/bem/COREAD_gex.ro")
    commonCell <- intersect(rownames(COREAD_mut), rownames(COREAD_gex))
    feature <- cbind(COREAD_mut[commonCell,], COREAD_gex[commonCell,])
  }
  if (feat=="v17_COREAD_mut.cnv"){
    load("data/pipeline/v17/bem/COREAD_mut.ro")
    load("data/pipeline/v17/bem/COREAD_cnv.ro")
    commonCell <- intersect(rownames(COREAD_mut), rownames(COREAD_cnv))
    feature <- cbind(COREAD_mut[commonCell,], COREAD_cnv[commonCell,])
  }
  if (feat=="v17_COREAD_mut.speed"){
    load("data/pipeline/v17/bem/COREAD_mut.ro")
    load("data/pipeline/v17/bem/COREAD_speed.ro")
    commonCell <- intersect(rownames(COREAD_mut), rownames(COREAD_speed))
    feature <- cbind(COREAD_mut[commonCell,], COREAD_speed[commonCell,])
  }
  if (feat=="v17_COREAD_methyl.gex"){
    load("data/pipeline/v17/bem/COREAD_methyl.ro")
    load("data/pipeline/v17/bem/COREAD_gex.ro")
    commonCell <- intersect(rownames(COREAD_methyl), rownames(COREAD_gex))
    feature <- cbind(COREAD_methyl[commonCell,], COREAD_gex[commonCell,])
  }
  if (feat=="v17_COREAD_methyl.cnv"){
    load("data/pipeline/v17/bem/COREAD_methyl.ro")
    load("data/pipeline/v17/bem/COREAD_cnv.ro")
    commonCell <- intersect(rownames(COREAD_methyl), rownames(COREAD_cnv))
    feature <- cbind(COREAD_methyl[commonCell,], COREAD_cnv[commonCell,])
  }
  if (feat=="v17_COREAD_methyl.speed"){
    load("data/pipeline/v17/bem/COREAD_methyl.ro")
    load("data/pipeline/v17/bem/COREAD_speed.ro")
    commonCell <- intersect(rownames(COREAD_methyl), rownames(COREAD_speed))
    feature <- cbind(COREAD_methyl[commonCell,], COREAD_speed[commonCell,])
  }
  if (feat=="v17_COREAD_gex.cnv"){
    load("data/pipeline/v17/bem/COREAD_gex.ro")
    load("data/pipeline/v17/bem/COREAD_cnv.ro")
    commonCell <- intersect(rownames(COREAD_gex), rownames(COREAD_cnv))
    feature <- cbind(COREAD_gex[commonCell,], COREAD_cnv[commonCell,])
  }
  if (feat=="v17_COREAD_gex.speed"){
    load("data/pipeline/v17/bem/COREAD_gex.ro")
    load("data/pipeline/v17/bem/COREAD_speed.ro")
    commonCell <- intersect(rownames(COREAD_gex), rownames(COREAD_speed))
    feature <- cbind(COREAD_gex[commonCell,], COREAD_speed[commonCell,])
  }
  if (feat=="v17_COREAD_cnv.speed"){
    load("data/pipeline/v17/bem/COREAD_cnv.ro")
    load("data/pipeline/v17/bem/COREAD_speed.ro")
    commonCell <- intersect(rownames(COREAD_cnv), rownames(COREAD_speed))
    feature <- cbind(COREAD_cnv[commonCell,], COREAD_speed[commonCell,])
  }
  if (feat=="v17_COREAD_mut.methyl.gex"){
    load("data/pipeline/v17/bem/COREAD_mut.ro")
    load("data/pipeline/v17/bem/COREAD_methyl.ro")
    load("data/pipeline/v17/bem/COREAD_gex.ro")
    commonCell <- intersect(rownames(COREAD_mut), intersect(rownames(COREAD_methyl), rownames(COREAD_gex)))
    feature <- cbind(COREAD_mut[commonCell,], COREAD_methyl[commonCell,], COREAD_gex[commonCell,])
  }
  if (feat=="v17_COREAD_mut.methyl.cnv"){
    load("data/pipeline/v17/bem/COREAD_mut.ro")
    load("data/pipeline/v17/bem/COREAD_methyl.ro")
    load("data/pipeline/v17/bem/COREAD_cnv.ro")
    commonCell <- intersect(rownames(COREAD_mut), intersect(rownames(COREAD_methyl), rownames(COREAD_cnv)))
    feature <- cbind(COREAD_mut[commonCell,], COREAD_methyl[commonCell,], COREAD_cnv[commonCell,])
  }
  if (feat=="v17_COREAD_mut.methyl.speed"){
    load("data/pipeline/v17/bem/COREAD_mut.ro")
    load("data/pipeline/v17/bem/COREAD_methyl.ro")
    load("data/pipeline/v17/bem/COREAD_speed.ro")
    commonCell <- intersect(rownames(COREAD_mut), intersect(rownames(COREAD_methyl), rownames(COREAD_speed)))
    feature <- cbind(COREAD_mut[commonCell,], COREAD_methyl[commonCell,], COREAD_speed[commonCell,])
  }
  if (feat=="v17_COREAD_mut.gex.cnv"){
    load("data/pipeline/v17/bem/COREAD_mut.ro")
    load("data/pipeline/v17/bem/COREAD_gex.ro")
    load("data/pipeline/v17/bem/COREAD_cnv.ro")
    commonCell <- intersect(rownames(COREAD_mut), intersect(rownames(COREAD_gex), rownames(COREAD_cnv)))
    feature <- cbind(COREAD_mut[commonCell,], COREAD_gex[commonCell,], COREAD_cnv[commonCell,])
  }
  if (feat=="v17_COREAD_mut.gex.speed"){
    load("data/pipeline/v17/bem/COREAD_mut.ro")
    load("data/pipeline/v17/bem/COREAD_gex.ro")
    load("data/pipeline/v17/bem/COREAD_speed.ro")
    commonCell <- intersect(rownames(COREAD_mut), intersect(rownames(COREAD_gex), rownames(COREAD_speed)))
    feature <- cbind(COREAD_mut[commonCell,], COREAD_gex[commonCell,], COREAD_speed[commonCell,])
  }
  if (feat=="v17_COREAD_mut.cnv.speed"){
    load("data/pipeline/v17/bem/COREAD_mut.ro")
    load("data/pipeline/v17/bem/COREAD_cnv.ro")
    load("data/pipeline/v17/bem/COREAD_speed.ro")
    commonCell <- intersect(rownames(COREAD_mut), intersect(rownames(COREAD_cnv), rownames(COREAD_speed)))
    feature <- cbind(COREAD_mut[commonCell,], COREAD_cnv[commonCell,], COREAD_speed[commonCell,])
  }
  if (feat=="v17_COREAD_methyl.gex.cnv"){
    load("data/pipeline/v17/bem/COREAD_methyl.ro")
    load("data/pipeline/v17/bem/COREAD_gex.ro")
    load("data/pipeline/v17/bem/COREAD_cnv.ro")
    commonCell <- intersect(rownames(COREAD_methyl), intersect(rownames(COREAD_gex), rownames(COREAD_cnv)))
    feature <- cbind(COREAD_methyl[commonCell,], COREAD_gex[commonCell,], COREAD_cnv[commonCell,])
  }
  if (feat=="v17_COREAD_methyl.gex.speed"){
    load("data/pipeline/v17/bem/COREAD_methyl.ro")
    load("data/pipeline/v17/bem/COREAD_gex.ro")
    load("data/pipeline/v17/bem/COREAD_speed.ro")
    commonCell <- intersect(rownames(COREAD_methyl), intersect(rownames(COREAD_gex), rownames(COREAD_speed)))
    feature <- cbind(COREAD_methyl[commonCell,], COREAD_gex[commonCell,], COREAD_speed[commonCell,])
  }
  if (feat=="v17_COREAD_methyl.cnv.speed"){
    load("data/pipeline/v17/bem/COREAD_methyl.ro")
    load("data/pipeline/v17/bem/COREAD_cnv.ro")
    load("data/pipeline/v17/bem/COREAD_speed.ro")
    commonCell <- intersect(rownames(COREAD_methyl), intersect(rownames(COREAD_cnv), rownames(COREAD_speed)))
    feature <- cbind(COREAD_methyl[commonCell,], COREAD_cnv[commonCell,], COREAD_speed[commonCell,])
  }
  if (feat=="v17_COREAD_gex.cnv.speed"){
    load("data/pipeline/v17/bem/COREAD_gex.ro")
    load("data/pipeline/v17/bem/COREAD_cnv.ro")
    load("data/pipeline/v17/bem/COREAD_speed.ro")
    commonCell <- intersect(rownames(COREAD_gex), intersect(rownames(COREAD_cnv), rownames(COREAD_speed)))
    feature <- cbind(COREAD_gex[commonCell,], COREAD_cnv[commonCell,], COREAD_speed[commonCell,])
  }
  if (feat=="v17_COREAD_mut.methyl.gex.cnv"){
    load("data/pipeline/v17/bem/COREAD_mut.ro")
    load("data/pipeline/v17/bem/COREAD_methyl.ro")
    load("data/pipeline/v17/bem/COREAD_gex.ro")
    load("data/pipeline/v17/bem/COREAD_cnv.ro")
    commonCell <- intersect(rownames(COREAD_mut), intersect(rownames(COREAD_methyl), intersect(rownames(COREAD_gex), rownames(COREAD_cnv))))
    feature <- cbind(COREAD_mut[commonCell,], COREAD_methyl[commonCell,], COREAD_gex[commonCell,], COREAD_cnv[commonCell,])
  }
  if (feat=="v17_COREAD_mut.methyl.gex.speed"){
    load("data/pipeline/v17/bem/COREAD_mut.ro")
    load("data/pipeline/v17/bem/COREAD_methyl.ro")
    load("data/pipeline/v17/bem/COREAD_gex.ro")
    load("data/pipeline/v17/bem/COREAD_speed.ro")
    commonCell <- intersect(rownames(COREAD_mut), intersect(rownames(COREAD_methyl), intersect(rownames(COREAD_gex), rownames(COREAD_speed))))
    feature <- cbind(COREAD_mut[commonCell,], COREAD_methyl[commonCell,], COREAD_gex[commonCell,], COREAD_speed[commonCell,])
  }
  if (feat=="v17_COREAD_mut.methyl.cnv.speed"){
    load("data/pipeline/v17/bem/COREAD_mut.ro")
    load("data/pipeline/v17/bem/COREAD_methyl.ro")
    load("data/pipeline/v17/bem/COREAD_cnv.ro")
    load("data/pipeline/v17/bem/COREAD_speed.ro")
    commonCell <- intersect(rownames(COREAD_mut), intersect(rownames(COREAD_methyl), intersect(rownames(COREAD_cnv), rownames(COREAD_speed))))
    feature <- cbind(COREAD_mut[commonCell,], COREAD_methyl[commonCell,], COREAD_cnv[commonCell,], COREAD_speed[commonCell,])
  }
  if (feat=="v17_COREAD_mut.gex.cnv.speed"){
    load("data/pipeline/v17/bem/COREAD_mut.ro")
    load("data/pipeline/v17/bem/COREAD_gex.ro")
    load("data/pipeline/v17/bem/COREAD_cnv.ro")
    load("data/pipeline/v17/bem/COREAD_speed.ro")
    commonCell <- intersect(rownames(COREAD_mut), intersect(rownames(COREAD_gex), intersect(rownames(COREAD_cnv), rownames(COREAD_speed))))
    feature <- cbind(COREAD_mut[commonCell,], COREAD_gex[commonCell,], COREAD_cnv[commonCell,], COREAD_speed[commonCell,])
  }
  if (feat=="v17_COREAD_methyl.gex.cnv.speed"){
    load("data/pipeline/v17/bem/COREAD_methyl.ro")
    load("data/pipeline/v17/bem/COREAD_gex.ro")
    load("data/pipeline/v17/bem/COREAD_cnv.ro")
    load("data/pipeline/v17/bem/COREAD_speed.ro")
    commonCell <- intersect(rownames(COREAD_methyl), intersect(rownames(COREAD_gex), intersect(rownames(COREAD_cnv), rownames(COREAD_speed))))
    feature <- cbind(COREAD_methyl[commonCell,], COREAD_gex[commonCell,], COREAD_cnv[commonCell,], COREAD_speed[commonCell,])
  }
  if (feat=="v17_COREAD_mut.methyl.gex.cnv.speed"){
    load("data/pipeline/v17/bem/COREAD_mut.ro")
    load("data/pipeline/v17/bem/COREAD_methyl.ro")
    load("data/pipeline/v17/bem/COREAD_gex.ro")
    load("data/pipeline/v17/bem/COREAD_cnv.ro")
    load("data/pipeline/v17/bem/COREAD_speed.ro")
    commonCell <- intersect(rownames(COREAD_mut), intersect(rownames(COREAD_methyl), intersect(rownames(COREAD_gex), intersect(rownames(COREAD_cnv), rownames(COREAD_speed)))))
    feature <- cbind(COREAD_mut[commonCell,], COREAD_methyl[commonCell,], COREAD_gex[commonCell,], COREAD_cnv[commonCell,], COREAD_speed[commonCell,])
  }
  if (feat=="v17_DLBC_mut"){
    load("data/pipeline/v17/bem/DLBC_mut.ro")
    commonCell <- rownames(DLBC_mut)
    feature <- cbind(DLBC_mut[commonCell,])
  }
  if (feat=="v17_DLBC_methyl"){
    load("data/pipeline/v17/bem/DLBC_methyl.ro")
    commonCell <- rownames(DLBC_methyl)
    feature <- cbind(DLBC_methyl[commonCell,])
  }
  if (feat=="v17_DLBC_gex"){
    load("data/pipeline/v17/bem/DLBC_gex.ro")
    commonCell <- rownames(DLBC_gex)
    feature <- cbind(DLBC_gex[commonCell,])
  }
  if (feat=="v17_DLBC_cnv"){
    load("data/pipeline/v17/bem/DLBC_cnv.ro")
    commonCell <- rownames(DLBC_cnv)
    feature <- cbind(DLBC_cnv[commonCell,])
  }
  if (feat=="v17_DLBC_speed"){
    load("data/pipeline/v17/bem/DLBC_speed.ro")
    commonCell <- rownames(DLBC_speed)
    feature <- cbind(DLBC_speed[commonCell,])
  }
  if (feat=="v17_DLBC_mut.methyl"){
    load("data/pipeline/v17/bem/DLBC_mut.ro")
    load("data/pipeline/v17/bem/DLBC_methyl.ro")
    commonCell <- intersect(rownames(DLBC_mut), rownames(DLBC_methyl))
    feature <- cbind(DLBC_mut[commonCell,], DLBC_methyl[commonCell,])
  }
  if (feat=="v17_DLBC_mut.gex"){
    load("data/pipeline/v17/bem/DLBC_mut.ro")
    load("data/pipeline/v17/bem/DLBC_gex.ro")
    commonCell <- intersect(rownames(DLBC_mut), rownames(DLBC_gex))
    feature <- cbind(DLBC_mut[commonCell,], DLBC_gex[commonCell,])
  }
  if (feat=="v17_DLBC_mut.cnv"){
    load("data/pipeline/v17/bem/DLBC_mut.ro")
    load("data/pipeline/v17/bem/DLBC_cnv.ro")
    commonCell <- intersect(rownames(DLBC_mut), rownames(DLBC_cnv))
    feature <- cbind(DLBC_mut[commonCell,], DLBC_cnv[commonCell,])
  }
  if (feat=="v17_DLBC_mut.speed"){
    load("data/pipeline/v17/bem/DLBC_mut.ro")
    load("data/pipeline/v17/bem/DLBC_speed.ro")
    commonCell <- intersect(rownames(DLBC_mut), rownames(DLBC_speed))
    feature <- cbind(DLBC_mut[commonCell,], DLBC_speed[commonCell,])
  }
  if (feat=="v17_DLBC_methyl.gex"){
    load("data/pipeline/v17/bem/DLBC_methyl.ro")
    load("data/pipeline/v17/bem/DLBC_gex.ro")
    commonCell <- intersect(rownames(DLBC_methyl), rownames(DLBC_gex))
    feature <- cbind(DLBC_methyl[commonCell,], DLBC_gex[commonCell,])
  }
  if (feat=="v17_DLBC_methyl.cnv"){
    load("data/pipeline/v17/bem/DLBC_methyl.ro")
    load("data/pipeline/v17/bem/DLBC_cnv.ro")
    commonCell <- intersect(rownames(DLBC_methyl), rownames(DLBC_cnv))
    feature <- cbind(DLBC_methyl[commonCell,], DLBC_cnv[commonCell,])
  }
  if (feat=="v17_DLBC_methyl.speed"){
    load("data/pipeline/v17/bem/DLBC_methyl.ro")
    load("data/pipeline/v17/bem/DLBC_speed.ro")
    commonCell <- intersect(rownames(DLBC_methyl), rownames(DLBC_speed))
    feature <- cbind(DLBC_methyl[commonCell,], DLBC_speed[commonCell,])
  }
  if (feat=="v17_DLBC_gex.cnv"){
    load("data/pipeline/v17/bem/DLBC_gex.ro")
    load("data/pipeline/v17/bem/DLBC_cnv.ro")
    commonCell <- intersect(rownames(DLBC_gex), rownames(DLBC_cnv))
    feature <- cbind(DLBC_gex[commonCell,], DLBC_cnv[commonCell,])
  }
  if (feat=="v17_DLBC_gex.speed"){
    load("data/pipeline/v17/bem/DLBC_gex.ro")
    load("data/pipeline/v17/bem/DLBC_speed.ro")
    commonCell <- intersect(rownames(DLBC_gex), rownames(DLBC_speed))
    feature <- cbind(DLBC_gex[commonCell,], DLBC_speed[commonCell,])
  }
  if (feat=="v17_DLBC_cnv.speed"){
    load("data/pipeline/v17/bem/DLBC_cnv.ro")
    load("data/pipeline/v17/bem/DLBC_speed.ro")
    commonCell <- intersect(rownames(DLBC_cnv), rownames(DLBC_speed))
    feature <- cbind(DLBC_cnv[commonCell,], DLBC_speed[commonCell,])
  }
  if (feat=="v17_DLBC_mut.methyl.gex"){
    load("data/pipeline/v17/bem/DLBC_mut.ro")
    load("data/pipeline/v17/bem/DLBC_methyl.ro")
    load("data/pipeline/v17/bem/DLBC_gex.ro")
    commonCell <- intersect(rownames(DLBC_mut), intersect(rownames(DLBC_methyl), rownames(DLBC_gex)))
    feature <- cbind(DLBC_mut[commonCell,], DLBC_methyl[commonCell,], DLBC_gex[commonCell,])
  }
  if (feat=="v17_DLBC_mut.methyl.cnv"){
    load("data/pipeline/v17/bem/DLBC_mut.ro")
    load("data/pipeline/v17/bem/DLBC_methyl.ro")
    load("data/pipeline/v17/bem/DLBC_cnv.ro")
    commonCell <- intersect(rownames(DLBC_mut), intersect(rownames(DLBC_methyl), rownames(DLBC_cnv)))
    feature <- cbind(DLBC_mut[commonCell,], DLBC_methyl[commonCell,], DLBC_cnv[commonCell,])
  }
  if (feat=="v17_DLBC_mut.methyl.speed"){
    load("data/pipeline/v17/bem/DLBC_mut.ro")
    load("data/pipeline/v17/bem/DLBC_methyl.ro")
    load("data/pipeline/v17/bem/DLBC_speed.ro")
    commonCell <- intersect(rownames(DLBC_mut), intersect(rownames(DLBC_methyl), rownames(DLBC_speed)))
    feature <- cbind(DLBC_mut[commonCell,], DLBC_methyl[commonCell,], DLBC_speed[commonCell,])
  }
  if (feat=="v17_DLBC_mut.gex.cnv"){
    load("data/pipeline/v17/bem/DLBC_mut.ro")
    load("data/pipeline/v17/bem/DLBC_gex.ro")
    load("data/pipeline/v17/bem/DLBC_cnv.ro")
    commonCell <- intersect(rownames(DLBC_mut), intersect(rownames(DLBC_gex), rownames(DLBC_cnv)))
    feature <- cbind(DLBC_mut[commonCell,], DLBC_gex[commonCell,], DLBC_cnv[commonCell,])
  }
  if (feat=="v17_DLBC_mut.gex.speed"){
    load("data/pipeline/v17/bem/DLBC_mut.ro")
    load("data/pipeline/v17/bem/DLBC_gex.ro")
    load("data/pipeline/v17/bem/DLBC_speed.ro")
    commonCell <- intersect(rownames(DLBC_mut), intersect(rownames(DLBC_gex), rownames(DLBC_speed)))
    feature <- cbind(DLBC_mut[commonCell,], DLBC_gex[commonCell,], DLBC_speed[commonCell,])
  }
  if (feat=="v17_DLBC_mut.cnv.speed"){
    load("data/pipeline/v17/bem/DLBC_mut.ro")
    load("data/pipeline/v17/bem/DLBC_cnv.ro")
    load("data/pipeline/v17/bem/DLBC_speed.ro")
    commonCell <- intersect(rownames(DLBC_mut), intersect(rownames(DLBC_cnv), rownames(DLBC_speed)))
    feature <- cbind(DLBC_mut[commonCell,], DLBC_cnv[commonCell,], DLBC_speed[commonCell,])
  }
  if (feat=="v17_DLBC_methyl.gex.cnv"){
    load("data/pipeline/v17/bem/DLBC_methyl.ro")
    load("data/pipeline/v17/bem/DLBC_gex.ro")
    load("data/pipeline/v17/bem/DLBC_cnv.ro")
    commonCell <- intersect(rownames(DLBC_methyl), intersect(rownames(DLBC_gex), rownames(DLBC_cnv)))
    feature <- cbind(DLBC_methyl[commonCell,], DLBC_gex[commonCell,], DLBC_cnv[commonCell,])
  }
  if (feat=="v17_DLBC_methyl.gex.speed"){
    load("data/pipeline/v17/bem/DLBC_methyl.ro")
    load("data/pipeline/v17/bem/DLBC_gex.ro")
    load("data/pipeline/v17/bem/DLBC_speed.ro")
    commonCell <- intersect(rownames(DLBC_methyl), intersect(rownames(DLBC_gex), rownames(DLBC_speed)))
    feature <- cbind(DLBC_methyl[commonCell,], DLBC_gex[commonCell,], DLBC_speed[commonCell,])
  }
  if (feat=="v17_DLBC_methyl.cnv.speed"){
    load("data/pipeline/v17/bem/DLBC_methyl.ro")
    load("data/pipeline/v17/bem/DLBC_cnv.ro")
    load("data/pipeline/v17/bem/DLBC_speed.ro")
    commonCell <- intersect(rownames(DLBC_methyl), intersect(rownames(DLBC_cnv), rownames(DLBC_speed)))
    feature <- cbind(DLBC_methyl[commonCell,], DLBC_cnv[commonCell,], DLBC_speed[commonCell,])
  }
  if (feat=="v17_DLBC_gex.cnv.speed"){
    load("data/pipeline/v17/bem/DLBC_gex.ro")
    load("data/pipeline/v17/bem/DLBC_cnv.ro")
    load("data/pipeline/v17/bem/DLBC_speed.ro")
    commonCell <- intersect(rownames(DLBC_gex), intersect(rownames(DLBC_cnv), rownames(DLBC_speed)))
    feature <- cbind(DLBC_gex[commonCell,], DLBC_cnv[commonCell,], DLBC_speed[commonCell,])
  }
  if (feat=="v17_DLBC_mut.methyl.gex.cnv"){
    load("data/pipeline/v17/bem/DLBC_mut.ro")
    load("data/pipeline/v17/bem/DLBC_methyl.ro")
    load("data/pipeline/v17/bem/DLBC_gex.ro")
    load("data/pipeline/v17/bem/DLBC_cnv.ro")
    commonCell <- intersect(rownames(DLBC_mut), intersect(rownames(DLBC_methyl), intersect(rownames(DLBC_gex), rownames(DLBC_cnv))))
    feature <- cbind(DLBC_mut[commonCell,], DLBC_methyl[commonCell,], DLBC_gex[commonCell,], DLBC_cnv[commonCell,])
  }
  if (feat=="v17_DLBC_mut.methyl.gex.speed"){
    load("data/pipeline/v17/bem/DLBC_mut.ro")
    load("data/pipeline/v17/bem/DLBC_methyl.ro")
    load("data/pipeline/v17/bem/DLBC_gex.ro")
    load("data/pipeline/v17/bem/DLBC_speed.ro")
    commonCell <- intersect(rownames(DLBC_mut), intersect(rownames(DLBC_methyl), intersect(rownames(DLBC_gex), rownames(DLBC_speed))))
    feature <- cbind(DLBC_mut[commonCell,], DLBC_methyl[commonCell,], DLBC_gex[commonCell,], DLBC_speed[commonCell,])
  }
  if (feat=="v17_DLBC_mut.methyl.cnv.speed"){
    load("data/pipeline/v17/bem/DLBC_mut.ro")
    load("data/pipeline/v17/bem/DLBC_methyl.ro")
    load("data/pipeline/v17/bem/DLBC_cnv.ro")
    load("data/pipeline/v17/bem/DLBC_speed.ro")
    commonCell <- intersect(rownames(DLBC_mut), intersect(rownames(DLBC_methyl), intersect(rownames(DLBC_cnv), rownames(DLBC_speed))))
    feature <- cbind(DLBC_mut[commonCell,], DLBC_methyl[commonCell,], DLBC_cnv[commonCell,], DLBC_speed[commonCell,])
  }
  if (feat=="v17_DLBC_mut.gex.cnv.speed"){
    load("data/pipeline/v17/bem/DLBC_mut.ro")
    load("data/pipeline/v17/bem/DLBC_gex.ro")
    load("data/pipeline/v17/bem/DLBC_cnv.ro")
    load("data/pipeline/v17/bem/DLBC_speed.ro")
    commonCell <- intersect(rownames(DLBC_mut), intersect(rownames(DLBC_gex), intersect(rownames(DLBC_cnv), rownames(DLBC_speed))))
    feature <- cbind(DLBC_mut[commonCell,], DLBC_gex[commonCell,], DLBC_cnv[commonCell,], DLBC_speed[commonCell,])
  }
  if (feat=="v17_DLBC_methyl.gex.cnv.speed"){
    load("data/pipeline/v17/bem/DLBC_methyl.ro")
    load("data/pipeline/v17/bem/DLBC_gex.ro")
    load("data/pipeline/v17/bem/DLBC_cnv.ro")
    load("data/pipeline/v17/bem/DLBC_speed.ro")
    commonCell <- intersect(rownames(DLBC_methyl), intersect(rownames(DLBC_gex), intersect(rownames(DLBC_cnv), rownames(DLBC_speed))))
    feature <- cbind(DLBC_methyl[commonCell,], DLBC_gex[commonCell,], DLBC_cnv[commonCell,], DLBC_speed[commonCell,])
  }
  if (feat=="v17_DLBC_mut.methyl.gex.cnv.speed"){
    load("data/pipeline/v17/bem/DLBC_mut.ro")
    load("data/pipeline/v17/bem/DLBC_methyl.ro")
    load("data/pipeline/v17/bem/DLBC_gex.ro")
    load("data/pipeline/v17/bem/DLBC_cnv.ro")
    load("data/pipeline/v17/bem/DLBC_speed.ro")
    commonCell <- intersect(rownames(DLBC_mut), intersect(rownames(DLBC_methyl), intersect(rownames(DLBC_gex), intersect(rownames(DLBC_cnv), rownames(DLBC_speed)))))
    feature <- cbind(DLBC_mut[commonCell,], DLBC_methyl[commonCell,], DLBC_gex[commonCell,], DLBC_cnv[commonCell,], DLBC_speed[commonCell,])
  }
  if (feat=="v17_ESCA_mut"){
    load("data/pipeline/v17/bem/ESCA_mut.ro")
    commonCell <- rownames(ESCA_mut)
    feature <- cbind(ESCA_mut[commonCell,])
  }
  if (feat=="v17_ESCA_methyl"){
    load("data/pipeline/v17/bem/ESCA_methyl.ro")
    commonCell <- rownames(ESCA_methyl)
    feature <- cbind(ESCA_methyl[commonCell,])
  }
  if (feat=="v17_ESCA_gex"){
    load("data/pipeline/v17/bem/ESCA_gex.ro")
    commonCell <- rownames(ESCA_gex)
    feature <- cbind(ESCA_gex[commonCell,])
  }
  if (feat=="v17_ESCA_cnv"){
    load("data/pipeline/v17/bem/ESCA_cnv.ro")
    commonCell <- rownames(ESCA_cnv)
    feature <- cbind(ESCA_cnv[commonCell,])
  }
  if (feat=="v17_ESCA_speed"){
    load("data/pipeline/v17/bem/ESCA_speed.ro")
    commonCell <- rownames(ESCA_speed)
    feature <- cbind(ESCA_speed[commonCell,])
  }
  if (feat=="v17_ESCA_mut.methyl"){
    load("data/pipeline/v17/bem/ESCA_mut.ro")
    load("data/pipeline/v17/bem/ESCA_methyl.ro")
    commonCell <- intersect(rownames(ESCA_mut), rownames(ESCA_methyl))
    feature <- cbind(ESCA_mut[commonCell,], ESCA_methyl[commonCell,])
  }
  if (feat=="v17_ESCA_mut.gex"){
    load("data/pipeline/v17/bem/ESCA_mut.ro")
    load("data/pipeline/v17/bem/ESCA_gex.ro")
    commonCell <- intersect(rownames(ESCA_mut), rownames(ESCA_gex))
    feature <- cbind(ESCA_mut[commonCell,], ESCA_gex[commonCell,])
  }
  if (feat=="v17_ESCA_mut.cnv"){
    load("data/pipeline/v17/bem/ESCA_mut.ro")
    load("data/pipeline/v17/bem/ESCA_cnv.ro")
    commonCell <- intersect(rownames(ESCA_mut), rownames(ESCA_cnv))
    feature <- cbind(ESCA_mut[commonCell,], ESCA_cnv[commonCell,])
  }
  if (feat=="v17_ESCA_mut.speed"){
    load("data/pipeline/v17/bem/ESCA_mut.ro")
    load("data/pipeline/v17/bem/ESCA_speed.ro")
    commonCell <- intersect(rownames(ESCA_mut), rownames(ESCA_speed))
    feature <- cbind(ESCA_mut[commonCell,], ESCA_speed[commonCell,])
  }
  if (feat=="v17_ESCA_methyl.gex"){
    load("data/pipeline/v17/bem/ESCA_methyl.ro")
    load("data/pipeline/v17/bem/ESCA_gex.ro")
    commonCell <- intersect(rownames(ESCA_methyl), rownames(ESCA_gex))
    feature <- cbind(ESCA_methyl[commonCell,], ESCA_gex[commonCell,])
  }
  if (feat=="v17_ESCA_methyl.cnv"){
    load("data/pipeline/v17/bem/ESCA_methyl.ro")
    load("data/pipeline/v17/bem/ESCA_cnv.ro")
    commonCell <- intersect(rownames(ESCA_methyl), rownames(ESCA_cnv))
    feature <- cbind(ESCA_methyl[commonCell,], ESCA_cnv[commonCell,])
  }
  if (feat=="v17_ESCA_methyl.speed"){
    load("data/pipeline/v17/bem/ESCA_methyl.ro")
    load("data/pipeline/v17/bem/ESCA_speed.ro")
    commonCell <- intersect(rownames(ESCA_methyl), rownames(ESCA_speed))
    feature <- cbind(ESCA_methyl[commonCell,], ESCA_speed[commonCell,])
  }
  if (feat=="v17_ESCA_gex.cnv"){
    load("data/pipeline/v17/bem/ESCA_gex.ro")
    load("data/pipeline/v17/bem/ESCA_cnv.ro")
    commonCell <- intersect(rownames(ESCA_gex), rownames(ESCA_cnv))
    feature <- cbind(ESCA_gex[commonCell,], ESCA_cnv[commonCell,])
  }
  if (feat=="v17_ESCA_gex.speed"){
    load("data/pipeline/v17/bem/ESCA_gex.ro")
    load("data/pipeline/v17/bem/ESCA_speed.ro")
    commonCell <- intersect(rownames(ESCA_gex), rownames(ESCA_speed))
    feature <- cbind(ESCA_gex[commonCell,], ESCA_speed[commonCell,])
  }
  if (feat=="v17_ESCA_cnv.speed"){
    load("data/pipeline/v17/bem/ESCA_cnv.ro")
    load("data/pipeline/v17/bem/ESCA_speed.ro")
    commonCell <- intersect(rownames(ESCA_cnv), rownames(ESCA_speed))
    feature <- cbind(ESCA_cnv[commonCell,], ESCA_speed[commonCell,])
  }
  if (feat=="v17_ESCA_mut.methyl.gex"){
    load("data/pipeline/v17/bem/ESCA_mut.ro")
    load("data/pipeline/v17/bem/ESCA_methyl.ro")
    load("data/pipeline/v17/bem/ESCA_gex.ro")
    commonCell <- intersect(rownames(ESCA_mut), intersect(rownames(ESCA_methyl), rownames(ESCA_gex)))
    feature <- cbind(ESCA_mut[commonCell,], ESCA_methyl[commonCell,], ESCA_gex[commonCell,])
  }
  if (feat=="v17_ESCA_mut.methyl.cnv"){
    load("data/pipeline/v17/bem/ESCA_mut.ro")
    load("data/pipeline/v17/bem/ESCA_methyl.ro")
    load("data/pipeline/v17/bem/ESCA_cnv.ro")
    commonCell <- intersect(rownames(ESCA_mut), intersect(rownames(ESCA_methyl), rownames(ESCA_cnv)))
    feature <- cbind(ESCA_mut[commonCell,], ESCA_methyl[commonCell,], ESCA_cnv[commonCell,])
  }
  if (feat=="v17_ESCA_mut.methyl.speed"){
    load("data/pipeline/v17/bem/ESCA_mut.ro")
    load("data/pipeline/v17/bem/ESCA_methyl.ro")
    load("data/pipeline/v17/bem/ESCA_speed.ro")
    commonCell <- intersect(rownames(ESCA_mut), intersect(rownames(ESCA_methyl), rownames(ESCA_speed)))
    feature <- cbind(ESCA_mut[commonCell,], ESCA_methyl[commonCell,], ESCA_speed[commonCell,])
  }
  if (feat=="v17_ESCA_mut.gex.cnv"){
    load("data/pipeline/v17/bem/ESCA_mut.ro")
    load("data/pipeline/v17/bem/ESCA_gex.ro")
    load("data/pipeline/v17/bem/ESCA_cnv.ro")
    commonCell <- intersect(rownames(ESCA_mut), intersect(rownames(ESCA_gex), rownames(ESCA_cnv)))
    feature <- cbind(ESCA_mut[commonCell,], ESCA_gex[commonCell,], ESCA_cnv[commonCell,])
  }
  if (feat=="v17_ESCA_mut.gex.speed"){
    load("data/pipeline/v17/bem/ESCA_mut.ro")
    load("data/pipeline/v17/bem/ESCA_gex.ro")
    load("data/pipeline/v17/bem/ESCA_speed.ro")
    commonCell <- intersect(rownames(ESCA_mut), intersect(rownames(ESCA_gex), rownames(ESCA_speed)))
    feature <- cbind(ESCA_mut[commonCell,], ESCA_gex[commonCell,], ESCA_speed[commonCell,])
  }
  if (feat=="v17_ESCA_mut.cnv.speed"){
    load("data/pipeline/v17/bem/ESCA_mut.ro")
    load("data/pipeline/v17/bem/ESCA_cnv.ro")
    load("data/pipeline/v17/bem/ESCA_speed.ro")
    commonCell <- intersect(rownames(ESCA_mut), intersect(rownames(ESCA_cnv), rownames(ESCA_speed)))
    feature <- cbind(ESCA_mut[commonCell,], ESCA_cnv[commonCell,], ESCA_speed[commonCell,])
  }
  if (feat=="v17_ESCA_methyl.gex.cnv"){
    load("data/pipeline/v17/bem/ESCA_methyl.ro")
    load("data/pipeline/v17/bem/ESCA_gex.ro")
    load("data/pipeline/v17/bem/ESCA_cnv.ro")
    commonCell <- intersect(rownames(ESCA_methyl), intersect(rownames(ESCA_gex), rownames(ESCA_cnv)))
    feature <- cbind(ESCA_methyl[commonCell,], ESCA_gex[commonCell,], ESCA_cnv[commonCell,])
  }
  if (feat=="v17_ESCA_methyl.gex.speed"){
    load("data/pipeline/v17/bem/ESCA_methyl.ro")
    load("data/pipeline/v17/bem/ESCA_gex.ro")
    load("data/pipeline/v17/bem/ESCA_speed.ro")
    commonCell <- intersect(rownames(ESCA_methyl), intersect(rownames(ESCA_gex), rownames(ESCA_speed)))
    feature <- cbind(ESCA_methyl[commonCell,], ESCA_gex[commonCell,], ESCA_speed[commonCell,])
  }
  if (feat=="v17_ESCA_methyl.cnv.speed"){
    load("data/pipeline/v17/bem/ESCA_methyl.ro")
    load("data/pipeline/v17/bem/ESCA_cnv.ro")
    load("data/pipeline/v17/bem/ESCA_speed.ro")
    commonCell <- intersect(rownames(ESCA_methyl), intersect(rownames(ESCA_cnv), rownames(ESCA_speed)))
    feature <- cbind(ESCA_methyl[commonCell,], ESCA_cnv[commonCell,], ESCA_speed[commonCell,])
  }
  if (feat=="v17_ESCA_gex.cnv.speed"){
    load("data/pipeline/v17/bem/ESCA_gex.ro")
    load("data/pipeline/v17/bem/ESCA_cnv.ro")
    load("data/pipeline/v17/bem/ESCA_speed.ro")
    commonCell <- intersect(rownames(ESCA_gex), intersect(rownames(ESCA_cnv), rownames(ESCA_speed)))
    feature <- cbind(ESCA_gex[commonCell,], ESCA_cnv[commonCell,], ESCA_speed[commonCell,])
  }
  if (feat=="v17_ESCA_mut.methyl.gex.cnv"){
    load("data/pipeline/v17/bem/ESCA_mut.ro")
    load("data/pipeline/v17/bem/ESCA_methyl.ro")
    load("data/pipeline/v17/bem/ESCA_gex.ro")
    load("data/pipeline/v17/bem/ESCA_cnv.ro")
    commonCell <- intersect(rownames(ESCA_mut), intersect(rownames(ESCA_methyl), intersect(rownames(ESCA_gex), rownames(ESCA_cnv))))
    feature <- cbind(ESCA_mut[commonCell,], ESCA_methyl[commonCell,], ESCA_gex[commonCell,], ESCA_cnv[commonCell,])
  }
  if (feat=="v17_ESCA_mut.methyl.gex.speed"){
    load("data/pipeline/v17/bem/ESCA_mut.ro")
    load("data/pipeline/v17/bem/ESCA_methyl.ro")
    load("data/pipeline/v17/bem/ESCA_gex.ro")
    load("data/pipeline/v17/bem/ESCA_speed.ro")
    commonCell <- intersect(rownames(ESCA_mut), intersect(rownames(ESCA_methyl), intersect(rownames(ESCA_gex), rownames(ESCA_speed))))
    feature <- cbind(ESCA_mut[commonCell,], ESCA_methyl[commonCell,], ESCA_gex[commonCell,], ESCA_speed[commonCell,])
  }
  if (feat=="v17_ESCA_mut.methyl.cnv.speed"){
    load("data/pipeline/v17/bem/ESCA_mut.ro")
    load("data/pipeline/v17/bem/ESCA_methyl.ro")
    load("data/pipeline/v17/bem/ESCA_cnv.ro")
    load("data/pipeline/v17/bem/ESCA_speed.ro")
    commonCell <- intersect(rownames(ESCA_mut), intersect(rownames(ESCA_methyl), intersect(rownames(ESCA_cnv), rownames(ESCA_speed))))
    feature <- cbind(ESCA_mut[commonCell,], ESCA_methyl[commonCell,], ESCA_cnv[commonCell,], ESCA_speed[commonCell,])
  }
  if (feat=="v17_ESCA_mut.gex.cnv.speed"){
    load("data/pipeline/v17/bem/ESCA_mut.ro")
    load("data/pipeline/v17/bem/ESCA_gex.ro")
    load("data/pipeline/v17/bem/ESCA_cnv.ro")
    load("data/pipeline/v17/bem/ESCA_speed.ro")
    commonCell <- intersect(rownames(ESCA_mut), intersect(rownames(ESCA_gex), intersect(rownames(ESCA_cnv), rownames(ESCA_speed))))
    feature <- cbind(ESCA_mut[commonCell,], ESCA_gex[commonCell,], ESCA_cnv[commonCell,], ESCA_speed[commonCell,])
  }
  if (feat=="v17_ESCA_methyl.gex.cnv.speed"){
    load("data/pipeline/v17/bem/ESCA_methyl.ro")
    load("data/pipeline/v17/bem/ESCA_gex.ro")
    load("data/pipeline/v17/bem/ESCA_cnv.ro")
    load("data/pipeline/v17/bem/ESCA_speed.ro")
    commonCell <- intersect(rownames(ESCA_methyl), intersect(rownames(ESCA_gex), intersect(rownames(ESCA_cnv), rownames(ESCA_speed))))
    feature <- cbind(ESCA_methyl[commonCell,], ESCA_gex[commonCell,], ESCA_cnv[commonCell,], ESCA_speed[commonCell,])
  }
  if (feat=="v17_ESCA_mut.methyl.gex.cnv.speed"){
    load("data/pipeline/v17/bem/ESCA_mut.ro")
    load("data/pipeline/v17/bem/ESCA_methyl.ro")
    load("data/pipeline/v17/bem/ESCA_gex.ro")
    load("data/pipeline/v17/bem/ESCA_cnv.ro")
    load("data/pipeline/v17/bem/ESCA_speed.ro")
    commonCell <- intersect(rownames(ESCA_mut), intersect(rownames(ESCA_methyl), intersect(rownames(ESCA_gex), intersect(rownames(ESCA_cnv), rownames(ESCA_speed)))))
    feature <- cbind(ESCA_mut[commonCell,], ESCA_methyl[commonCell,], ESCA_gex[commonCell,], ESCA_cnv[commonCell,], ESCA_speed[commonCell,])
  }
  if (feat=="v17_GBM_mut"){
    load("data/pipeline/v17/bem/GBM_mut.ro")
    commonCell <- rownames(GBM_mut)
    feature <- cbind(GBM_mut[commonCell,])
  }
  if (feat=="v17_GBM_methyl"){
    load("data/pipeline/v17/bem/GBM_methyl.ro")
    commonCell <- rownames(GBM_methyl)
    feature <- cbind(GBM_methyl[commonCell,])
  }
  if (feat=="v17_GBM_gex"){
    load("data/pipeline/v17/bem/GBM_gex.ro")
    commonCell <- rownames(GBM_gex)
    feature <- cbind(GBM_gex[commonCell,])
  }
  if (feat=="v17_GBM_cnv"){
    load("data/pipeline/v17/bem/GBM_cnv.ro")
    commonCell <- rownames(GBM_cnv)
    feature <- cbind(GBM_cnv[commonCell,])
  }
  if (feat=="v17_GBM_speed"){
    load("data/pipeline/v17/bem/GBM_speed.ro")
    commonCell <- rownames(GBM_speed)
    feature <- cbind(GBM_speed[commonCell,])
  }
  if (feat=="v17_GBM_mut.methyl"){
    load("data/pipeline/v17/bem/GBM_mut.ro")
    load("data/pipeline/v17/bem/GBM_methyl.ro")
    commonCell <- intersect(rownames(GBM_mut), rownames(GBM_methyl))
    feature <- cbind(GBM_mut[commonCell,], GBM_methyl[commonCell,])
  }
  if (feat=="v17_GBM_mut.gex"){
    load("data/pipeline/v17/bem/GBM_mut.ro")
    load("data/pipeline/v17/bem/GBM_gex.ro")
    commonCell <- intersect(rownames(GBM_mut), rownames(GBM_gex))
    feature <- cbind(GBM_mut[commonCell,], GBM_gex[commonCell,])
  }
  if (feat=="v17_GBM_mut.cnv"){
    load("data/pipeline/v17/bem/GBM_mut.ro")
    load("data/pipeline/v17/bem/GBM_cnv.ro")
    commonCell <- intersect(rownames(GBM_mut), rownames(GBM_cnv))
    feature <- cbind(GBM_mut[commonCell,], GBM_cnv[commonCell,])
  }
  if (feat=="v17_GBM_mut.speed"){
    load("data/pipeline/v17/bem/GBM_mut.ro")
    load("data/pipeline/v17/bem/GBM_speed.ro")
    commonCell <- intersect(rownames(GBM_mut), rownames(GBM_speed))
    feature <- cbind(GBM_mut[commonCell,], GBM_speed[commonCell,])
  }
  if (feat=="v17_GBM_methyl.gex"){
    load("data/pipeline/v17/bem/GBM_methyl.ro")
    load("data/pipeline/v17/bem/GBM_gex.ro")
    commonCell <- intersect(rownames(GBM_methyl), rownames(GBM_gex))
    feature <- cbind(GBM_methyl[commonCell,], GBM_gex[commonCell,])
  }
  if (feat=="v17_GBM_methyl.cnv"){
    load("data/pipeline/v17/bem/GBM_methyl.ro")
    load("data/pipeline/v17/bem/GBM_cnv.ro")
    commonCell <- intersect(rownames(GBM_methyl), rownames(GBM_cnv))
    feature <- cbind(GBM_methyl[commonCell,], GBM_cnv[commonCell,])
  }
  if (feat=="v17_GBM_methyl.speed"){
    load("data/pipeline/v17/bem/GBM_methyl.ro")
    load("data/pipeline/v17/bem/GBM_speed.ro")
    commonCell <- intersect(rownames(GBM_methyl), rownames(GBM_speed))
    feature <- cbind(GBM_methyl[commonCell,], GBM_speed[commonCell,])
  }
  if (feat=="v17_GBM_gex.cnv"){
    load("data/pipeline/v17/bem/GBM_gex.ro")
    load("data/pipeline/v17/bem/GBM_cnv.ro")
    commonCell <- intersect(rownames(GBM_gex), rownames(GBM_cnv))
    feature <- cbind(GBM_gex[commonCell,], GBM_cnv[commonCell,])
  }
  if (feat=="v17_GBM_gex.speed"){
    load("data/pipeline/v17/bem/GBM_gex.ro")
    load("data/pipeline/v17/bem/GBM_speed.ro")
    commonCell <- intersect(rownames(GBM_gex), rownames(GBM_speed))
    feature <- cbind(GBM_gex[commonCell,], GBM_speed[commonCell,])
  }
  if (feat=="v17_GBM_cnv.speed"){
    load("data/pipeline/v17/bem/GBM_cnv.ro")
    load("data/pipeline/v17/bem/GBM_speed.ro")
    commonCell <- intersect(rownames(GBM_cnv), rownames(GBM_speed))
    feature <- cbind(GBM_cnv[commonCell,], GBM_speed[commonCell,])
  }
  if (feat=="v17_GBM_mut.methyl.gex"){
    load("data/pipeline/v17/bem/GBM_mut.ro")
    load("data/pipeline/v17/bem/GBM_methyl.ro")
    load("data/pipeline/v17/bem/GBM_gex.ro")
    commonCell <- intersect(rownames(GBM_mut), intersect(rownames(GBM_methyl), rownames(GBM_gex)))
    feature <- cbind(GBM_mut[commonCell,], GBM_methyl[commonCell,], GBM_gex[commonCell,])
  }
  if (feat=="v17_GBM_mut.methyl.cnv"){
    load("data/pipeline/v17/bem/GBM_mut.ro")
    load("data/pipeline/v17/bem/GBM_methyl.ro")
    load("data/pipeline/v17/bem/GBM_cnv.ro")
    commonCell <- intersect(rownames(GBM_mut), intersect(rownames(GBM_methyl), rownames(GBM_cnv)))
    feature <- cbind(GBM_mut[commonCell,], GBM_methyl[commonCell,], GBM_cnv[commonCell,])
  }
  if (feat=="v17_GBM_mut.methyl.speed"){
    load("data/pipeline/v17/bem/GBM_mut.ro")
    load("data/pipeline/v17/bem/GBM_methyl.ro")
    load("data/pipeline/v17/bem/GBM_speed.ro")
    commonCell <- intersect(rownames(GBM_mut), intersect(rownames(GBM_methyl), rownames(GBM_speed)))
    feature <- cbind(GBM_mut[commonCell,], GBM_methyl[commonCell,], GBM_speed[commonCell,])
  }
  if (feat=="v17_GBM_mut.gex.cnv"){
    load("data/pipeline/v17/bem/GBM_mut.ro")
    load("data/pipeline/v17/bem/GBM_gex.ro")
    load("data/pipeline/v17/bem/GBM_cnv.ro")
    commonCell <- intersect(rownames(GBM_mut), intersect(rownames(GBM_gex), rownames(GBM_cnv)))
    feature <- cbind(GBM_mut[commonCell,], GBM_gex[commonCell,], GBM_cnv[commonCell,])
  }
  if (feat=="v17_GBM_mut.gex.speed"){
    load("data/pipeline/v17/bem/GBM_mut.ro")
    load("data/pipeline/v17/bem/GBM_gex.ro")
    load("data/pipeline/v17/bem/GBM_speed.ro")
    commonCell <- intersect(rownames(GBM_mut), intersect(rownames(GBM_gex), rownames(GBM_speed)))
    feature <- cbind(GBM_mut[commonCell,], GBM_gex[commonCell,], GBM_speed[commonCell,])
  }
  if (feat=="v17_GBM_mut.cnv.speed"){
    load("data/pipeline/v17/bem/GBM_mut.ro")
    load("data/pipeline/v17/bem/GBM_cnv.ro")
    load("data/pipeline/v17/bem/GBM_speed.ro")
    commonCell <- intersect(rownames(GBM_mut), intersect(rownames(GBM_cnv), rownames(GBM_speed)))
    feature <- cbind(GBM_mut[commonCell,], GBM_cnv[commonCell,], GBM_speed[commonCell,])
  }
  if (feat=="v17_GBM_methyl.gex.cnv"){
    load("data/pipeline/v17/bem/GBM_methyl.ro")
    load("data/pipeline/v17/bem/GBM_gex.ro")
    load("data/pipeline/v17/bem/GBM_cnv.ro")
    commonCell <- intersect(rownames(GBM_methyl), intersect(rownames(GBM_gex), rownames(GBM_cnv)))
    feature <- cbind(GBM_methyl[commonCell,], GBM_gex[commonCell,], GBM_cnv[commonCell,])
  }
  if (feat=="v17_GBM_methyl.gex.speed"){
    load("data/pipeline/v17/bem/GBM_methyl.ro")
    load("data/pipeline/v17/bem/GBM_gex.ro")
    load("data/pipeline/v17/bem/GBM_speed.ro")
    commonCell <- intersect(rownames(GBM_methyl), intersect(rownames(GBM_gex), rownames(GBM_speed)))
    feature <- cbind(GBM_methyl[commonCell,], GBM_gex[commonCell,], GBM_speed[commonCell,])
  }
  if (feat=="v17_GBM_methyl.cnv.speed"){
    load("data/pipeline/v17/bem/GBM_methyl.ro")
    load("data/pipeline/v17/bem/GBM_cnv.ro")
    load("data/pipeline/v17/bem/GBM_speed.ro")
    commonCell <- intersect(rownames(GBM_methyl), intersect(rownames(GBM_cnv), rownames(GBM_speed)))
    feature <- cbind(GBM_methyl[commonCell,], GBM_cnv[commonCell,], GBM_speed[commonCell,])
  }
  if (feat=="v17_GBM_gex.cnv.speed"){
    load("data/pipeline/v17/bem/GBM_gex.ro")
    load("data/pipeline/v17/bem/GBM_cnv.ro")
    load("data/pipeline/v17/bem/GBM_speed.ro")
    commonCell <- intersect(rownames(GBM_gex), intersect(rownames(GBM_cnv), rownames(GBM_speed)))
    feature <- cbind(GBM_gex[commonCell,], GBM_cnv[commonCell,], GBM_speed[commonCell,])
  }
  if (feat=="v17_GBM_mut.methyl.gex.cnv"){
    load("data/pipeline/v17/bem/GBM_mut.ro")
    load("data/pipeline/v17/bem/GBM_methyl.ro")
    load("data/pipeline/v17/bem/GBM_gex.ro")
    load("data/pipeline/v17/bem/GBM_cnv.ro")
    commonCell <- intersect(rownames(GBM_mut), intersect(rownames(GBM_methyl), intersect(rownames(GBM_gex), rownames(GBM_cnv))))
    feature <- cbind(GBM_mut[commonCell,], GBM_methyl[commonCell,], GBM_gex[commonCell,], GBM_cnv[commonCell,])
  }
  if (feat=="v17_GBM_mut.methyl.gex.speed"){
    load("data/pipeline/v17/bem/GBM_mut.ro")
    load("data/pipeline/v17/bem/GBM_methyl.ro")
    load("data/pipeline/v17/bem/GBM_gex.ro")
    load("data/pipeline/v17/bem/GBM_speed.ro")
    commonCell <- intersect(rownames(GBM_mut), intersect(rownames(GBM_methyl), intersect(rownames(GBM_gex), rownames(GBM_speed))))
    feature <- cbind(GBM_mut[commonCell,], GBM_methyl[commonCell,], GBM_gex[commonCell,], GBM_speed[commonCell,])
  }
  if (feat=="v17_GBM_mut.methyl.cnv.speed"){
    load("data/pipeline/v17/bem/GBM_mut.ro")
    load("data/pipeline/v17/bem/GBM_methyl.ro")
    load("data/pipeline/v17/bem/GBM_cnv.ro")
    load("data/pipeline/v17/bem/GBM_speed.ro")
    commonCell <- intersect(rownames(GBM_mut), intersect(rownames(GBM_methyl), intersect(rownames(GBM_cnv), rownames(GBM_speed))))
    feature <- cbind(GBM_mut[commonCell,], GBM_methyl[commonCell,], GBM_cnv[commonCell,], GBM_speed[commonCell,])
  }
  if (feat=="v17_GBM_mut.gex.cnv.speed"){
    load("data/pipeline/v17/bem/GBM_mut.ro")
    load("data/pipeline/v17/bem/GBM_gex.ro")
    load("data/pipeline/v17/bem/GBM_cnv.ro")
    load("data/pipeline/v17/bem/GBM_speed.ro")
    commonCell <- intersect(rownames(GBM_mut), intersect(rownames(GBM_gex), intersect(rownames(GBM_cnv), rownames(GBM_speed))))
    feature <- cbind(GBM_mut[commonCell,], GBM_gex[commonCell,], GBM_cnv[commonCell,], GBM_speed[commonCell,])
  }
  if (feat=="v17_GBM_methyl.gex.cnv.speed"){
    load("data/pipeline/v17/bem/GBM_methyl.ro")
    load("data/pipeline/v17/bem/GBM_gex.ro")
    load("data/pipeline/v17/bem/GBM_cnv.ro")
    load("data/pipeline/v17/bem/GBM_speed.ro")
    commonCell <- intersect(rownames(GBM_methyl), intersect(rownames(GBM_gex), intersect(rownames(GBM_cnv), rownames(GBM_speed))))
    feature <- cbind(GBM_methyl[commonCell,], GBM_gex[commonCell,], GBM_cnv[commonCell,], GBM_speed[commonCell,])
  }
  if (feat=="v17_GBM_mut.methyl.gex.cnv.speed"){
    load("data/pipeline/v17/bem/GBM_mut.ro")
    load("data/pipeline/v17/bem/GBM_methyl.ro")
    load("data/pipeline/v17/bem/GBM_gex.ro")
    load("data/pipeline/v17/bem/GBM_cnv.ro")
    load("data/pipeline/v17/bem/GBM_speed.ro")
    commonCell <- intersect(rownames(GBM_mut), intersect(rownames(GBM_methyl), intersect(rownames(GBM_gex), intersect(rownames(GBM_cnv), rownames(GBM_speed)))))
    feature <- cbind(GBM_mut[commonCell,], GBM_methyl[commonCell,], GBM_gex[commonCell,], GBM_cnv[commonCell,], GBM_speed[commonCell,])
  }
  if (feat=="v17_HNSC_mut"){
    load("data/pipeline/v17/bem/HNSC_mut.ro")
    commonCell <- rownames(HNSC_mut)
    feature <- cbind(HNSC_mut[commonCell,])
  }
  if (feat=="v17_HNSC_methyl"){
    load("data/pipeline/v17/bem/HNSC_methyl.ro")
    commonCell <- rownames(HNSC_methyl)
    feature <- cbind(HNSC_methyl[commonCell,])
  }
  if (feat=="v17_HNSC_gex"){
    load("data/pipeline/v17/bem/HNSC_gex.ro")
    commonCell <- rownames(HNSC_gex)
    feature <- cbind(HNSC_gex[commonCell,])
  }
  if (feat=="v17_HNSC_cnv"){
    load("data/pipeline/v17/bem/HNSC_cnv.ro")
    commonCell <- rownames(HNSC_cnv)
    feature <- cbind(HNSC_cnv[commonCell,])
  }
  if (feat=="v17_HNSC_speed"){
    load("data/pipeline/v17/bem/HNSC_speed.ro")
    commonCell <- rownames(HNSC_speed)
    feature <- cbind(HNSC_speed[commonCell,])
  }
  if (feat=="v17_HNSC_mut.methyl"){
    load("data/pipeline/v17/bem/HNSC_mut.ro")
    load("data/pipeline/v17/bem/HNSC_methyl.ro")
    commonCell <- intersect(rownames(HNSC_mut), rownames(HNSC_methyl))
    feature <- cbind(HNSC_mut[commonCell,], HNSC_methyl[commonCell,])
  }
  if (feat=="v17_HNSC_mut.gex"){
    load("data/pipeline/v17/bem/HNSC_mut.ro")
    load("data/pipeline/v17/bem/HNSC_gex.ro")
    commonCell <- intersect(rownames(HNSC_mut), rownames(HNSC_gex))
    feature <- cbind(HNSC_mut[commonCell,], HNSC_gex[commonCell,])
  }
  if (feat=="v17_HNSC_mut.cnv"){
    load("data/pipeline/v17/bem/HNSC_mut.ro")
    load("data/pipeline/v17/bem/HNSC_cnv.ro")
    commonCell <- intersect(rownames(HNSC_mut), rownames(HNSC_cnv))
    feature <- cbind(HNSC_mut[commonCell,], HNSC_cnv[commonCell,])
  }
  if (feat=="v17_HNSC_mut.speed"){
    load("data/pipeline/v17/bem/HNSC_mut.ro")
    load("data/pipeline/v17/bem/HNSC_speed.ro")
    commonCell <- intersect(rownames(HNSC_mut), rownames(HNSC_speed))
    feature <- cbind(HNSC_mut[commonCell,], HNSC_speed[commonCell,])
  }
  if (feat=="v17_HNSC_methyl.gex"){
    load("data/pipeline/v17/bem/HNSC_methyl.ro")
    load("data/pipeline/v17/bem/HNSC_gex.ro")
    commonCell <- intersect(rownames(HNSC_methyl), rownames(HNSC_gex))
    feature <- cbind(HNSC_methyl[commonCell,], HNSC_gex[commonCell,])
  }
  if (feat=="v17_HNSC_methyl.cnv"){
    load("data/pipeline/v17/bem/HNSC_methyl.ro")
    load("data/pipeline/v17/bem/HNSC_cnv.ro")
    commonCell <- intersect(rownames(HNSC_methyl), rownames(HNSC_cnv))
    feature <- cbind(HNSC_methyl[commonCell,], HNSC_cnv[commonCell,])
  }
  if (feat=="v17_HNSC_methyl.speed"){
    load("data/pipeline/v17/bem/HNSC_methyl.ro")
    load("data/pipeline/v17/bem/HNSC_speed.ro")
    commonCell <- intersect(rownames(HNSC_methyl), rownames(HNSC_speed))
    feature <- cbind(HNSC_methyl[commonCell,], HNSC_speed[commonCell,])
  }
  if (feat=="v17_HNSC_gex.cnv"){
    load("data/pipeline/v17/bem/HNSC_gex.ro")
    load("data/pipeline/v17/bem/HNSC_cnv.ro")
    commonCell <- intersect(rownames(HNSC_gex), rownames(HNSC_cnv))
    feature <- cbind(HNSC_gex[commonCell,], HNSC_cnv[commonCell,])
  }
  if (feat=="v17_HNSC_gex.speed"){
    load("data/pipeline/v17/bem/HNSC_gex.ro")
    load("data/pipeline/v17/bem/HNSC_speed.ro")
    commonCell <- intersect(rownames(HNSC_gex), rownames(HNSC_speed))
    feature <- cbind(HNSC_gex[commonCell,], HNSC_speed[commonCell,])
  }
  if (feat=="v17_HNSC_cnv.speed"){
    load("data/pipeline/v17/bem/HNSC_cnv.ro")
    load("data/pipeline/v17/bem/HNSC_speed.ro")
    commonCell <- intersect(rownames(HNSC_cnv), rownames(HNSC_speed))
    feature <- cbind(HNSC_cnv[commonCell,], HNSC_speed[commonCell,])
  }
  if (feat=="v17_HNSC_mut.methyl.gex"){
    load("data/pipeline/v17/bem/HNSC_mut.ro")
    load("data/pipeline/v17/bem/HNSC_methyl.ro")
    load("data/pipeline/v17/bem/HNSC_gex.ro")
    commonCell <- intersect(rownames(HNSC_mut), intersect(rownames(HNSC_methyl), rownames(HNSC_gex)))
    feature <- cbind(HNSC_mut[commonCell,], HNSC_methyl[commonCell,], HNSC_gex[commonCell,])
  }
  if (feat=="v17_HNSC_mut.methyl.cnv"){
    load("data/pipeline/v17/bem/HNSC_mut.ro")
    load("data/pipeline/v17/bem/HNSC_methyl.ro")
    load("data/pipeline/v17/bem/HNSC_cnv.ro")
    commonCell <- intersect(rownames(HNSC_mut), intersect(rownames(HNSC_methyl), rownames(HNSC_cnv)))
    feature <- cbind(HNSC_mut[commonCell,], HNSC_methyl[commonCell,], HNSC_cnv[commonCell,])
  }
  if (feat=="v17_HNSC_mut.methyl.speed"){
    load("data/pipeline/v17/bem/HNSC_mut.ro")
    load("data/pipeline/v17/bem/HNSC_methyl.ro")
    load("data/pipeline/v17/bem/HNSC_speed.ro")
    commonCell <- intersect(rownames(HNSC_mut), intersect(rownames(HNSC_methyl), rownames(HNSC_speed)))
    feature <- cbind(HNSC_mut[commonCell,], HNSC_methyl[commonCell,], HNSC_speed[commonCell,])
  }
  if (feat=="v17_HNSC_mut.gex.cnv"){
    load("data/pipeline/v17/bem/HNSC_mut.ro")
    load("data/pipeline/v17/bem/HNSC_gex.ro")
    load("data/pipeline/v17/bem/HNSC_cnv.ro")
    commonCell <- intersect(rownames(HNSC_mut), intersect(rownames(HNSC_gex), rownames(HNSC_cnv)))
    feature <- cbind(HNSC_mut[commonCell,], HNSC_gex[commonCell,], HNSC_cnv[commonCell,])
  }
  if (feat=="v17_HNSC_mut.gex.speed"){
    load("data/pipeline/v17/bem/HNSC_mut.ro")
    load("data/pipeline/v17/bem/HNSC_gex.ro")
    load("data/pipeline/v17/bem/HNSC_speed.ro")
    commonCell <- intersect(rownames(HNSC_mut), intersect(rownames(HNSC_gex), rownames(HNSC_speed)))
    feature <- cbind(HNSC_mut[commonCell,], HNSC_gex[commonCell,], HNSC_speed[commonCell,])
  }
  if (feat=="v17_HNSC_mut.cnv.speed"){
    load("data/pipeline/v17/bem/HNSC_mut.ro")
    load("data/pipeline/v17/bem/HNSC_cnv.ro")
    load("data/pipeline/v17/bem/HNSC_speed.ro")
    commonCell <- intersect(rownames(HNSC_mut), intersect(rownames(HNSC_cnv), rownames(HNSC_speed)))
    feature <- cbind(HNSC_mut[commonCell,], HNSC_cnv[commonCell,], HNSC_speed[commonCell,])
  }
  if (feat=="v17_HNSC_methyl.gex.cnv"){
    load("data/pipeline/v17/bem/HNSC_methyl.ro")
    load("data/pipeline/v17/bem/HNSC_gex.ro")
    load("data/pipeline/v17/bem/HNSC_cnv.ro")
    commonCell <- intersect(rownames(HNSC_methyl), intersect(rownames(HNSC_gex), rownames(HNSC_cnv)))
    feature <- cbind(HNSC_methyl[commonCell,], HNSC_gex[commonCell,], HNSC_cnv[commonCell,])
  }
  if (feat=="v17_HNSC_methyl.gex.speed"){
    load("data/pipeline/v17/bem/HNSC_methyl.ro")
    load("data/pipeline/v17/bem/HNSC_gex.ro")
    load("data/pipeline/v17/bem/HNSC_speed.ro")
    commonCell <- intersect(rownames(HNSC_methyl), intersect(rownames(HNSC_gex), rownames(HNSC_speed)))
    feature <- cbind(HNSC_methyl[commonCell,], HNSC_gex[commonCell,], HNSC_speed[commonCell,])
  }
  if (feat=="v17_HNSC_methyl.cnv.speed"){
    load("data/pipeline/v17/bem/HNSC_methyl.ro")
    load("data/pipeline/v17/bem/HNSC_cnv.ro")
    load("data/pipeline/v17/bem/HNSC_speed.ro")
    commonCell <- intersect(rownames(HNSC_methyl), intersect(rownames(HNSC_cnv), rownames(HNSC_speed)))
    feature <- cbind(HNSC_methyl[commonCell,], HNSC_cnv[commonCell,], HNSC_speed[commonCell,])
  }
  if (feat=="v17_HNSC_gex.cnv.speed"){
    load("data/pipeline/v17/bem/HNSC_gex.ro")
    load("data/pipeline/v17/bem/HNSC_cnv.ro")
    load("data/pipeline/v17/bem/HNSC_speed.ro")
    commonCell <- intersect(rownames(HNSC_gex), intersect(rownames(HNSC_cnv), rownames(HNSC_speed)))
    feature <- cbind(HNSC_gex[commonCell,], HNSC_cnv[commonCell,], HNSC_speed[commonCell,])
  }
  if (feat=="v17_HNSC_mut.methyl.gex.cnv"){
    load("data/pipeline/v17/bem/HNSC_mut.ro")
    load("data/pipeline/v17/bem/HNSC_methyl.ro")
    load("data/pipeline/v17/bem/HNSC_gex.ro")
    load("data/pipeline/v17/bem/HNSC_cnv.ro")
    commonCell <- intersect(rownames(HNSC_mut), intersect(rownames(HNSC_methyl), intersect(rownames(HNSC_gex), rownames(HNSC_cnv))))
    feature <- cbind(HNSC_mut[commonCell,], HNSC_methyl[commonCell,], HNSC_gex[commonCell,], HNSC_cnv[commonCell,])
  }
  if (feat=="v17_HNSC_mut.methyl.gex.speed"){
    load("data/pipeline/v17/bem/HNSC_mut.ro")
    load("data/pipeline/v17/bem/HNSC_methyl.ro")
    load("data/pipeline/v17/bem/HNSC_gex.ro")
    load("data/pipeline/v17/bem/HNSC_speed.ro")
    commonCell <- intersect(rownames(HNSC_mut), intersect(rownames(HNSC_methyl), intersect(rownames(HNSC_gex), rownames(HNSC_speed))))
    feature <- cbind(HNSC_mut[commonCell,], HNSC_methyl[commonCell,], HNSC_gex[commonCell,], HNSC_speed[commonCell,])
  }
  if (feat=="v17_HNSC_mut.methyl.cnv.speed"){
    load("data/pipeline/v17/bem/HNSC_mut.ro")
    load("data/pipeline/v17/bem/HNSC_methyl.ro")
    load("data/pipeline/v17/bem/HNSC_cnv.ro")
    load("data/pipeline/v17/bem/HNSC_speed.ro")
    commonCell <- intersect(rownames(HNSC_mut), intersect(rownames(HNSC_methyl), intersect(rownames(HNSC_cnv), rownames(HNSC_speed))))
    feature <- cbind(HNSC_mut[commonCell,], HNSC_methyl[commonCell,], HNSC_cnv[commonCell,], HNSC_speed[commonCell,])
  }
  if (feat=="v17_HNSC_mut.gex.cnv.speed"){
    load("data/pipeline/v17/bem/HNSC_mut.ro")
    load("data/pipeline/v17/bem/HNSC_gex.ro")
    load("data/pipeline/v17/bem/HNSC_cnv.ro")
    load("data/pipeline/v17/bem/HNSC_speed.ro")
    commonCell <- intersect(rownames(HNSC_mut), intersect(rownames(HNSC_gex), intersect(rownames(HNSC_cnv), rownames(HNSC_speed))))
    feature <- cbind(HNSC_mut[commonCell,], HNSC_gex[commonCell,], HNSC_cnv[commonCell,], HNSC_speed[commonCell,])
  }
  if (feat=="v17_HNSC_methyl.gex.cnv.speed"){
    load("data/pipeline/v17/bem/HNSC_methyl.ro")
    load("data/pipeline/v17/bem/HNSC_gex.ro")
    load("data/pipeline/v17/bem/HNSC_cnv.ro")
    load("data/pipeline/v17/bem/HNSC_speed.ro")
    commonCell <- intersect(rownames(HNSC_methyl), intersect(rownames(HNSC_gex), intersect(rownames(HNSC_cnv), rownames(HNSC_speed))))
    feature <- cbind(HNSC_methyl[commonCell,], HNSC_gex[commonCell,], HNSC_cnv[commonCell,], HNSC_speed[commonCell,])
  }
  if (feat=="v17_HNSC_mut.methyl.gex.cnv.speed"){
    load("data/pipeline/v17/bem/HNSC_mut.ro")
    load("data/pipeline/v17/bem/HNSC_methyl.ro")
    load("data/pipeline/v17/bem/HNSC_gex.ro")
    load("data/pipeline/v17/bem/HNSC_cnv.ro")
    load("data/pipeline/v17/bem/HNSC_speed.ro")
    commonCell <- intersect(rownames(HNSC_mut), intersect(rownames(HNSC_methyl), intersect(rownames(HNSC_gex), intersect(rownames(HNSC_cnv), rownames(HNSC_speed)))))
    feature <- cbind(HNSC_mut[commonCell,], HNSC_methyl[commonCell,], HNSC_gex[commonCell,], HNSC_cnv[commonCell,], HNSC_speed[commonCell,])
  }
  if (feat=="v17_KIRC_mut"){
    load("data/pipeline/v17/bem/KIRC_mut.ro")
    commonCell <- rownames(KIRC_mut)
    feature <- cbind(KIRC_mut[commonCell,])
  }
  if (feat=="v17_KIRC_methyl"){
    load("data/pipeline/v17/bem/KIRC_methyl.ro")
    commonCell <- rownames(KIRC_methyl)
    feature <- cbind(KIRC_methyl[commonCell,])
  }
  if (feat=="v17_KIRC_gex"){
    load("data/pipeline/v17/bem/KIRC_gex.ro")
    commonCell <- rownames(KIRC_gex)
    feature <- cbind(KIRC_gex[commonCell,])
  }
  if (feat=="v17_KIRC_cnv"){
    load("data/pipeline/v17/bem/KIRC_cnv.ro")
    commonCell <- rownames(KIRC_cnv)
    feature <- cbind(KIRC_cnv[commonCell,])
  }
  if (feat=="v17_KIRC_speed"){
    load("data/pipeline/v17/bem/KIRC_speed.ro")
    commonCell <- rownames(KIRC_speed)
    feature <- cbind(KIRC_speed[commonCell,])
  }
  if (feat=="v17_KIRC_mut.methyl"){
    load("data/pipeline/v17/bem/KIRC_mut.ro")
    load("data/pipeline/v17/bem/KIRC_methyl.ro")
    commonCell <- intersect(rownames(KIRC_mut), rownames(KIRC_methyl))
    feature <- cbind(KIRC_mut[commonCell,], KIRC_methyl[commonCell,])
  }
  if (feat=="v17_KIRC_mut.gex"){
    load("data/pipeline/v17/bem/KIRC_mut.ro")
    load("data/pipeline/v17/bem/KIRC_gex.ro")
    commonCell <- intersect(rownames(KIRC_mut), rownames(KIRC_gex))
    feature <- cbind(KIRC_mut[commonCell,], KIRC_gex[commonCell,])
  }
  if (feat=="v17_KIRC_mut.cnv"){
    load("data/pipeline/v17/bem/KIRC_mut.ro")
    load("data/pipeline/v17/bem/KIRC_cnv.ro")
    commonCell <- intersect(rownames(KIRC_mut), rownames(KIRC_cnv))
    feature <- cbind(KIRC_mut[commonCell,], KIRC_cnv[commonCell,])
  }
  if (feat=="v17_KIRC_mut.speed"){
    load("data/pipeline/v17/bem/KIRC_mut.ro")
    load("data/pipeline/v17/bem/KIRC_speed.ro")
    commonCell <- intersect(rownames(KIRC_mut), rownames(KIRC_speed))
    feature <- cbind(KIRC_mut[commonCell,], KIRC_speed[commonCell,])
  }
  if (feat=="v17_KIRC_methyl.gex"){
    load("data/pipeline/v17/bem/KIRC_methyl.ro")
    load("data/pipeline/v17/bem/KIRC_gex.ro")
    commonCell <- intersect(rownames(KIRC_methyl), rownames(KIRC_gex))
    feature <- cbind(KIRC_methyl[commonCell,], KIRC_gex[commonCell,])
  }
  if (feat=="v17_KIRC_methyl.cnv"){
    load("data/pipeline/v17/bem/KIRC_methyl.ro")
    load("data/pipeline/v17/bem/KIRC_cnv.ro")
    commonCell <- intersect(rownames(KIRC_methyl), rownames(KIRC_cnv))
    feature <- cbind(KIRC_methyl[commonCell,], KIRC_cnv[commonCell,])
  }
  if (feat=="v17_KIRC_methyl.speed"){
    load("data/pipeline/v17/bem/KIRC_methyl.ro")
    load("data/pipeline/v17/bem/KIRC_speed.ro")
    commonCell <- intersect(rownames(KIRC_methyl), rownames(KIRC_speed))
    feature <- cbind(KIRC_methyl[commonCell,], KIRC_speed[commonCell,])
  }
  if (feat=="v17_KIRC_gex.cnv"){
    load("data/pipeline/v17/bem/KIRC_gex.ro")
    load("data/pipeline/v17/bem/KIRC_cnv.ro")
    commonCell <- intersect(rownames(KIRC_gex), rownames(KIRC_cnv))
    feature <- cbind(KIRC_gex[commonCell,], KIRC_cnv[commonCell,])
  }
  if (feat=="v17_KIRC_gex.speed"){
    load("data/pipeline/v17/bem/KIRC_gex.ro")
    load("data/pipeline/v17/bem/KIRC_speed.ro")
    commonCell <- intersect(rownames(KIRC_gex), rownames(KIRC_speed))
    feature <- cbind(KIRC_gex[commonCell,], KIRC_speed[commonCell,])
  }
  if (feat=="v17_KIRC_cnv.speed"){
    load("data/pipeline/v17/bem/KIRC_cnv.ro")
    load("data/pipeline/v17/bem/KIRC_speed.ro")
    commonCell <- intersect(rownames(KIRC_cnv), rownames(KIRC_speed))
    feature <- cbind(KIRC_cnv[commonCell,], KIRC_speed[commonCell,])
  }
  if (feat=="v17_KIRC_mut.methyl.gex"){
    load("data/pipeline/v17/bem/KIRC_mut.ro")
    load("data/pipeline/v17/bem/KIRC_methyl.ro")
    load("data/pipeline/v17/bem/KIRC_gex.ro")
    commonCell <- intersect(rownames(KIRC_mut), intersect(rownames(KIRC_methyl), rownames(KIRC_gex)))
    feature <- cbind(KIRC_mut[commonCell,], KIRC_methyl[commonCell,], KIRC_gex[commonCell,])
  }
  if (feat=="v17_KIRC_mut.methyl.cnv"){
    load("data/pipeline/v17/bem/KIRC_mut.ro")
    load("data/pipeline/v17/bem/KIRC_methyl.ro")
    load("data/pipeline/v17/bem/KIRC_cnv.ro")
    commonCell <- intersect(rownames(KIRC_mut), intersect(rownames(KIRC_methyl), rownames(KIRC_cnv)))
    feature <- cbind(KIRC_mut[commonCell,], KIRC_methyl[commonCell,], KIRC_cnv[commonCell,])
  }
  if (feat=="v17_KIRC_mut.methyl.speed"){
    load("data/pipeline/v17/bem/KIRC_mut.ro")
    load("data/pipeline/v17/bem/KIRC_methyl.ro")
    load("data/pipeline/v17/bem/KIRC_speed.ro")
    commonCell <- intersect(rownames(KIRC_mut), intersect(rownames(KIRC_methyl), rownames(KIRC_speed)))
    feature <- cbind(KIRC_mut[commonCell,], KIRC_methyl[commonCell,], KIRC_speed[commonCell,])
  }
  if (feat=="v17_KIRC_mut.gex.cnv"){
    load("data/pipeline/v17/bem/KIRC_mut.ro")
    load("data/pipeline/v17/bem/KIRC_gex.ro")
    load("data/pipeline/v17/bem/KIRC_cnv.ro")
    commonCell <- intersect(rownames(KIRC_mut), intersect(rownames(KIRC_gex), rownames(KIRC_cnv)))
    feature <- cbind(KIRC_mut[commonCell,], KIRC_gex[commonCell,], KIRC_cnv[commonCell,])
  }
  if (feat=="v17_KIRC_mut.gex.speed"){
    load("data/pipeline/v17/bem/KIRC_mut.ro")
    load("data/pipeline/v17/bem/KIRC_gex.ro")
    load("data/pipeline/v17/bem/KIRC_speed.ro")
    commonCell <- intersect(rownames(KIRC_mut), intersect(rownames(KIRC_gex), rownames(KIRC_speed)))
    feature <- cbind(KIRC_mut[commonCell,], KIRC_gex[commonCell,], KIRC_speed[commonCell,])
  }
  if (feat=="v17_KIRC_mut.cnv.speed"){
    load("data/pipeline/v17/bem/KIRC_mut.ro")
    load("data/pipeline/v17/bem/KIRC_cnv.ro")
    load("data/pipeline/v17/bem/KIRC_speed.ro")
    commonCell <- intersect(rownames(KIRC_mut), intersect(rownames(KIRC_cnv), rownames(KIRC_speed)))
    feature <- cbind(KIRC_mut[commonCell,], KIRC_cnv[commonCell,], KIRC_speed[commonCell,])
  }
  if (feat=="v17_KIRC_methyl.gex.cnv"){
    load("data/pipeline/v17/bem/KIRC_methyl.ro")
    load("data/pipeline/v17/bem/KIRC_gex.ro")
    load("data/pipeline/v17/bem/KIRC_cnv.ro")
    commonCell <- intersect(rownames(KIRC_methyl), intersect(rownames(KIRC_gex), rownames(KIRC_cnv)))
    feature <- cbind(KIRC_methyl[commonCell,], KIRC_gex[commonCell,], KIRC_cnv[commonCell,])
  }
  if (feat=="v17_KIRC_methyl.gex.speed"){
    load("data/pipeline/v17/bem/KIRC_methyl.ro")
    load("data/pipeline/v17/bem/KIRC_gex.ro")
    load("data/pipeline/v17/bem/KIRC_speed.ro")
    commonCell <- intersect(rownames(KIRC_methyl), intersect(rownames(KIRC_gex), rownames(KIRC_speed)))
    feature <- cbind(KIRC_methyl[commonCell,], KIRC_gex[commonCell,], KIRC_speed[commonCell,])
  }
  if (feat=="v17_KIRC_methyl.cnv.speed"){
    load("data/pipeline/v17/bem/KIRC_methyl.ro")
    load("data/pipeline/v17/bem/KIRC_cnv.ro")
    load("data/pipeline/v17/bem/KIRC_speed.ro")
    commonCell <- intersect(rownames(KIRC_methyl), intersect(rownames(KIRC_cnv), rownames(KIRC_speed)))
    feature <- cbind(KIRC_methyl[commonCell,], KIRC_cnv[commonCell,], KIRC_speed[commonCell,])
  }
  if (feat=="v17_KIRC_gex.cnv.speed"){
    load("data/pipeline/v17/bem/KIRC_gex.ro")
    load("data/pipeline/v17/bem/KIRC_cnv.ro")
    load("data/pipeline/v17/bem/KIRC_speed.ro")
    commonCell <- intersect(rownames(KIRC_gex), intersect(rownames(KIRC_cnv), rownames(KIRC_speed)))
    feature <- cbind(KIRC_gex[commonCell,], KIRC_cnv[commonCell,], KIRC_speed[commonCell,])
  }
  if (feat=="v17_KIRC_mut.methyl.gex.cnv"){
    load("data/pipeline/v17/bem/KIRC_mut.ro")
    load("data/pipeline/v17/bem/KIRC_methyl.ro")
    load("data/pipeline/v17/bem/KIRC_gex.ro")
    load("data/pipeline/v17/bem/KIRC_cnv.ro")
    commonCell <- intersect(rownames(KIRC_mut), intersect(rownames(KIRC_methyl), intersect(rownames(KIRC_gex), rownames(KIRC_cnv))))
    feature <- cbind(KIRC_mut[commonCell,], KIRC_methyl[commonCell,], KIRC_gex[commonCell,], KIRC_cnv[commonCell,])
  }
  if (feat=="v17_KIRC_mut.methyl.gex.speed"){
    load("data/pipeline/v17/bem/KIRC_mut.ro")
    load("data/pipeline/v17/bem/KIRC_methyl.ro")
    load("data/pipeline/v17/bem/KIRC_gex.ro")
    load("data/pipeline/v17/bem/KIRC_speed.ro")
    commonCell <- intersect(rownames(KIRC_mut), intersect(rownames(KIRC_methyl), intersect(rownames(KIRC_gex), rownames(KIRC_speed))))
    feature <- cbind(KIRC_mut[commonCell,], KIRC_methyl[commonCell,], KIRC_gex[commonCell,], KIRC_speed[commonCell,])
  }
  if (feat=="v17_KIRC_mut.methyl.cnv.speed"){
    load("data/pipeline/v17/bem/KIRC_mut.ro")
    load("data/pipeline/v17/bem/KIRC_methyl.ro")
    load("data/pipeline/v17/bem/KIRC_cnv.ro")
    load("data/pipeline/v17/bem/KIRC_speed.ro")
    commonCell <- intersect(rownames(KIRC_mut), intersect(rownames(KIRC_methyl), intersect(rownames(KIRC_cnv), rownames(KIRC_speed))))
    feature <- cbind(KIRC_mut[commonCell,], KIRC_methyl[commonCell,], KIRC_cnv[commonCell,], KIRC_speed[commonCell,])
  }
  if (feat=="v17_KIRC_mut.gex.cnv.speed"){
    load("data/pipeline/v17/bem/KIRC_mut.ro")
    load("data/pipeline/v17/bem/KIRC_gex.ro")
    load("data/pipeline/v17/bem/KIRC_cnv.ro")
    load("data/pipeline/v17/bem/KIRC_speed.ro")
    commonCell <- intersect(rownames(KIRC_mut), intersect(rownames(KIRC_gex), intersect(rownames(KIRC_cnv), rownames(KIRC_speed))))
    feature <- cbind(KIRC_mut[commonCell,], KIRC_gex[commonCell,], KIRC_cnv[commonCell,], KIRC_speed[commonCell,])
  }
  if (feat=="v17_KIRC_methyl.gex.cnv.speed"){
    load("data/pipeline/v17/bem/KIRC_methyl.ro")
    load("data/pipeline/v17/bem/KIRC_gex.ro")
    load("data/pipeline/v17/bem/KIRC_cnv.ro")
    load("data/pipeline/v17/bem/KIRC_speed.ro")
    commonCell <- intersect(rownames(KIRC_methyl), intersect(rownames(KIRC_gex), intersect(rownames(KIRC_cnv), rownames(KIRC_speed))))
    feature <- cbind(KIRC_methyl[commonCell,], KIRC_gex[commonCell,], KIRC_cnv[commonCell,], KIRC_speed[commonCell,])
  }
  if (feat=="v17_KIRC_mut.methyl.gex.cnv.speed"){
    load("data/pipeline/v17/bem/KIRC_mut.ro")
    load("data/pipeline/v17/bem/KIRC_methyl.ro")
    load("data/pipeline/v17/bem/KIRC_gex.ro")
    load("data/pipeline/v17/bem/KIRC_cnv.ro")
    load("data/pipeline/v17/bem/KIRC_speed.ro")
    commonCell <- intersect(rownames(KIRC_mut), intersect(rownames(KIRC_methyl), intersect(rownames(KIRC_gex), intersect(rownames(KIRC_cnv), rownames(KIRC_speed)))))
    feature <- cbind(KIRC_mut[commonCell,], KIRC_methyl[commonCell,], KIRC_gex[commonCell,], KIRC_cnv[commonCell,], KIRC_speed[commonCell,])
  }
  if (feat=="v17_LAML_mut"){
    load("data/pipeline/v17/bem/LAML_mut.ro")
    commonCell <- rownames(LAML_mut)
    feature <- cbind(LAML_mut[commonCell,])
  }
  if (feat=="v17_LAML_methyl"){
    load("data/pipeline/v17/bem/LAML_methyl.ro")
    commonCell <- rownames(LAML_methyl)
    feature <- cbind(LAML_methyl[commonCell,])
  }
  if (feat=="v17_LAML_gex"){
    load("data/pipeline/v17/bem/LAML_gex.ro")
    commonCell <- rownames(LAML_gex)
    feature <- cbind(LAML_gex[commonCell,])
  }
  if (feat=="v17_LAML_cnv"){
    load("data/pipeline/v17/bem/LAML_cnv.ro")
    commonCell <- rownames(LAML_cnv)
    feature <- cbind(LAML_cnv[commonCell,])
  }
  if (feat=="v17_LAML_speed"){
    load("data/pipeline/v17/bem/LAML_speed.ro")
    commonCell <- rownames(LAML_speed)
    feature <- cbind(LAML_speed[commonCell,])
  }
  if (feat=="v17_LAML_mut.methyl"){
    load("data/pipeline/v17/bem/LAML_mut.ro")
    load("data/pipeline/v17/bem/LAML_methyl.ro")
    commonCell <- intersect(rownames(LAML_mut), rownames(LAML_methyl))
    feature <- cbind(LAML_mut[commonCell,], LAML_methyl[commonCell,])
  }
  if (feat=="v17_LAML_mut.gex"){
    load("data/pipeline/v17/bem/LAML_mut.ro")
    load("data/pipeline/v17/bem/LAML_gex.ro")
    commonCell <- intersect(rownames(LAML_mut), rownames(LAML_gex))
    feature <- cbind(LAML_mut[commonCell,], LAML_gex[commonCell,])
  }
  if (feat=="v17_LAML_mut.cnv"){
    load("data/pipeline/v17/bem/LAML_mut.ro")
    load("data/pipeline/v17/bem/LAML_cnv.ro")
    commonCell <- intersect(rownames(LAML_mut), rownames(LAML_cnv))
    feature <- cbind(LAML_mut[commonCell,], LAML_cnv[commonCell,])
  }
  if (feat=="v17_LAML_mut.speed"){
    load("data/pipeline/v17/bem/LAML_mut.ro")
    load("data/pipeline/v17/bem/LAML_speed.ro")
    commonCell <- intersect(rownames(LAML_mut), rownames(LAML_speed))
    feature <- cbind(LAML_mut[commonCell,], LAML_speed[commonCell,])
  }
  if (feat=="v17_LAML_methyl.gex"){
    load("data/pipeline/v17/bem/LAML_methyl.ro")
    load("data/pipeline/v17/bem/LAML_gex.ro")
    commonCell <- intersect(rownames(LAML_methyl), rownames(LAML_gex))
    feature <- cbind(LAML_methyl[commonCell,], LAML_gex[commonCell,])
  }
  if (feat=="v17_LAML_methyl.cnv"){
    load("data/pipeline/v17/bem/LAML_methyl.ro")
    load("data/pipeline/v17/bem/LAML_cnv.ro")
    commonCell <- intersect(rownames(LAML_methyl), rownames(LAML_cnv))
    feature <- cbind(LAML_methyl[commonCell,], LAML_cnv[commonCell,])
  }
  if (feat=="v17_LAML_methyl.speed"){
    load("data/pipeline/v17/bem/LAML_methyl.ro")
    load("data/pipeline/v17/bem/LAML_speed.ro")
    commonCell <- intersect(rownames(LAML_methyl), rownames(LAML_speed))
    feature <- cbind(LAML_methyl[commonCell,], LAML_speed[commonCell,])
  }
  if (feat=="v17_LAML_gex.cnv"){
    load("data/pipeline/v17/bem/LAML_gex.ro")
    load("data/pipeline/v17/bem/LAML_cnv.ro")
    commonCell <- intersect(rownames(LAML_gex), rownames(LAML_cnv))
    feature <- cbind(LAML_gex[commonCell,], LAML_cnv[commonCell,])
  }
  if (feat=="v17_LAML_gex.speed"){
    load("data/pipeline/v17/bem/LAML_gex.ro")
    load("data/pipeline/v17/bem/LAML_speed.ro")
    commonCell <- intersect(rownames(LAML_gex), rownames(LAML_speed))
    feature <- cbind(LAML_gex[commonCell,], LAML_speed[commonCell,])
  }
  if (feat=="v17_LAML_cnv.speed"){
    load("data/pipeline/v17/bem/LAML_cnv.ro")
    load("data/pipeline/v17/bem/LAML_speed.ro")
    commonCell <- intersect(rownames(LAML_cnv), rownames(LAML_speed))
    feature <- cbind(LAML_cnv[commonCell,], LAML_speed[commonCell,])
  }
  if (feat=="v17_LAML_mut.methyl.gex"){
    load("data/pipeline/v17/bem/LAML_mut.ro")
    load("data/pipeline/v17/bem/LAML_methyl.ro")
    load("data/pipeline/v17/bem/LAML_gex.ro")
    commonCell <- intersect(rownames(LAML_mut), intersect(rownames(LAML_methyl), rownames(LAML_gex)))
    feature <- cbind(LAML_mut[commonCell,], LAML_methyl[commonCell,], LAML_gex[commonCell,])
  }
  if (feat=="v17_LAML_mut.methyl.cnv"){
    load("data/pipeline/v17/bem/LAML_mut.ro")
    load("data/pipeline/v17/bem/LAML_methyl.ro")
    load("data/pipeline/v17/bem/LAML_cnv.ro")
    commonCell <- intersect(rownames(LAML_mut), intersect(rownames(LAML_methyl), rownames(LAML_cnv)))
    feature <- cbind(LAML_mut[commonCell,], LAML_methyl[commonCell,], LAML_cnv[commonCell,])
  }
  if (feat=="v17_LAML_mut.methyl.speed"){
    load("data/pipeline/v17/bem/LAML_mut.ro")
    load("data/pipeline/v17/bem/LAML_methyl.ro")
    load("data/pipeline/v17/bem/LAML_speed.ro")
    commonCell <- intersect(rownames(LAML_mut), intersect(rownames(LAML_methyl), rownames(LAML_speed)))
    feature <- cbind(LAML_mut[commonCell,], LAML_methyl[commonCell,], LAML_speed[commonCell,])
  }
  if (feat=="v17_LAML_mut.gex.cnv"){
    load("data/pipeline/v17/bem/LAML_mut.ro")
    load("data/pipeline/v17/bem/LAML_gex.ro")
    load("data/pipeline/v17/bem/LAML_cnv.ro")
    commonCell <- intersect(rownames(LAML_mut), intersect(rownames(LAML_gex), rownames(LAML_cnv)))
    feature <- cbind(LAML_mut[commonCell,], LAML_gex[commonCell,], LAML_cnv[commonCell,])
  }
  if (feat=="v17_LAML_mut.gex.speed"){
    load("data/pipeline/v17/bem/LAML_mut.ro")
    load("data/pipeline/v17/bem/LAML_gex.ro")
    load("data/pipeline/v17/bem/LAML_speed.ro")
    commonCell <- intersect(rownames(LAML_mut), intersect(rownames(LAML_gex), rownames(LAML_speed)))
    feature <- cbind(LAML_mut[commonCell,], LAML_gex[commonCell,], LAML_speed[commonCell,])
  }
  if (feat=="v17_LAML_mut.cnv.speed"){
    load("data/pipeline/v17/bem/LAML_mut.ro")
    load("data/pipeline/v17/bem/LAML_cnv.ro")
    load("data/pipeline/v17/bem/LAML_speed.ro")
    commonCell <- intersect(rownames(LAML_mut), intersect(rownames(LAML_cnv), rownames(LAML_speed)))
    feature <- cbind(LAML_mut[commonCell,], LAML_cnv[commonCell,], LAML_speed[commonCell,])
  }
  if (feat=="v17_LAML_methyl.gex.cnv"){
    load("data/pipeline/v17/bem/LAML_methyl.ro")
    load("data/pipeline/v17/bem/LAML_gex.ro")
    load("data/pipeline/v17/bem/LAML_cnv.ro")
    commonCell <- intersect(rownames(LAML_methyl), intersect(rownames(LAML_gex), rownames(LAML_cnv)))
    feature <- cbind(LAML_methyl[commonCell,], LAML_gex[commonCell,], LAML_cnv[commonCell,])
  }
  if (feat=="v17_LAML_methyl.gex.speed"){
    load("data/pipeline/v17/bem/LAML_methyl.ro")
    load("data/pipeline/v17/bem/LAML_gex.ro")
    load("data/pipeline/v17/bem/LAML_speed.ro")
    commonCell <- intersect(rownames(LAML_methyl), intersect(rownames(LAML_gex), rownames(LAML_speed)))
    feature <- cbind(LAML_methyl[commonCell,], LAML_gex[commonCell,], LAML_speed[commonCell,])
  }
  if (feat=="v17_LAML_methyl.cnv.speed"){
    load("data/pipeline/v17/bem/LAML_methyl.ro")
    load("data/pipeline/v17/bem/LAML_cnv.ro")
    load("data/pipeline/v17/bem/LAML_speed.ro")
    commonCell <- intersect(rownames(LAML_methyl), intersect(rownames(LAML_cnv), rownames(LAML_speed)))
    feature <- cbind(LAML_methyl[commonCell,], LAML_cnv[commonCell,], LAML_speed[commonCell,])
  }
  if (feat=="v17_LAML_gex.cnv.speed"){
    load("data/pipeline/v17/bem/LAML_gex.ro")
    load("data/pipeline/v17/bem/LAML_cnv.ro")
    load("data/pipeline/v17/bem/LAML_speed.ro")
    commonCell <- intersect(rownames(LAML_gex), intersect(rownames(LAML_cnv), rownames(LAML_speed)))
    feature <- cbind(LAML_gex[commonCell,], LAML_cnv[commonCell,], LAML_speed[commonCell,])
  }
  if (feat=="v17_LAML_mut.methyl.gex.cnv"){
    load("data/pipeline/v17/bem/LAML_mut.ro")
    load("data/pipeline/v17/bem/LAML_methyl.ro")
    load("data/pipeline/v17/bem/LAML_gex.ro")
    load("data/pipeline/v17/bem/LAML_cnv.ro")
    commonCell <- intersect(rownames(LAML_mut), intersect(rownames(LAML_methyl), intersect(rownames(LAML_gex), rownames(LAML_cnv))))
    feature <- cbind(LAML_mut[commonCell,], LAML_methyl[commonCell,], LAML_gex[commonCell,], LAML_cnv[commonCell,])
  }
  if (feat=="v17_LAML_mut.methyl.gex.speed"){
    load("data/pipeline/v17/bem/LAML_mut.ro")
    load("data/pipeline/v17/bem/LAML_methyl.ro")
    load("data/pipeline/v17/bem/LAML_gex.ro")
    load("data/pipeline/v17/bem/LAML_speed.ro")
    commonCell <- intersect(rownames(LAML_mut), intersect(rownames(LAML_methyl), intersect(rownames(LAML_gex), rownames(LAML_speed))))
    feature <- cbind(LAML_mut[commonCell,], LAML_methyl[commonCell,], LAML_gex[commonCell,], LAML_speed[commonCell,])
  }
  if (feat=="v17_LAML_mut.methyl.cnv.speed"){
    load("data/pipeline/v17/bem/LAML_mut.ro")
    load("data/pipeline/v17/bem/LAML_methyl.ro")
    load("data/pipeline/v17/bem/LAML_cnv.ro")
    load("data/pipeline/v17/bem/LAML_speed.ro")
    commonCell <- intersect(rownames(LAML_mut), intersect(rownames(LAML_methyl), intersect(rownames(LAML_cnv), rownames(LAML_speed))))
    feature <- cbind(LAML_mut[commonCell,], LAML_methyl[commonCell,], LAML_cnv[commonCell,], LAML_speed[commonCell,])
  }
  if (feat=="v17_LAML_mut.gex.cnv.speed"){
    load("data/pipeline/v17/bem/LAML_mut.ro")
    load("data/pipeline/v17/bem/LAML_gex.ro")
    load("data/pipeline/v17/bem/LAML_cnv.ro")
    load("data/pipeline/v17/bem/LAML_speed.ro")
    commonCell <- intersect(rownames(LAML_mut), intersect(rownames(LAML_gex), intersect(rownames(LAML_cnv), rownames(LAML_speed))))
    feature <- cbind(LAML_mut[commonCell,], LAML_gex[commonCell,], LAML_cnv[commonCell,], LAML_speed[commonCell,])
  }
  if (feat=="v17_LAML_methyl.gex.cnv.speed"){
    load("data/pipeline/v17/bem/LAML_methyl.ro")
    load("data/pipeline/v17/bem/LAML_gex.ro")
    load("data/pipeline/v17/bem/LAML_cnv.ro")
    load("data/pipeline/v17/bem/LAML_speed.ro")
    commonCell <- intersect(rownames(LAML_methyl), intersect(rownames(LAML_gex), intersect(rownames(LAML_cnv), rownames(LAML_speed))))
    feature <- cbind(LAML_methyl[commonCell,], LAML_gex[commonCell,], LAML_cnv[commonCell,], LAML_speed[commonCell,])
  }
  if (feat=="v17_LAML_mut.methyl.gex.cnv.speed"){
    load("data/pipeline/v17/bem/LAML_mut.ro")
    load("data/pipeline/v17/bem/LAML_methyl.ro")
    load("data/pipeline/v17/bem/LAML_gex.ro")
    load("data/pipeline/v17/bem/LAML_cnv.ro")
    load("data/pipeline/v17/bem/LAML_speed.ro")
    commonCell <- intersect(rownames(LAML_mut), intersect(rownames(LAML_methyl), intersect(rownames(LAML_gex), intersect(rownames(LAML_cnv), rownames(LAML_speed)))))
    feature <- cbind(LAML_mut[commonCell,], LAML_methyl[commonCell,], LAML_gex[commonCell,], LAML_cnv[commonCell,], LAML_speed[commonCell,])
  }
  if (feat=="v17_LGG_mut"){
    load("data/pipeline/v17/bem/LGG_mut.ro")
    commonCell <- rownames(LGG_mut)
    feature <- cbind(LGG_mut[commonCell,])
  }
  if (feat=="v17_LGG_methyl"){
    load("data/pipeline/v17/bem/LGG_methyl.ro")
    commonCell <- rownames(LGG_methyl)
    feature <- cbind(LGG_methyl[commonCell,])
  }
  if (feat=="v17_LGG_gex"){
    load("data/pipeline/v17/bem/LGG_gex.ro")
    commonCell <- rownames(LGG_gex)
    feature <- cbind(LGG_gex[commonCell,])
  }
  if (feat=="v17_LGG_cnv"){
    load("data/pipeline/v17/bem/LGG_cnv.ro")
    commonCell <- rownames(LGG_cnv)
    feature <- cbind(LGG_cnv[commonCell,])
  }
  if (feat=="v17_LGG_speed"){
    load("data/pipeline/v17/bem/LGG_speed.ro")
    commonCell <- rownames(LGG_speed)
    feature <- cbind(LGG_speed[commonCell,])
  }
  if (feat=="v17_LGG_mut.methyl"){
    load("data/pipeline/v17/bem/LGG_mut.ro")
    load("data/pipeline/v17/bem/LGG_methyl.ro")
    commonCell <- intersect(rownames(LGG_mut), rownames(LGG_methyl))
    feature <- cbind(LGG_mut[commonCell,], LGG_methyl[commonCell,])
  }
  if (feat=="v17_LGG_mut.gex"){
    load("data/pipeline/v17/bem/LGG_mut.ro")
    load("data/pipeline/v17/bem/LGG_gex.ro")
    commonCell <- intersect(rownames(LGG_mut), rownames(LGG_gex))
    feature <- cbind(LGG_mut[commonCell,], LGG_gex[commonCell,])
  }
  if (feat=="v17_LGG_mut.cnv"){
    load("data/pipeline/v17/bem/LGG_mut.ro")
    load("data/pipeline/v17/bem/LGG_cnv.ro")
    commonCell <- intersect(rownames(LGG_mut), rownames(LGG_cnv))
    feature <- cbind(LGG_mut[commonCell,], LGG_cnv[commonCell,])
  }
  if (feat=="v17_LGG_mut.speed"){
    load("data/pipeline/v17/bem/LGG_mut.ro")
    load("data/pipeline/v17/bem/LGG_speed.ro")
    commonCell <- intersect(rownames(LGG_mut), rownames(LGG_speed))
    feature <- cbind(LGG_mut[commonCell,], LGG_speed[commonCell,])
  }
  if (feat=="v17_LGG_methyl.gex"){
    load("data/pipeline/v17/bem/LGG_methyl.ro")
    load("data/pipeline/v17/bem/LGG_gex.ro")
    commonCell <- intersect(rownames(LGG_methyl), rownames(LGG_gex))
    feature <- cbind(LGG_methyl[commonCell,], LGG_gex[commonCell,])
  }
  if (feat=="v17_LGG_methyl.cnv"){
    load("data/pipeline/v17/bem/LGG_methyl.ro")
    load("data/pipeline/v17/bem/LGG_cnv.ro")
    commonCell <- intersect(rownames(LGG_methyl), rownames(LGG_cnv))
    feature <- cbind(LGG_methyl[commonCell,], LGG_cnv[commonCell,])
  }
  if (feat=="v17_LGG_methyl.speed"){
    load("data/pipeline/v17/bem/LGG_methyl.ro")
    load("data/pipeline/v17/bem/LGG_speed.ro")
    commonCell <- intersect(rownames(LGG_methyl), rownames(LGG_speed))
    feature <- cbind(LGG_methyl[commonCell,], LGG_speed[commonCell,])
  }
  if (feat=="v17_LGG_gex.cnv"){
    load("data/pipeline/v17/bem/LGG_gex.ro")
    load("data/pipeline/v17/bem/LGG_cnv.ro")
    commonCell <- intersect(rownames(LGG_gex), rownames(LGG_cnv))
    feature <- cbind(LGG_gex[commonCell,], LGG_cnv[commonCell,])
  }
  if (feat=="v17_LGG_gex.speed"){
    load("data/pipeline/v17/bem/LGG_gex.ro")
    load("data/pipeline/v17/bem/LGG_speed.ro")
    commonCell <- intersect(rownames(LGG_gex), rownames(LGG_speed))
    feature <- cbind(LGG_gex[commonCell,], LGG_speed[commonCell,])
  }
  if (feat=="v17_LGG_cnv.speed"){
    load("data/pipeline/v17/bem/LGG_cnv.ro")
    load("data/pipeline/v17/bem/LGG_speed.ro")
    commonCell <- intersect(rownames(LGG_cnv), rownames(LGG_speed))
    feature <- cbind(LGG_cnv[commonCell,], LGG_speed[commonCell,])
  }
  if (feat=="v17_LGG_mut.methyl.gex"){
    load("data/pipeline/v17/bem/LGG_mut.ro")
    load("data/pipeline/v17/bem/LGG_methyl.ro")
    load("data/pipeline/v17/bem/LGG_gex.ro")
    commonCell <- intersect(rownames(LGG_mut), intersect(rownames(LGG_methyl), rownames(LGG_gex)))
    feature <- cbind(LGG_mut[commonCell,], LGG_methyl[commonCell,], LGG_gex[commonCell,])
  }
  if (feat=="v17_LGG_mut.methyl.cnv"){
    load("data/pipeline/v17/bem/LGG_mut.ro")
    load("data/pipeline/v17/bem/LGG_methyl.ro")
    load("data/pipeline/v17/bem/LGG_cnv.ro")
    commonCell <- intersect(rownames(LGG_mut), intersect(rownames(LGG_methyl), rownames(LGG_cnv)))
    feature <- cbind(LGG_mut[commonCell,], LGG_methyl[commonCell,], LGG_cnv[commonCell,])
  }
  if (feat=="v17_LGG_mut.methyl.speed"){
    load("data/pipeline/v17/bem/LGG_mut.ro")
    load("data/pipeline/v17/bem/LGG_methyl.ro")
    load("data/pipeline/v17/bem/LGG_speed.ro")
    commonCell <- intersect(rownames(LGG_mut), intersect(rownames(LGG_methyl), rownames(LGG_speed)))
    feature <- cbind(LGG_mut[commonCell,], LGG_methyl[commonCell,], LGG_speed[commonCell,])
  }
  if (feat=="v17_LGG_mut.gex.cnv"){
    load("data/pipeline/v17/bem/LGG_mut.ro")
    load("data/pipeline/v17/bem/LGG_gex.ro")
    load("data/pipeline/v17/bem/LGG_cnv.ro")
    commonCell <- intersect(rownames(LGG_mut), intersect(rownames(LGG_gex), rownames(LGG_cnv)))
    feature <- cbind(LGG_mut[commonCell,], LGG_gex[commonCell,], LGG_cnv[commonCell,])
  }
  if (feat=="v17_LGG_mut.gex.speed"){
    load("data/pipeline/v17/bem/LGG_mut.ro")
    load("data/pipeline/v17/bem/LGG_gex.ro")
    load("data/pipeline/v17/bem/LGG_speed.ro")
    commonCell <- intersect(rownames(LGG_mut), intersect(rownames(LGG_gex), rownames(LGG_speed)))
    feature <- cbind(LGG_mut[commonCell,], LGG_gex[commonCell,], LGG_speed[commonCell,])
  }
  if (feat=="v17_LGG_mut.cnv.speed"){
    load("data/pipeline/v17/bem/LGG_mut.ro")
    load("data/pipeline/v17/bem/LGG_cnv.ro")
    load("data/pipeline/v17/bem/LGG_speed.ro")
    commonCell <- intersect(rownames(LGG_mut), intersect(rownames(LGG_cnv), rownames(LGG_speed)))
    feature <- cbind(LGG_mut[commonCell,], LGG_cnv[commonCell,], LGG_speed[commonCell,])
  }
  if (feat=="v17_LGG_methyl.gex.cnv"){
    load("data/pipeline/v17/bem/LGG_methyl.ro")
    load("data/pipeline/v17/bem/LGG_gex.ro")
    load("data/pipeline/v17/bem/LGG_cnv.ro")
    commonCell <- intersect(rownames(LGG_methyl), intersect(rownames(LGG_gex), rownames(LGG_cnv)))
    feature <- cbind(LGG_methyl[commonCell,], LGG_gex[commonCell,], LGG_cnv[commonCell,])
  }
  if (feat=="v17_LGG_methyl.gex.speed"){
    load("data/pipeline/v17/bem/LGG_methyl.ro")
    load("data/pipeline/v17/bem/LGG_gex.ro")
    load("data/pipeline/v17/bem/LGG_speed.ro")
    commonCell <- intersect(rownames(LGG_methyl), intersect(rownames(LGG_gex), rownames(LGG_speed)))
    feature <- cbind(LGG_methyl[commonCell,], LGG_gex[commonCell,], LGG_speed[commonCell,])
  }
  if (feat=="v17_LGG_methyl.cnv.speed"){
    load("data/pipeline/v17/bem/LGG_methyl.ro")
    load("data/pipeline/v17/bem/LGG_cnv.ro")
    load("data/pipeline/v17/bem/LGG_speed.ro")
    commonCell <- intersect(rownames(LGG_methyl), intersect(rownames(LGG_cnv), rownames(LGG_speed)))
    feature <- cbind(LGG_methyl[commonCell,], LGG_cnv[commonCell,], LGG_speed[commonCell,])
  }
  if (feat=="v17_LGG_gex.cnv.speed"){
    load("data/pipeline/v17/bem/LGG_gex.ro")
    load("data/pipeline/v17/bem/LGG_cnv.ro")
    load("data/pipeline/v17/bem/LGG_speed.ro")
    commonCell <- intersect(rownames(LGG_gex), intersect(rownames(LGG_cnv), rownames(LGG_speed)))
    feature <- cbind(LGG_gex[commonCell,], LGG_cnv[commonCell,], LGG_speed[commonCell,])
  }
  if (feat=="v17_LGG_mut.methyl.gex.cnv"){
    load("data/pipeline/v17/bem/LGG_mut.ro")
    load("data/pipeline/v17/bem/LGG_methyl.ro")
    load("data/pipeline/v17/bem/LGG_gex.ro")
    load("data/pipeline/v17/bem/LGG_cnv.ro")
    commonCell <- intersect(rownames(LGG_mut), intersect(rownames(LGG_methyl), intersect(rownames(LGG_gex), rownames(LGG_cnv))))
    feature <- cbind(LGG_mut[commonCell,], LGG_methyl[commonCell,], LGG_gex[commonCell,], LGG_cnv[commonCell,])
  }
  if (feat=="v17_LGG_mut.methyl.gex.speed"){
    load("data/pipeline/v17/bem/LGG_mut.ro")
    load("data/pipeline/v17/bem/LGG_methyl.ro")
    load("data/pipeline/v17/bem/LGG_gex.ro")
    load("data/pipeline/v17/bem/LGG_speed.ro")
    commonCell <- intersect(rownames(LGG_mut), intersect(rownames(LGG_methyl), intersect(rownames(LGG_gex), rownames(LGG_speed))))
    feature <- cbind(LGG_mut[commonCell,], LGG_methyl[commonCell,], LGG_gex[commonCell,], LGG_speed[commonCell,])
  }
  if (feat=="v17_LGG_mut.methyl.cnv.speed"){
    load("data/pipeline/v17/bem/LGG_mut.ro")
    load("data/pipeline/v17/bem/LGG_methyl.ro")
    load("data/pipeline/v17/bem/LGG_cnv.ro")
    load("data/pipeline/v17/bem/LGG_speed.ro")
    commonCell <- intersect(rownames(LGG_mut), intersect(rownames(LGG_methyl), intersect(rownames(LGG_cnv), rownames(LGG_speed))))
    feature <- cbind(LGG_mut[commonCell,], LGG_methyl[commonCell,], LGG_cnv[commonCell,], LGG_speed[commonCell,])
  }
  if (feat=="v17_LGG_mut.gex.cnv.speed"){
    load("data/pipeline/v17/bem/LGG_mut.ro")
    load("data/pipeline/v17/bem/LGG_gex.ro")
    load("data/pipeline/v17/bem/LGG_cnv.ro")
    load("data/pipeline/v17/bem/LGG_speed.ro")
    commonCell <- intersect(rownames(LGG_mut), intersect(rownames(LGG_gex), intersect(rownames(LGG_cnv), rownames(LGG_speed))))
    feature <- cbind(LGG_mut[commonCell,], LGG_gex[commonCell,], LGG_cnv[commonCell,], LGG_speed[commonCell,])
  }
  if (feat=="v17_LGG_methyl.gex.cnv.speed"){
    load("data/pipeline/v17/bem/LGG_methyl.ro")
    load("data/pipeline/v17/bem/LGG_gex.ro")
    load("data/pipeline/v17/bem/LGG_cnv.ro")
    load("data/pipeline/v17/bem/LGG_speed.ro")
    commonCell <- intersect(rownames(LGG_methyl), intersect(rownames(LGG_gex), intersect(rownames(LGG_cnv), rownames(LGG_speed))))
    feature <- cbind(LGG_methyl[commonCell,], LGG_gex[commonCell,], LGG_cnv[commonCell,], LGG_speed[commonCell,])
  }
  if (feat=="v17_LGG_mut.methyl.gex.cnv.speed"){
    load("data/pipeline/v17/bem/LGG_mut.ro")
    load("data/pipeline/v17/bem/LGG_methyl.ro")
    load("data/pipeline/v17/bem/LGG_gex.ro")
    load("data/pipeline/v17/bem/LGG_cnv.ro")
    load("data/pipeline/v17/bem/LGG_speed.ro")
    commonCell <- intersect(rownames(LGG_mut), intersect(rownames(LGG_methyl), intersect(rownames(LGG_gex), intersect(rownames(LGG_cnv), rownames(LGG_speed)))))
    feature <- cbind(LGG_mut[commonCell,], LGG_methyl[commonCell,], LGG_gex[commonCell,], LGG_cnv[commonCell,], LGG_speed[commonCell,])
  }
  if (feat=="v17_LIHC_mut"){
    load("data/pipeline/v17/bem/LIHC_mut.ro")
    commonCell <- rownames(LIHC_mut)
    feature <- cbind(LIHC_mut[commonCell,])
  }
  if (feat=="v17_LIHC_methyl"){
    load("data/pipeline/v17/bem/LIHC_methyl.ro")
    commonCell <- rownames(LIHC_methyl)
    feature <- cbind(LIHC_methyl[commonCell,])
  }
  if (feat=="v17_LIHC_gex"){
    load("data/pipeline/v17/bem/LIHC_gex.ro")
    commonCell <- rownames(LIHC_gex)
    feature <- cbind(LIHC_gex[commonCell,])
  }
  if (feat=="v17_LIHC_cnv"){
    load("data/pipeline/v17/bem/LIHC_cnv.ro")
    commonCell <- rownames(LIHC_cnv)
    feature <- cbind(LIHC_cnv[commonCell,])
  }
  if (feat=="v17_LIHC_speed"){
    load("data/pipeline/v17/bem/LIHC_speed.ro")
    commonCell <- rownames(LIHC_speed)
    feature <- cbind(LIHC_speed[commonCell,])
  }
  if (feat=="v17_LIHC_mut.methyl"){
    load("data/pipeline/v17/bem/LIHC_mut.ro")
    load("data/pipeline/v17/bem/LIHC_methyl.ro")
    commonCell <- intersect(rownames(LIHC_mut), rownames(LIHC_methyl))
    feature <- cbind(LIHC_mut[commonCell,], LIHC_methyl[commonCell,])
  }
  if (feat=="v17_LIHC_mut.gex"){
    load("data/pipeline/v17/bem/LIHC_mut.ro")
    load("data/pipeline/v17/bem/LIHC_gex.ro")
    commonCell <- intersect(rownames(LIHC_mut), rownames(LIHC_gex))
    feature <- cbind(LIHC_mut[commonCell,], LIHC_gex[commonCell,])
  }
  if (feat=="v17_LIHC_mut.cnv"){
    load("data/pipeline/v17/bem/LIHC_mut.ro")
    load("data/pipeline/v17/bem/LIHC_cnv.ro")
    commonCell <- intersect(rownames(LIHC_mut), rownames(LIHC_cnv))
    feature <- cbind(LIHC_mut[commonCell,], LIHC_cnv[commonCell,])
  }
  if (feat=="v17_LIHC_mut.speed"){
    load("data/pipeline/v17/bem/LIHC_mut.ro")
    load("data/pipeline/v17/bem/LIHC_speed.ro")
    commonCell <- intersect(rownames(LIHC_mut), rownames(LIHC_speed))
    feature <- cbind(LIHC_mut[commonCell,], LIHC_speed[commonCell,])
  }
  if (feat=="v17_LIHC_methyl.gex"){
    load("data/pipeline/v17/bem/LIHC_methyl.ro")
    load("data/pipeline/v17/bem/LIHC_gex.ro")
    commonCell <- intersect(rownames(LIHC_methyl), rownames(LIHC_gex))
    feature <- cbind(LIHC_methyl[commonCell,], LIHC_gex[commonCell,])
  }
  if (feat=="v17_LIHC_methyl.cnv"){
    load("data/pipeline/v17/bem/LIHC_methyl.ro")
    load("data/pipeline/v17/bem/LIHC_cnv.ro")
    commonCell <- intersect(rownames(LIHC_methyl), rownames(LIHC_cnv))
    feature <- cbind(LIHC_methyl[commonCell,], LIHC_cnv[commonCell,])
  }
  if (feat=="v17_LIHC_methyl.speed"){
    load("data/pipeline/v17/bem/LIHC_methyl.ro")
    load("data/pipeline/v17/bem/LIHC_speed.ro")
    commonCell <- intersect(rownames(LIHC_methyl), rownames(LIHC_speed))
    feature <- cbind(LIHC_methyl[commonCell,], LIHC_speed[commonCell,])
  }
  if (feat=="v17_LIHC_gex.cnv"){
    load("data/pipeline/v17/bem/LIHC_gex.ro")
    load("data/pipeline/v17/bem/LIHC_cnv.ro")
    commonCell <- intersect(rownames(LIHC_gex), rownames(LIHC_cnv))
    feature <- cbind(LIHC_gex[commonCell,], LIHC_cnv[commonCell,])
  }
  if (feat=="v17_LIHC_gex.speed"){
    load("data/pipeline/v17/bem/LIHC_gex.ro")
    load("data/pipeline/v17/bem/LIHC_speed.ro")
    commonCell <- intersect(rownames(LIHC_gex), rownames(LIHC_speed))
    feature <- cbind(LIHC_gex[commonCell,], LIHC_speed[commonCell,])
  }
  if (feat=="v17_LIHC_cnv.speed"){
    load("data/pipeline/v17/bem/LIHC_cnv.ro")
    load("data/pipeline/v17/bem/LIHC_speed.ro")
    commonCell <- intersect(rownames(LIHC_cnv), rownames(LIHC_speed))
    feature <- cbind(LIHC_cnv[commonCell,], LIHC_speed[commonCell,])
  }
  if (feat=="v17_LIHC_mut.methyl.gex"){
    load("data/pipeline/v17/bem/LIHC_mut.ro")
    load("data/pipeline/v17/bem/LIHC_methyl.ro")
    load("data/pipeline/v17/bem/LIHC_gex.ro")
    commonCell <- intersect(rownames(LIHC_mut), intersect(rownames(LIHC_methyl), rownames(LIHC_gex)))
    feature <- cbind(LIHC_mut[commonCell,], LIHC_methyl[commonCell,], LIHC_gex[commonCell,])
  }
  if (feat=="v17_LIHC_mut.methyl.cnv"){
    load("data/pipeline/v17/bem/LIHC_mut.ro")
    load("data/pipeline/v17/bem/LIHC_methyl.ro")
    load("data/pipeline/v17/bem/LIHC_cnv.ro")
    commonCell <- intersect(rownames(LIHC_mut), intersect(rownames(LIHC_methyl), rownames(LIHC_cnv)))
    feature <- cbind(LIHC_mut[commonCell,], LIHC_methyl[commonCell,], LIHC_cnv[commonCell,])
  }
  if (feat=="v17_LIHC_mut.methyl.speed"){
    load("data/pipeline/v17/bem/LIHC_mut.ro")
    load("data/pipeline/v17/bem/LIHC_methyl.ro")
    load("data/pipeline/v17/bem/LIHC_speed.ro")
    commonCell <- intersect(rownames(LIHC_mut), intersect(rownames(LIHC_methyl), rownames(LIHC_speed)))
    feature <- cbind(LIHC_mut[commonCell,], LIHC_methyl[commonCell,], LIHC_speed[commonCell,])
  }
  if (feat=="v17_LIHC_mut.gex.cnv"){
    load("data/pipeline/v17/bem/LIHC_mut.ro")
    load("data/pipeline/v17/bem/LIHC_gex.ro")
    load("data/pipeline/v17/bem/LIHC_cnv.ro")
    commonCell <- intersect(rownames(LIHC_mut), intersect(rownames(LIHC_gex), rownames(LIHC_cnv)))
    feature <- cbind(LIHC_mut[commonCell,], LIHC_gex[commonCell,], LIHC_cnv[commonCell,])
  }
  if (feat=="v17_LIHC_mut.gex.speed"){
    load("data/pipeline/v17/bem/LIHC_mut.ro")
    load("data/pipeline/v17/bem/LIHC_gex.ro")
    load("data/pipeline/v17/bem/LIHC_speed.ro")
    commonCell <- intersect(rownames(LIHC_mut), intersect(rownames(LIHC_gex), rownames(LIHC_speed)))
    feature <- cbind(LIHC_mut[commonCell,], LIHC_gex[commonCell,], LIHC_speed[commonCell,])
  }
  if (feat=="v17_LIHC_mut.cnv.speed"){
    load("data/pipeline/v17/bem/LIHC_mut.ro")
    load("data/pipeline/v17/bem/LIHC_cnv.ro")
    load("data/pipeline/v17/bem/LIHC_speed.ro")
    commonCell <- intersect(rownames(LIHC_mut), intersect(rownames(LIHC_cnv), rownames(LIHC_speed)))
    feature <- cbind(LIHC_mut[commonCell,], LIHC_cnv[commonCell,], LIHC_speed[commonCell,])
  }
  if (feat=="v17_LIHC_methyl.gex.cnv"){
    load("data/pipeline/v17/bem/LIHC_methyl.ro")
    load("data/pipeline/v17/bem/LIHC_gex.ro")
    load("data/pipeline/v17/bem/LIHC_cnv.ro")
    commonCell <- intersect(rownames(LIHC_methyl), intersect(rownames(LIHC_gex), rownames(LIHC_cnv)))
    feature <- cbind(LIHC_methyl[commonCell,], LIHC_gex[commonCell,], LIHC_cnv[commonCell,])
  }
  if (feat=="v17_LIHC_methyl.gex.speed"){
    load("data/pipeline/v17/bem/LIHC_methyl.ro")
    load("data/pipeline/v17/bem/LIHC_gex.ro")
    load("data/pipeline/v17/bem/LIHC_speed.ro")
    commonCell <- intersect(rownames(LIHC_methyl), intersect(rownames(LIHC_gex), rownames(LIHC_speed)))
    feature <- cbind(LIHC_methyl[commonCell,], LIHC_gex[commonCell,], LIHC_speed[commonCell,])
  }
  if (feat=="v17_LIHC_methyl.cnv.speed"){
    load("data/pipeline/v17/bem/LIHC_methyl.ro")
    load("data/pipeline/v17/bem/LIHC_cnv.ro")
    load("data/pipeline/v17/bem/LIHC_speed.ro")
    commonCell <- intersect(rownames(LIHC_methyl), intersect(rownames(LIHC_cnv), rownames(LIHC_speed)))
    feature <- cbind(LIHC_methyl[commonCell,], LIHC_cnv[commonCell,], LIHC_speed[commonCell,])
  }
  if (feat=="v17_LIHC_gex.cnv.speed"){
    load("data/pipeline/v17/bem/LIHC_gex.ro")
    load("data/pipeline/v17/bem/LIHC_cnv.ro")
    load("data/pipeline/v17/bem/LIHC_speed.ro")
    commonCell <- intersect(rownames(LIHC_gex), intersect(rownames(LIHC_cnv), rownames(LIHC_speed)))
    feature <- cbind(LIHC_gex[commonCell,], LIHC_cnv[commonCell,], LIHC_speed[commonCell,])
  }
  if (feat=="v17_LIHC_mut.methyl.gex.cnv"){
    load("data/pipeline/v17/bem/LIHC_mut.ro")
    load("data/pipeline/v17/bem/LIHC_methyl.ro")
    load("data/pipeline/v17/bem/LIHC_gex.ro")
    load("data/pipeline/v17/bem/LIHC_cnv.ro")
    commonCell <- intersect(rownames(LIHC_mut), intersect(rownames(LIHC_methyl), intersect(rownames(LIHC_gex), rownames(LIHC_cnv))))
    feature <- cbind(LIHC_mut[commonCell,], LIHC_methyl[commonCell,], LIHC_gex[commonCell,], LIHC_cnv[commonCell,])
  }
  if (feat=="v17_LIHC_mut.methyl.gex.speed"){
    load("data/pipeline/v17/bem/LIHC_mut.ro")
    load("data/pipeline/v17/bem/LIHC_methyl.ro")
    load("data/pipeline/v17/bem/LIHC_gex.ro")
    load("data/pipeline/v17/bem/LIHC_speed.ro")
    commonCell <- intersect(rownames(LIHC_mut), intersect(rownames(LIHC_methyl), intersect(rownames(LIHC_gex), rownames(LIHC_speed))))
    feature <- cbind(LIHC_mut[commonCell,], LIHC_methyl[commonCell,], LIHC_gex[commonCell,], LIHC_speed[commonCell,])
  }
  if (feat=="v17_LIHC_mut.methyl.cnv.speed"){
    load("data/pipeline/v17/bem/LIHC_mut.ro")
    load("data/pipeline/v17/bem/LIHC_methyl.ro")
    load("data/pipeline/v17/bem/LIHC_cnv.ro")
    load("data/pipeline/v17/bem/LIHC_speed.ro")
    commonCell <- intersect(rownames(LIHC_mut), intersect(rownames(LIHC_methyl), intersect(rownames(LIHC_cnv), rownames(LIHC_speed))))
    feature <- cbind(LIHC_mut[commonCell,], LIHC_methyl[commonCell,], LIHC_cnv[commonCell,], LIHC_speed[commonCell,])
  }
  if (feat=="v17_LIHC_mut.gex.cnv.speed"){
    load("data/pipeline/v17/bem/LIHC_mut.ro")
    load("data/pipeline/v17/bem/LIHC_gex.ro")
    load("data/pipeline/v17/bem/LIHC_cnv.ro")
    load("data/pipeline/v17/bem/LIHC_speed.ro")
    commonCell <- intersect(rownames(LIHC_mut), intersect(rownames(LIHC_gex), intersect(rownames(LIHC_cnv), rownames(LIHC_speed))))
    feature <- cbind(LIHC_mut[commonCell,], LIHC_gex[commonCell,], LIHC_cnv[commonCell,], LIHC_speed[commonCell,])
  }
  if (feat=="v17_LIHC_methyl.gex.cnv.speed"){
    load("data/pipeline/v17/bem/LIHC_methyl.ro")
    load("data/pipeline/v17/bem/LIHC_gex.ro")
    load("data/pipeline/v17/bem/LIHC_cnv.ro")
    load("data/pipeline/v17/bem/LIHC_speed.ro")
    commonCell <- intersect(rownames(LIHC_methyl), intersect(rownames(LIHC_gex), intersect(rownames(LIHC_cnv), rownames(LIHC_speed))))
    feature <- cbind(LIHC_methyl[commonCell,], LIHC_gex[commonCell,], LIHC_cnv[commonCell,], LIHC_speed[commonCell,])
  }
  if (feat=="v17_LIHC_mut.methyl.gex.cnv.speed"){
    load("data/pipeline/v17/bem/LIHC_mut.ro")
    load("data/pipeline/v17/bem/LIHC_methyl.ro")
    load("data/pipeline/v17/bem/LIHC_gex.ro")
    load("data/pipeline/v17/bem/LIHC_cnv.ro")
    load("data/pipeline/v17/bem/LIHC_speed.ro")
    commonCell <- intersect(rownames(LIHC_mut), intersect(rownames(LIHC_methyl), intersect(rownames(LIHC_gex), intersect(rownames(LIHC_cnv), rownames(LIHC_speed)))))
    feature <- cbind(LIHC_mut[commonCell,], LIHC_methyl[commonCell,], LIHC_gex[commonCell,], LIHC_cnv[commonCell,], LIHC_speed[commonCell,])
  }
  if (feat=="v17_LUAD_mut"){
    load("data/pipeline/v17/bem/LUAD_mut.ro")
    commonCell <- rownames(LUAD_mut)
    feature <- cbind(LUAD_mut[commonCell,])
  }
  if (feat=="v17_LUAD_methyl"){
    load("data/pipeline/v17/bem/LUAD_methyl.ro")
    commonCell <- rownames(LUAD_methyl)
    feature <- cbind(LUAD_methyl[commonCell,])
  }
  if (feat=="v17_LUAD_gex"){
    load("data/pipeline/v17/bem/LUAD_gex.ro")
    commonCell <- rownames(LUAD_gex)
    feature <- cbind(LUAD_gex[commonCell,])
  }
  if (feat=="v17_LUAD_cnv"){
    load("data/pipeline/v17/bem/LUAD_cnv.ro")
    commonCell <- rownames(LUAD_cnv)
    feature <- cbind(LUAD_cnv[commonCell,])
  }
  if (feat=="v17_LUAD_speed"){
    load("data/pipeline/v17/bem/LUAD_speed.ro")
    commonCell <- rownames(LUAD_speed)
    feature <- cbind(LUAD_speed[commonCell,])
  }
  if (feat=="v17_LUAD_mut.methyl"){
    load("data/pipeline/v17/bem/LUAD_mut.ro")
    load("data/pipeline/v17/bem/LUAD_methyl.ro")
    commonCell <- intersect(rownames(LUAD_mut), rownames(LUAD_methyl))
    feature <- cbind(LUAD_mut[commonCell,], LUAD_methyl[commonCell,])
  }
  if (feat=="v17_LUAD_mut.gex"){
    load("data/pipeline/v17/bem/LUAD_mut.ro")
    load("data/pipeline/v17/bem/LUAD_gex.ro")
    commonCell <- intersect(rownames(LUAD_mut), rownames(LUAD_gex))
    feature <- cbind(LUAD_mut[commonCell,], LUAD_gex[commonCell,])
  }
  if (feat=="v17_LUAD_mut.cnv"){
    load("data/pipeline/v17/bem/LUAD_mut.ro")
    load("data/pipeline/v17/bem/LUAD_cnv.ro")
    commonCell <- intersect(rownames(LUAD_mut), rownames(LUAD_cnv))
    feature <- cbind(LUAD_mut[commonCell,], LUAD_cnv[commonCell,])
  }
  if (feat=="v17_LUAD_mut.speed"){
    load("data/pipeline/v17/bem/LUAD_mut.ro")
    load("data/pipeline/v17/bem/LUAD_speed.ro")
    commonCell <- intersect(rownames(LUAD_mut), rownames(LUAD_speed))
    feature <- cbind(LUAD_mut[commonCell,], LUAD_speed[commonCell,])
  }
  if (feat=="v17_LUAD_methyl.gex"){
    load("data/pipeline/v17/bem/LUAD_methyl.ro")
    load("data/pipeline/v17/bem/LUAD_gex.ro")
    commonCell <- intersect(rownames(LUAD_methyl), rownames(LUAD_gex))
    feature <- cbind(LUAD_methyl[commonCell,], LUAD_gex[commonCell,])
  }
  if (feat=="v17_LUAD_methyl.cnv"){
    load("data/pipeline/v17/bem/LUAD_methyl.ro")
    load("data/pipeline/v17/bem/LUAD_cnv.ro")
    commonCell <- intersect(rownames(LUAD_methyl), rownames(LUAD_cnv))
    feature <- cbind(LUAD_methyl[commonCell,], LUAD_cnv[commonCell,])
  }
  if (feat=="v17_LUAD_methyl.speed"){
    load("data/pipeline/v17/bem/LUAD_methyl.ro")
    load("data/pipeline/v17/bem/LUAD_speed.ro")
    commonCell <- intersect(rownames(LUAD_methyl), rownames(LUAD_speed))
    feature <- cbind(LUAD_methyl[commonCell,], LUAD_speed[commonCell,])
  }
  if (feat=="v17_LUAD_gex.cnv"){
    load("data/pipeline/v17/bem/LUAD_gex.ro")
    load("data/pipeline/v17/bem/LUAD_cnv.ro")
    commonCell <- intersect(rownames(LUAD_gex), rownames(LUAD_cnv))
    feature <- cbind(LUAD_gex[commonCell,], LUAD_cnv[commonCell,])
  }
  if (feat=="v17_LUAD_gex.speed"){
    load("data/pipeline/v17/bem/LUAD_gex.ro")
    load("data/pipeline/v17/bem/LUAD_speed.ro")
    commonCell <- intersect(rownames(LUAD_gex), rownames(LUAD_speed))
    feature <- cbind(LUAD_gex[commonCell,], LUAD_speed[commonCell,])
  }
  if (feat=="v17_LUAD_cnv.speed"){
    load("data/pipeline/v17/bem/LUAD_cnv.ro")
    load("data/pipeline/v17/bem/LUAD_speed.ro")
    commonCell <- intersect(rownames(LUAD_cnv), rownames(LUAD_speed))
    feature <- cbind(LUAD_cnv[commonCell,], LUAD_speed[commonCell,])
  }
  if (feat=="v17_LUAD_mut.methyl.gex"){
    load("data/pipeline/v17/bem/LUAD_mut.ro")
    load("data/pipeline/v17/bem/LUAD_methyl.ro")
    load("data/pipeline/v17/bem/LUAD_gex.ro")
    commonCell <- intersect(rownames(LUAD_mut), intersect(rownames(LUAD_methyl), rownames(LUAD_gex)))
    feature <- cbind(LUAD_mut[commonCell,], LUAD_methyl[commonCell,], LUAD_gex[commonCell,])
  }
  if (feat=="v17_LUAD_mut.methyl.cnv"){
    load("data/pipeline/v17/bem/LUAD_mut.ro")
    load("data/pipeline/v17/bem/LUAD_methyl.ro")
    load("data/pipeline/v17/bem/LUAD_cnv.ro")
    commonCell <- intersect(rownames(LUAD_mut), intersect(rownames(LUAD_methyl), rownames(LUAD_cnv)))
    feature <- cbind(LUAD_mut[commonCell,], LUAD_methyl[commonCell,], LUAD_cnv[commonCell,])
  }
  if (feat=="v17_LUAD_mut.methyl.speed"){
    load("data/pipeline/v17/bem/LUAD_mut.ro")
    load("data/pipeline/v17/bem/LUAD_methyl.ro")
    load("data/pipeline/v17/bem/LUAD_speed.ro")
    commonCell <- intersect(rownames(LUAD_mut), intersect(rownames(LUAD_methyl), rownames(LUAD_speed)))
    feature <- cbind(LUAD_mut[commonCell,], LUAD_methyl[commonCell,], LUAD_speed[commonCell,])
  }
  if (feat=="v17_LUAD_mut.gex.cnv"){
    load("data/pipeline/v17/bem/LUAD_mut.ro")
    load("data/pipeline/v17/bem/LUAD_gex.ro")
    load("data/pipeline/v17/bem/LUAD_cnv.ro")
    commonCell <- intersect(rownames(LUAD_mut), intersect(rownames(LUAD_gex), rownames(LUAD_cnv)))
    feature <- cbind(LUAD_mut[commonCell,], LUAD_gex[commonCell,], LUAD_cnv[commonCell,])
  }
  if (feat=="v17_LUAD_mut.gex.speed"){
    load("data/pipeline/v17/bem/LUAD_mut.ro")
    load("data/pipeline/v17/bem/LUAD_gex.ro")
    load("data/pipeline/v17/bem/LUAD_speed.ro")
    commonCell <- intersect(rownames(LUAD_mut), intersect(rownames(LUAD_gex), rownames(LUAD_speed)))
    feature <- cbind(LUAD_mut[commonCell,], LUAD_gex[commonCell,], LUAD_speed[commonCell,])
  }
  if (feat=="v17_LUAD_mut.cnv.speed"){
    load("data/pipeline/v17/bem/LUAD_mut.ro")
    load("data/pipeline/v17/bem/LUAD_cnv.ro")
    load("data/pipeline/v17/bem/LUAD_speed.ro")
    commonCell <- intersect(rownames(LUAD_mut), intersect(rownames(LUAD_cnv), rownames(LUAD_speed)))
    feature <- cbind(LUAD_mut[commonCell,], LUAD_cnv[commonCell,], LUAD_speed[commonCell,])
  }
  if (feat=="v17_LUAD_methyl.gex.cnv"){
    load("data/pipeline/v17/bem/LUAD_methyl.ro")
    load("data/pipeline/v17/bem/LUAD_gex.ro")
    load("data/pipeline/v17/bem/LUAD_cnv.ro")
    commonCell <- intersect(rownames(LUAD_methyl), intersect(rownames(LUAD_gex), rownames(LUAD_cnv)))
    feature <- cbind(LUAD_methyl[commonCell,], LUAD_gex[commonCell,], LUAD_cnv[commonCell,])
  }
  if (feat=="v17_LUAD_methyl.gex.speed"){
    load("data/pipeline/v17/bem/LUAD_methyl.ro")
    load("data/pipeline/v17/bem/LUAD_gex.ro")
    load("data/pipeline/v17/bem/LUAD_speed.ro")
    commonCell <- intersect(rownames(LUAD_methyl), intersect(rownames(LUAD_gex), rownames(LUAD_speed)))
    feature <- cbind(LUAD_methyl[commonCell,], LUAD_gex[commonCell,], LUAD_speed[commonCell,])
  }
  if (feat=="v17_LUAD_methyl.cnv.speed"){
    load("data/pipeline/v17/bem/LUAD_methyl.ro")
    load("data/pipeline/v17/bem/LUAD_cnv.ro")
    load("data/pipeline/v17/bem/LUAD_speed.ro")
    commonCell <- intersect(rownames(LUAD_methyl), intersect(rownames(LUAD_cnv), rownames(LUAD_speed)))
    feature <- cbind(LUAD_methyl[commonCell,], LUAD_cnv[commonCell,], LUAD_speed[commonCell,])
  }
  if (feat=="v17_LUAD_gex.cnv.speed"){
    load("data/pipeline/v17/bem/LUAD_gex.ro")
    load("data/pipeline/v17/bem/LUAD_cnv.ro")
    load("data/pipeline/v17/bem/LUAD_speed.ro")
    commonCell <- intersect(rownames(LUAD_gex), intersect(rownames(LUAD_cnv), rownames(LUAD_speed)))
    feature <- cbind(LUAD_gex[commonCell,], LUAD_cnv[commonCell,], LUAD_speed[commonCell,])
  }
  if (feat=="v17_LUAD_mut.methyl.gex.cnv"){
    load("data/pipeline/v17/bem/LUAD_mut.ro")
    load("data/pipeline/v17/bem/LUAD_methyl.ro")
    load("data/pipeline/v17/bem/LUAD_gex.ro")
    load("data/pipeline/v17/bem/LUAD_cnv.ro")
    commonCell <- intersect(rownames(LUAD_mut), intersect(rownames(LUAD_methyl), intersect(rownames(LUAD_gex), rownames(LUAD_cnv))))
    feature <- cbind(LUAD_mut[commonCell,], LUAD_methyl[commonCell,], LUAD_gex[commonCell,], LUAD_cnv[commonCell,])
  }
  if (feat=="v17_LUAD_mut.methyl.gex.speed"){
    load("data/pipeline/v17/bem/LUAD_mut.ro")
    load("data/pipeline/v17/bem/LUAD_methyl.ro")
    load("data/pipeline/v17/bem/LUAD_gex.ro")
    load("data/pipeline/v17/bem/LUAD_speed.ro")
    commonCell <- intersect(rownames(LUAD_mut), intersect(rownames(LUAD_methyl), intersect(rownames(LUAD_gex), rownames(LUAD_speed))))
    feature <- cbind(LUAD_mut[commonCell,], LUAD_methyl[commonCell,], LUAD_gex[commonCell,], LUAD_speed[commonCell,])
  }
  if (feat=="v17_LUAD_mut.methyl.cnv.speed"){
    load("data/pipeline/v17/bem/LUAD_mut.ro")
    load("data/pipeline/v17/bem/LUAD_methyl.ro")
    load("data/pipeline/v17/bem/LUAD_cnv.ro")
    load("data/pipeline/v17/bem/LUAD_speed.ro")
    commonCell <- intersect(rownames(LUAD_mut), intersect(rownames(LUAD_methyl), intersect(rownames(LUAD_cnv), rownames(LUAD_speed))))
    feature <- cbind(LUAD_mut[commonCell,], LUAD_methyl[commonCell,], LUAD_cnv[commonCell,], LUAD_speed[commonCell,])
  }
  if (feat=="v17_LUAD_mut.gex.cnv.speed"){
    load("data/pipeline/v17/bem/LUAD_mut.ro")
    load("data/pipeline/v17/bem/LUAD_gex.ro")
    load("data/pipeline/v17/bem/LUAD_cnv.ro")
    load("data/pipeline/v17/bem/LUAD_speed.ro")
    commonCell <- intersect(rownames(LUAD_mut), intersect(rownames(LUAD_gex), intersect(rownames(LUAD_cnv), rownames(LUAD_speed))))
    feature <- cbind(LUAD_mut[commonCell,], LUAD_gex[commonCell,], LUAD_cnv[commonCell,], LUAD_speed[commonCell,])
  }
  if (feat=="v17_LUAD_methyl.gex.cnv.speed"){
    load("data/pipeline/v17/bem/LUAD_methyl.ro")
    load("data/pipeline/v17/bem/LUAD_gex.ro")
    load("data/pipeline/v17/bem/LUAD_cnv.ro")
    load("data/pipeline/v17/bem/LUAD_speed.ro")
    commonCell <- intersect(rownames(LUAD_methyl), intersect(rownames(LUAD_gex), intersect(rownames(LUAD_cnv), rownames(LUAD_speed))))
    feature <- cbind(LUAD_methyl[commonCell,], LUAD_gex[commonCell,], LUAD_cnv[commonCell,], LUAD_speed[commonCell,])
  }
  if (feat=="v17_LUAD_mut.methyl.gex.cnv.speed"){
    load("data/pipeline/v17/bem/LUAD_mut.ro")
    load("data/pipeline/v17/bem/LUAD_methyl.ro")
    load("data/pipeline/v17/bem/LUAD_gex.ro")
    load("data/pipeline/v17/bem/LUAD_cnv.ro")
    load("data/pipeline/v17/bem/LUAD_speed.ro")
    commonCell <- intersect(rownames(LUAD_mut), intersect(rownames(LUAD_methyl), intersect(rownames(LUAD_gex), intersect(rownames(LUAD_cnv), rownames(LUAD_speed)))))
    feature <- cbind(LUAD_mut[commonCell,], LUAD_methyl[commonCell,], LUAD_gex[commonCell,], LUAD_cnv[commonCell,], LUAD_speed[commonCell,])
  }
  if (feat=="v17_LUSC_mut"){
    load("data/pipeline/v17/bem/LUSC_mut.ro")
    commonCell <- rownames(LUSC_mut)
    feature <- cbind(LUSC_mut[commonCell,])
  }
  if (feat=="v17_LUSC_methyl"){
    load("data/pipeline/v17/bem/LUSC_methyl.ro")
    commonCell <- rownames(LUSC_methyl)
    feature <- cbind(LUSC_methyl[commonCell,])
  }
  if (feat=="v17_LUSC_gex"){
    load("data/pipeline/v17/bem/LUSC_gex.ro")
    commonCell <- rownames(LUSC_gex)
    feature <- cbind(LUSC_gex[commonCell,])
  }
  if (feat=="v17_LUSC_cnv"){
    load("data/pipeline/v17/bem/LUSC_cnv.ro")
    commonCell <- rownames(LUSC_cnv)
    feature <- cbind(LUSC_cnv[commonCell,])
  }
  if (feat=="v17_LUSC_speed"){
    load("data/pipeline/v17/bem/LUSC_speed.ro")
    commonCell <- rownames(LUSC_speed)
    feature <- cbind(LUSC_speed[commonCell,])
  }
  if (feat=="v17_LUSC_mut.methyl"){
    load("data/pipeline/v17/bem/LUSC_mut.ro")
    load("data/pipeline/v17/bem/LUSC_methyl.ro")
    commonCell <- intersect(rownames(LUSC_mut), rownames(LUSC_methyl))
    feature <- cbind(LUSC_mut[commonCell,], LUSC_methyl[commonCell,])
  }
  if (feat=="v17_LUSC_mut.gex"){
    load("data/pipeline/v17/bem/LUSC_mut.ro")
    load("data/pipeline/v17/bem/LUSC_gex.ro")
    commonCell <- intersect(rownames(LUSC_mut), rownames(LUSC_gex))
    feature <- cbind(LUSC_mut[commonCell,], LUSC_gex[commonCell,])
  }
  if (feat=="v17_LUSC_mut.cnv"){
    load("data/pipeline/v17/bem/LUSC_mut.ro")
    load("data/pipeline/v17/bem/LUSC_cnv.ro")
    commonCell <- intersect(rownames(LUSC_mut), rownames(LUSC_cnv))
    feature <- cbind(LUSC_mut[commonCell,], LUSC_cnv[commonCell,])
  }
  if (feat=="v17_LUSC_mut.speed"){
    load("data/pipeline/v17/bem/LUSC_mut.ro")
    load("data/pipeline/v17/bem/LUSC_speed.ro")
    commonCell <- intersect(rownames(LUSC_mut), rownames(LUSC_speed))
    feature <- cbind(LUSC_mut[commonCell,], LUSC_speed[commonCell,])
  }
  if (feat=="v17_LUSC_methyl.gex"){
    load("data/pipeline/v17/bem/LUSC_methyl.ro")
    load("data/pipeline/v17/bem/LUSC_gex.ro")
    commonCell <- intersect(rownames(LUSC_methyl), rownames(LUSC_gex))
    feature <- cbind(LUSC_methyl[commonCell,], LUSC_gex[commonCell,])
  }
  if (feat=="v17_LUSC_methyl.cnv"){
    load("data/pipeline/v17/bem/LUSC_methyl.ro")
    load("data/pipeline/v17/bem/LUSC_cnv.ro")
    commonCell <- intersect(rownames(LUSC_methyl), rownames(LUSC_cnv))
    feature <- cbind(LUSC_methyl[commonCell,], LUSC_cnv[commonCell,])
  }
  if (feat=="v17_LUSC_methyl.speed"){
    load("data/pipeline/v17/bem/LUSC_methyl.ro")
    load("data/pipeline/v17/bem/LUSC_speed.ro")
    commonCell <- intersect(rownames(LUSC_methyl), rownames(LUSC_speed))
    feature <- cbind(LUSC_methyl[commonCell,], LUSC_speed[commonCell,])
  }
  if (feat=="v17_LUSC_gex.cnv"){
    load("data/pipeline/v17/bem/LUSC_gex.ro")
    load("data/pipeline/v17/bem/LUSC_cnv.ro")
    commonCell <- intersect(rownames(LUSC_gex), rownames(LUSC_cnv))
    feature <- cbind(LUSC_gex[commonCell,], LUSC_cnv[commonCell,])
  }
  if (feat=="v17_LUSC_gex.speed"){
    load("data/pipeline/v17/bem/LUSC_gex.ro")
    load("data/pipeline/v17/bem/LUSC_speed.ro")
    commonCell <- intersect(rownames(LUSC_gex), rownames(LUSC_speed))
    feature <- cbind(LUSC_gex[commonCell,], LUSC_speed[commonCell,])
  }
  if (feat=="v17_LUSC_cnv.speed"){
    load("data/pipeline/v17/bem/LUSC_cnv.ro")
    load("data/pipeline/v17/bem/LUSC_speed.ro")
    commonCell <- intersect(rownames(LUSC_cnv), rownames(LUSC_speed))
    feature <- cbind(LUSC_cnv[commonCell,], LUSC_speed[commonCell,])
  }
  if (feat=="v17_LUSC_mut.methyl.gex"){
    load("data/pipeline/v17/bem/LUSC_mut.ro")
    load("data/pipeline/v17/bem/LUSC_methyl.ro")
    load("data/pipeline/v17/bem/LUSC_gex.ro")
    commonCell <- intersect(rownames(LUSC_mut), intersect(rownames(LUSC_methyl), rownames(LUSC_gex)))
    feature <- cbind(LUSC_mut[commonCell,], LUSC_methyl[commonCell,], LUSC_gex[commonCell,])
  }
  if (feat=="v17_LUSC_mut.methyl.cnv"){
    load("data/pipeline/v17/bem/LUSC_mut.ro")
    load("data/pipeline/v17/bem/LUSC_methyl.ro")
    load("data/pipeline/v17/bem/LUSC_cnv.ro")
    commonCell <- intersect(rownames(LUSC_mut), intersect(rownames(LUSC_methyl), rownames(LUSC_cnv)))
    feature <- cbind(LUSC_mut[commonCell,], LUSC_methyl[commonCell,], LUSC_cnv[commonCell,])
  }
  if (feat=="v17_LUSC_mut.methyl.speed"){
    load("data/pipeline/v17/bem/LUSC_mut.ro")
    load("data/pipeline/v17/bem/LUSC_methyl.ro")
    load("data/pipeline/v17/bem/LUSC_speed.ro")
    commonCell <- intersect(rownames(LUSC_mut), intersect(rownames(LUSC_methyl), rownames(LUSC_speed)))
    feature <- cbind(LUSC_mut[commonCell,], LUSC_methyl[commonCell,], LUSC_speed[commonCell,])
  }
  if (feat=="v17_LUSC_mut.gex.cnv"){
    load("data/pipeline/v17/bem/LUSC_mut.ro")
    load("data/pipeline/v17/bem/LUSC_gex.ro")
    load("data/pipeline/v17/bem/LUSC_cnv.ro")
    commonCell <- intersect(rownames(LUSC_mut), intersect(rownames(LUSC_gex), rownames(LUSC_cnv)))
    feature <- cbind(LUSC_mut[commonCell,], LUSC_gex[commonCell,], LUSC_cnv[commonCell,])
  }
  if (feat=="v17_LUSC_mut.gex.speed"){
    load("data/pipeline/v17/bem/LUSC_mut.ro")
    load("data/pipeline/v17/bem/LUSC_gex.ro")
    load("data/pipeline/v17/bem/LUSC_speed.ro")
    commonCell <- intersect(rownames(LUSC_mut), intersect(rownames(LUSC_gex), rownames(LUSC_speed)))
    feature <- cbind(LUSC_mut[commonCell,], LUSC_gex[commonCell,], LUSC_speed[commonCell,])
  }
  if (feat=="v17_LUSC_mut.cnv.speed"){
    load("data/pipeline/v17/bem/LUSC_mut.ro")
    load("data/pipeline/v17/bem/LUSC_cnv.ro")
    load("data/pipeline/v17/bem/LUSC_speed.ro")
    commonCell <- intersect(rownames(LUSC_mut), intersect(rownames(LUSC_cnv), rownames(LUSC_speed)))
    feature <- cbind(LUSC_mut[commonCell,], LUSC_cnv[commonCell,], LUSC_speed[commonCell,])
  }
  if (feat=="v17_LUSC_methyl.gex.cnv"){
    load("data/pipeline/v17/bem/LUSC_methyl.ro")
    load("data/pipeline/v17/bem/LUSC_gex.ro")
    load("data/pipeline/v17/bem/LUSC_cnv.ro")
    commonCell <- intersect(rownames(LUSC_methyl), intersect(rownames(LUSC_gex), rownames(LUSC_cnv)))
    feature <- cbind(LUSC_methyl[commonCell,], LUSC_gex[commonCell,], LUSC_cnv[commonCell,])
  }
  if (feat=="v17_LUSC_methyl.gex.speed"){
    load("data/pipeline/v17/bem/LUSC_methyl.ro")
    load("data/pipeline/v17/bem/LUSC_gex.ro")
    load("data/pipeline/v17/bem/LUSC_speed.ro")
    commonCell <- intersect(rownames(LUSC_methyl), intersect(rownames(LUSC_gex), rownames(LUSC_speed)))
    feature <- cbind(LUSC_methyl[commonCell,], LUSC_gex[commonCell,], LUSC_speed[commonCell,])
  }
  if (feat=="v17_LUSC_methyl.cnv.speed"){
    load("data/pipeline/v17/bem/LUSC_methyl.ro")
    load("data/pipeline/v17/bem/LUSC_cnv.ro")
    load("data/pipeline/v17/bem/LUSC_speed.ro")
    commonCell <- intersect(rownames(LUSC_methyl), intersect(rownames(LUSC_cnv), rownames(LUSC_speed)))
    feature <- cbind(LUSC_methyl[commonCell,], LUSC_cnv[commonCell,], LUSC_speed[commonCell,])
  }
  if (feat=="v17_LUSC_gex.cnv.speed"){
    load("data/pipeline/v17/bem/LUSC_gex.ro")
    load("data/pipeline/v17/bem/LUSC_cnv.ro")
    load("data/pipeline/v17/bem/LUSC_speed.ro")
    commonCell <- intersect(rownames(LUSC_gex), intersect(rownames(LUSC_cnv), rownames(LUSC_speed)))
    feature <- cbind(LUSC_gex[commonCell,], LUSC_cnv[commonCell,], LUSC_speed[commonCell,])
  }
  if (feat=="v17_LUSC_mut.methyl.gex.cnv"){
    load("data/pipeline/v17/bem/LUSC_mut.ro")
    load("data/pipeline/v17/bem/LUSC_methyl.ro")
    load("data/pipeline/v17/bem/LUSC_gex.ro")
    load("data/pipeline/v17/bem/LUSC_cnv.ro")
    commonCell <- intersect(rownames(LUSC_mut), intersect(rownames(LUSC_methyl), intersect(rownames(LUSC_gex), rownames(LUSC_cnv))))
    feature <- cbind(LUSC_mut[commonCell,], LUSC_methyl[commonCell,], LUSC_gex[commonCell,], LUSC_cnv[commonCell,])
  }
  if (feat=="v17_LUSC_mut.methyl.gex.speed"){
    load("data/pipeline/v17/bem/LUSC_mut.ro")
    load("data/pipeline/v17/bem/LUSC_methyl.ro")
    load("data/pipeline/v17/bem/LUSC_gex.ro")
    load("data/pipeline/v17/bem/LUSC_speed.ro")
    commonCell <- intersect(rownames(LUSC_mut), intersect(rownames(LUSC_methyl), intersect(rownames(LUSC_gex), rownames(LUSC_speed))))
    feature <- cbind(LUSC_mut[commonCell,], LUSC_methyl[commonCell,], LUSC_gex[commonCell,], LUSC_speed[commonCell,])
  }
  if (feat=="v17_LUSC_mut.methyl.cnv.speed"){
    load("data/pipeline/v17/bem/LUSC_mut.ro")
    load("data/pipeline/v17/bem/LUSC_methyl.ro")
    load("data/pipeline/v17/bem/LUSC_cnv.ro")
    load("data/pipeline/v17/bem/LUSC_speed.ro")
    commonCell <- intersect(rownames(LUSC_mut), intersect(rownames(LUSC_methyl), intersect(rownames(LUSC_cnv), rownames(LUSC_speed))))
    feature <- cbind(LUSC_mut[commonCell,], LUSC_methyl[commonCell,], LUSC_cnv[commonCell,], LUSC_speed[commonCell,])
  }
  if (feat=="v17_LUSC_mut.gex.cnv.speed"){
    load("data/pipeline/v17/bem/LUSC_mut.ro")
    load("data/pipeline/v17/bem/LUSC_gex.ro")
    load("data/pipeline/v17/bem/LUSC_cnv.ro")
    load("data/pipeline/v17/bem/LUSC_speed.ro")
    commonCell <- intersect(rownames(LUSC_mut), intersect(rownames(LUSC_gex), intersect(rownames(LUSC_cnv), rownames(LUSC_speed))))
    feature <- cbind(LUSC_mut[commonCell,], LUSC_gex[commonCell,], LUSC_cnv[commonCell,], LUSC_speed[commonCell,])
  }
  if (feat=="v17_LUSC_methyl.gex.cnv.speed"){
    load("data/pipeline/v17/bem/LUSC_methyl.ro")
    load("data/pipeline/v17/bem/LUSC_gex.ro")
    load("data/pipeline/v17/bem/LUSC_cnv.ro")
    load("data/pipeline/v17/bem/LUSC_speed.ro")
    commonCell <- intersect(rownames(LUSC_methyl), intersect(rownames(LUSC_gex), intersect(rownames(LUSC_cnv), rownames(LUSC_speed))))
    feature <- cbind(LUSC_methyl[commonCell,], LUSC_gex[commonCell,], LUSC_cnv[commonCell,], LUSC_speed[commonCell,])
  }
  if (feat=="v17_LUSC_mut.methyl.gex.cnv.speed"){
    load("data/pipeline/v17/bem/LUSC_mut.ro")
    load("data/pipeline/v17/bem/LUSC_methyl.ro")
    load("data/pipeline/v17/bem/LUSC_gex.ro")
    load("data/pipeline/v17/bem/LUSC_cnv.ro")
    load("data/pipeline/v17/bem/LUSC_speed.ro")
    commonCell <- intersect(rownames(LUSC_mut), intersect(rownames(LUSC_methyl), intersect(rownames(LUSC_gex), intersect(rownames(LUSC_cnv), rownames(LUSC_speed)))))
    feature <- cbind(LUSC_mut[commonCell,], LUSC_methyl[commonCell,], LUSC_gex[commonCell,], LUSC_cnv[commonCell,], LUSC_speed[commonCell,])
  }
  if (feat=="v17_OV_mut"){
    load("data/pipeline/v17/bem/OV_mut.ro")
    commonCell <- rownames(OV_mut)
    feature <- cbind(OV_mut[commonCell,])
  }
  if (feat=="v17_OV_methyl"){
    load("data/pipeline/v17/bem/OV_methyl.ro")
    commonCell <- rownames(OV_methyl)
    feature <- cbind(OV_methyl[commonCell,])
  }
  if (feat=="v17_OV_gex"){
    load("data/pipeline/v17/bem/OV_gex.ro")
    commonCell <- rownames(OV_gex)
    feature <- cbind(OV_gex[commonCell,])
  }
  if (feat=="v17_OV_cnv"){
    load("data/pipeline/v17/bem/OV_cnv.ro")
    commonCell <- rownames(OV_cnv)
    feature <- cbind(OV_cnv[commonCell,])
  }
  if (feat=="v17_OV_speed"){
    load("data/pipeline/v17/bem/OV_speed.ro")
    commonCell <- rownames(OV_speed)
    feature <- cbind(OV_speed[commonCell,])
  }
  if (feat=="v17_OV_mut.methyl"){
    load("data/pipeline/v17/bem/OV_mut.ro")
    load("data/pipeline/v17/bem/OV_methyl.ro")
    commonCell <- intersect(rownames(OV_mut), rownames(OV_methyl))
    feature <- cbind(OV_mut[commonCell,], OV_methyl[commonCell,])
  }
  if (feat=="v17_OV_mut.gex"){
    load("data/pipeline/v17/bem/OV_mut.ro")
    load("data/pipeline/v17/bem/OV_gex.ro")
    commonCell <- intersect(rownames(OV_mut), rownames(OV_gex))
    feature <- cbind(OV_mut[commonCell,], OV_gex[commonCell,])
  }
  if (feat=="v17_OV_mut.cnv"){
    load("data/pipeline/v17/bem/OV_mut.ro")
    load("data/pipeline/v17/bem/OV_cnv.ro")
    commonCell <- intersect(rownames(OV_mut), rownames(OV_cnv))
    feature <- cbind(OV_mut[commonCell,], OV_cnv[commonCell,])
  }
  if (feat=="v17_OV_mut.speed"){
    load("data/pipeline/v17/bem/OV_mut.ro")
    load("data/pipeline/v17/bem/OV_speed.ro")
    commonCell <- intersect(rownames(OV_mut), rownames(OV_speed))
    feature <- cbind(OV_mut[commonCell,], OV_speed[commonCell,])
  }
  if (feat=="v17_OV_methyl.gex"){
    load("data/pipeline/v17/bem/OV_methyl.ro")
    load("data/pipeline/v17/bem/OV_gex.ro")
    commonCell <- intersect(rownames(OV_methyl), rownames(OV_gex))
    feature <- cbind(OV_methyl[commonCell,], OV_gex[commonCell,])
  }
  if (feat=="v17_OV_methyl.cnv"){
    load("data/pipeline/v17/bem/OV_methyl.ro")
    load("data/pipeline/v17/bem/OV_cnv.ro")
    commonCell <- intersect(rownames(OV_methyl), rownames(OV_cnv))
    feature <- cbind(OV_methyl[commonCell,], OV_cnv[commonCell,])
  }
  if (feat=="v17_OV_methyl.speed"){
    load("data/pipeline/v17/bem/OV_methyl.ro")
    load("data/pipeline/v17/bem/OV_speed.ro")
    commonCell <- intersect(rownames(OV_methyl), rownames(OV_speed))
    feature <- cbind(OV_methyl[commonCell,], OV_speed[commonCell,])
  }
  if (feat=="v17_OV_gex.cnv"){
    load("data/pipeline/v17/bem/OV_gex.ro")
    load("data/pipeline/v17/bem/OV_cnv.ro")
    commonCell <- intersect(rownames(OV_gex), rownames(OV_cnv))
    feature <- cbind(OV_gex[commonCell,], OV_cnv[commonCell,])
  }
  if (feat=="v17_OV_gex.speed"){
    load("data/pipeline/v17/bem/OV_gex.ro")
    load("data/pipeline/v17/bem/OV_speed.ro")
    commonCell <- intersect(rownames(OV_gex), rownames(OV_speed))
    feature <- cbind(OV_gex[commonCell,], OV_speed[commonCell,])
  }
  if (feat=="v17_OV_cnv.speed"){
    load("data/pipeline/v17/bem/OV_cnv.ro")
    load("data/pipeline/v17/bem/OV_speed.ro")
    commonCell <- intersect(rownames(OV_cnv), rownames(OV_speed))
    feature <- cbind(OV_cnv[commonCell,], OV_speed[commonCell,])
  }
  if (feat=="v17_OV_mut.methyl.gex"){
    load("data/pipeline/v17/bem/OV_mut.ro")
    load("data/pipeline/v17/bem/OV_methyl.ro")
    load("data/pipeline/v17/bem/OV_gex.ro")
    commonCell <- intersect(rownames(OV_mut), intersect(rownames(OV_methyl), rownames(OV_gex)))
    feature <- cbind(OV_mut[commonCell,], OV_methyl[commonCell,], OV_gex[commonCell,])
  }
  if (feat=="v17_OV_mut.methyl.cnv"){
    load("data/pipeline/v17/bem/OV_mut.ro")
    load("data/pipeline/v17/bem/OV_methyl.ro")
    load("data/pipeline/v17/bem/OV_cnv.ro")
    commonCell <- intersect(rownames(OV_mut), intersect(rownames(OV_methyl), rownames(OV_cnv)))
    feature <- cbind(OV_mut[commonCell,], OV_methyl[commonCell,], OV_cnv[commonCell,])
  }
  if (feat=="v17_OV_mut.methyl.speed"){
    load("data/pipeline/v17/bem/OV_mut.ro")
    load("data/pipeline/v17/bem/OV_methyl.ro")
    load("data/pipeline/v17/bem/OV_speed.ro")
    commonCell <- intersect(rownames(OV_mut), intersect(rownames(OV_methyl), rownames(OV_speed)))
    feature <- cbind(OV_mut[commonCell,], OV_methyl[commonCell,], OV_speed[commonCell,])
  }
  if (feat=="v17_OV_mut.gex.cnv"){
    load("data/pipeline/v17/bem/OV_mut.ro")
    load("data/pipeline/v17/bem/OV_gex.ro")
    load("data/pipeline/v17/bem/OV_cnv.ro")
    commonCell <- intersect(rownames(OV_mut), intersect(rownames(OV_gex), rownames(OV_cnv)))
    feature <- cbind(OV_mut[commonCell,], OV_gex[commonCell,], OV_cnv[commonCell,])
  }
  if (feat=="v17_OV_mut.gex.speed"){
    load("data/pipeline/v17/bem/OV_mut.ro")
    load("data/pipeline/v17/bem/OV_gex.ro")
    load("data/pipeline/v17/bem/OV_speed.ro")
    commonCell <- intersect(rownames(OV_mut), intersect(rownames(OV_gex), rownames(OV_speed)))
    feature <- cbind(OV_mut[commonCell,], OV_gex[commonCell,], OV_speed[commonCell,])
  }
  if (feat=="v17_OV_mut.cnv.speed"){
    load("data/pipeline/v17/bem/OV_mut.ro")
    load("data/pipeline/v17/bem/OV_cnv.ro")
    load("data/pipeline/v17/bem/OV_speed.ro")
    commonCell <- intersect(rownames(OV_mut), intersect(rownames(OV_cnv), rownames(OV_speed)))
    feature <- cbind(OV_mut[commonCell,], OV_cnv[commonCell,], OV_speed[commonCell,])
  }
  if (feat=="v17_OV_methyl.gex.cnv"){
    load("data/pipeline/v17/bem/OV_methyl.ro")
    load("data/pipeline/v17/bem/OV_gex.ro")
    load("data/pipeline/v17/bem/OV_cnv.ro")
    commonCell <- intersect(rownames(OV_methyl), intersect(rownames(OV_gex), rownames(OV_cnv)))
    feature <- cbind(OV_methyl[commonCell,], OV_gex[commonCell,], OV_cnv[commonCell,])
  }
  if (feat=="v17_OV_methyl.gex.speed"){
    load("data/pipeline/v17/bem/OV_methyl.ro")
    load("data/pipeline/v17/bem/OV_gex.ro")
    load("data/pipeline/v17/bem/OV_speed.ro")
    commonCell <- intersect(rownames(OV_methyl), intersect(rownames(OV_gex), rownames(OV_speed)))
    feature <- cbind(OV_methyl[commonCell,], OV_gex[commonCell,], OV_speed[commonCell,])
  }
  if (feat=="v17_OV_methyl.cnv.speed"){
    load("data/pipeline/v17/bem/OV_methyl.ro")
    load("data/pipeline/v17/bem/OV_cnv.ro")
    load("data/pipeline/v17/bem/OV_speed.ro")
    commonCell <- intersect(rownames(OV_methyl), intersect(rownames(OV_cnv), rownames(OV_speed)))
    feature <- cbind(OV_methyl[commonCell,], OV_cnv[commonCell,], OV_speed[commonCell,])
  }
  if (feat=="v17_OV_gex.cnv.speed"){
    load("data/pipeline/v17/bem/OV_gex.ro")
    load("data/pipeline/v17/bem/OV_cnv.ro")
    load("data/pipeline/v17/bem/OV_speed.ro")
    commonCell <- intersect(rownames(OV_gex), intersect(rownames(OV_cnv), rownames(OV_speed)))
    feature <- cbind(OV_gex[commonCell,], OV_cnv[commonCell,], OV_speed[commonCell,])
  }
  if (feat=="v17_OV_mut.methyl.gex.cnv"){
    load("data/pipeline/v17/bem/OV_mut.ro")
    load("data/pipeline/v17/bem/OV_methyl.ro")
    load("data/pipeline/v17/bem/OV_gex.ro")
    load("data/pipeline/v17/bem/OV_cnv.ro")
    commonCell <- intersect(rownames(OV_mut), intersect(rownames(OV_methyl), intersect(rownames(OV_gex), rownames(OV_cnv))))
    feature <- cbind(OV_mut[commonCell,], OV_methyl[commonCell,], OV_gex[commonCell,], OV_cnv[commonCell,])
  }
  if (feat=="v17_OV_mut.methyl.gex.speed"){
    load("data/pipeline/v17/bem/OV_mut.ro")
    load("data/pipeline/v17/bem/OV_methyl.ro")
    load("data/pipeline/v17/bem/OV_gex.ro")
    load("data/pipeline/v17/bem/OV_speed.ro")
    commonCell <- intersect(rownames(OV_mut), intersect(rownames(OV_methyl), intersect(rownames(OV_gex), rownames(OV_speed))))
    feature <- cbind(OV_mut[commonCell,], OV_methyl[commonCell,], OV_gex[commonCell,], OV_speed[commonCell,])
  }
  if (feat=="v17_OV_mut.methyl.cnv.speed"){
    load("data/pipeline/v17/bem/OV_mut.ro")
    load("data/pipeline/v17/bem/OV_methyl.ro")
    load("data/pipeline/v17/bem/OV_cnv.ro")
    load("data/pipeline/v17/bem/OV_speed.ro")
    commonCell <- intersect(rownames(OV_mut), intersect(rownames(OV_methyl), intersect(rownames(OV_cnv), rownames(OV_speed))))
    feature <- cbind(OV_mut[commonCell,], OV_methyl[commonCell,], OV_cnv[commonCell,], OV_speed[commonCell,])
  }
  if (feat=="v17_OV_mut.gex.cnv.speed"){
    load("data/pipeline/v17/bem/OV_mut.ro")
    load("data/pipeline/v17/bem/OV_gex.ro")
    load("data/pipeline/v17/bem/OV_cnv.ro")
    load("data/pipeline/v17/bem/OV_speed.ro")
    commonCell <- intersect(rownames(OV_mut), intersect(rownames(OV_gex), intersect(rownames(OV_cnv), rownames(OV_speed))))
    feature <- cbind(OV_mut[commonCell,], OV_gex[commonCell,], OV_cnv[commonCell,], OV_speed[commonCell,])
  }
  if (feat=="v17_OV_methyl.gex.cnv.speed"){
    load("data/pipeline/v17/bem/OV_methyl.ro")
    load("data/pipeline/v17/bem/OV_gex.ro")
    load("data/pipeline/v17/bem/OV_cnv.ro")
    load("data/pipeline/v17/bem/OV_speed.ro")
    commonCell <- intersect(rownames(OV_methyl), intersect(rownames(OV_gex), intersect(rownames(OV_cnv), rownames(OV_speed))))
    feature <- cbind(OV_methyl[commonCell,], OV_gex[commonCell,], OV_cnv[commonCell,], OV_speed[commonCell,])
  }
  if (feat=="v17_OV_mut.methyl.gex.cnv.speed"){
    load("data/pipeline/v17/bem/OV_mut.ro")
    load("data/pipeline/v17/bem/OV_methyl.ro")
    load("data/pipeline/v17/bem/OV_gex.ro")
    load("data/pipeline/v17/bem/OV_cnv.ro")
    load("data/pipeline/v17/bem/OV_speed.ro")
    commonCell <- intersect(rownames(OV_mut), intersect(rownames(OV_methyl), intersect(rownames(OV_gex), intersect(rownames(OV_cnv), rownames(OV_speed)))))
    feature <- cbind(OV_mut[commonCell,], OV_methyl[commonCell,], OV_gex[commonCell,], OV_cnv[commonCell,], OV_speed[commonCell,])
  }
  if (feat=="v17_PAAD_mut"){
    load("data/pipeline/v17/bem/PAAD_mut.ro")
    commonCell <- rownames(PAAD_mut)
    feature <- cbind(PAAD_mut[commonCell,])
  }
  if (feat=="v17_PAAD_methyl"){
    load("data/pipeline/v17/bem/PAAD_methyl.ro")
    commonCell <- rownames(PAAD_methyl)
    feature <- cbind(PAAD_methyl[commonCell,])
  }
  if (feat=="v17_PAAD_gex"){
    load("data/pipeline/v17/bem/PAAD_gex.ro")
    commonCell <- rownames(PAAD_gex)
    feature <- cbind(PAAD_gex[commonCell,])
  }
  if (feat=="v17_PAAD_cnv"){
    load("data/pipeline/v17/bem/PAAD_cnv.ro")
    commonCell <- rownames(PAAD_cnv)
    feature <- cbind(PAAD_cnv[commonCell,])
  }
  if (feat=="v17_PAAD_speed"){
    load("data/pipeline/v17/bem/PAAD_speed.ro")
    commonCell <- rownames(PAAD_speed)
    feature <- cbind(PAAD_speed[commonCell,])
  }
  if (feat=="v17_PAAD_mut.methyl"){
    load("data/pipeline/v17/bem/PAAD_mut.ro")
    load("data/pipeline/v17/bem/PAAD_methyl.ro")
    commonCell <- intersect(rownames(PAAD_mut), rownames(PAAD_methyl))
    feature <- cbind(PAAD_mut[commonCell,], PAAD_methyl[commonCell,])
  }
  if (feat=="v17_PAAD_mut.gex"){
    load("data/pipeline/v17/bem/PAAD_mut.ro")
    load("data/pipeline/v17/bem/PAAD_gex.ro")
    commonCell <- intersect(rownames(PAAD_mut), rownames(PAAD_gex))
    feature <- cbind(PAAD_mut[commonCell,], PAAD_gex[commonCell,])
  }
  if (feat=="v17_PAAD_mut.cnv"){
    load("data/pipeline/v17/bem/PAAD_mut.ro")
    load("data/pipeline/v17/bem/PAAD_cnv.ro")
    commonCell <- intersect(rownames(PAAD_mut), rownames(PAAD_cnv))
    feature <- cbind(PAAD_mut[commonCell,], PAAD_cnv[commonCell,])
  }
  if (feat=="v17_PAAD_mut.speed"){
    load("data/pipeline/v17/bem/PAAD_mut.ro")
    load("data/pipeline/v17/bem/PAAD_speed.ro")
    commonCell <- intersect(rownames(PAAD_mut), rownames(PAAD_speed))
    feature <- cbind(PAAD_mut[commonCell,], PAAD_speed[commonCell,])
  }
  if (feat=="v17_PAAD_methyl.gex"){
    load("data/pipeline/v17/bem/PAAD_methyl.ro")
    load("data/pipeline/v17/bem/PAAD_gex.ro")
    commonCell <- intersect(rownames(PAAD_methyl), rownames(PAAD_gex))
    feature <- cbind(PAAD_methyl[commonCell,], PAAD_gex[commonCell,])
  }
  if (feat=="v17_PAAD_methyl.cnv"){
    load("data/pipeline/v17/bem/PAAD_methyl.ro")
    load("data/pipeline/v17/bem/PAAD_cnv.ro")
    commonCell <- intersect(rownames(PAAD_methyl), rownames(PAAD_cnv))
    feature <- cbind(PAAD_methyl[commonCell,], PAAD_cnv[commonCell,])
  }
  if (feat=="v17_PAAD_methyl.speed"){
    load("data/pipeline/v17/bem/PAAD_methyl.ro")
    load("data/pipeline/v17/bem/PAAD_speed.ro")
    commonCell <- intersect(rownames(PAAD_methyl), rownames(PAAD_speed))
    feature <- cbind(PAAD_methyl[commonCell,], PAAD_speed[commonCell,])
  }
  if (feat=="v17_PAAD_gex.cnv"){
    load("data/pipeline/v17/bem/PAAD_gex.ro")
    load("data/pipeline/v17/bem/PAAD_cnv.ro")
    commonCell <- intersect(rownames(PAAD_gex), rownames(PAAD_cnv))
    feature <- cbind(PAAD_gex[commonCell,], PAAD_cnv[commonCell,])
  }
  if (feat=="v17_PAAD_gex.speed"){
    load("data/pipeline/v17/bem/PAAD_gex.ro")
    load("data/pipeline/v17/bem/PAAD_speed.ro")
    commonCell <- intersect(rownames(PAAD_gex), rownames(PAAD_speed))
    feature <- cbind(PAAD_gex[commonCell,], PAAD_speed[commonCell,])
  }
  if (feat=="v17_PAAD_cnv.speed"){
    load("data/pipeline/v17/bem/PAAD_cnv.ro")
    load("data/pipeline/v17/bem/PAAD_speed.ro")
    commonCell <- intersect(rownames(PAAD_cnv), rownames(PAAD_speed))
    feature <- cbind(PAAD_cnv[commonCell,], PAAD_speed[commonCell,])
  }
  if (feat=="v17_PAAD_mut.methyl.gex"){
    load("data/pipeline/v17/bem/PAAD_mut.ro")
    load("data/pipeline/v17/bem/PAAD_methyl.ro")
    load("data/pipeline/v17/bem/PAAD_gex.ro")
    commonCell <- intersect(rownames(PAAD_mut), intersect(rownames(PAAD_methyl), rownames(PAAD_gex)))
    feature <- cbind(PAAD_mut[commonCell,], PAAD_methyl[commonCell,], PAAD_gex[commonCell,])
  }
  if (feat=="v17_PAAD_mut.methyl.cnv"){
    load("data/pipeline/v17/bem/PAAD_mut.ro")
    load("data/pipeline/v17/bem/PAAD_methyl.ro")
    load("data/pipeline/v17/bem/PAAD_cnv.ro")
    commonCell <- intersect(rownames(PAAD_mut), intersect(rownames(PAAD_methyl), rownames(PAAD_cnv)))
    feature <- cbind(PAAD_mut[commonCell,], PAAD_methyl[commonCell,], PAAD_cnv[commonCell,])
  }
  if (feat=="v17_PAAD_mut.methyl.speed"){
    load("data/pipeline/v17/bem/PAAD_mut.ro")
    load("data/pipeline/v17/bem/PAAD_methyl.ro")
    load("data/pipeline/v17/bem/PAAD_speed.ro")
    commonCell <- intersect(rownames(PAAD_mut), intersect(rownames(PAAD_methyl), rownames(PAAD_speed)))
    feature <- cbind(PAAD_mut[commonCell,], PAAD_methyl[commonCell,], PAAD_speed[commonCell,])
  }
  if (feat=="v17_PAAD_mut.gex.cnv"){
    load("data/pipeline/v17/bem/PAAD_mut.ro")
    load("data/pipeline/v17/bem/PAAD_gex.ro")
    load("data/pipeline/v17/bem/PAAD_cnv.ro")
    commonCell <- intersect(rownames(PAAD_mut), intersect(rownames(PAAD_gex), rownames(PAAD_cnv)))
    feature <- cbind(PAAD_mut[commonCell,], PAAD_gex[commonCell,], PAAD_cnv[commonCell,])
  }
  if (feat=="v17_PAAD_mut.gex.speed"){
    load("data/pipeline/v17/bem/PAAD_mut.ro")
    load("data/pipeline/v17/bem/PAAD_gex.ro")
    load("data/pipeline/v17/bem/PAAD_speed.ro")
    commonCell <- intersect(rownames(PAAD_mut), intersect(rownames(PAAD_gex), rownames(PAAD_speed)))
    feature <- cbind(PAAD_mut[commonCell,], PAAD_gex[commonCell,], PAAD_speed[commonCell,])
  }
  if (feat=="v17_PAAD_mut.cnv.speed"){
    load("data/pipeline/v17/bem/PAAD_mut.ro")
    load("data/pipeline/v17/bem/PAAD_cnv.ro")
    load("data/pipeline/v17/bem/PAAD_speed.ro")
    commonCell <- intersect(rownames(PAAD_mut), intersect(rownames(PAAD_cnv), rownames(PAAD_speed)))
    feature <- cbind(PAAD_mut[commonCell,], PAAD_cnv[commonCell,], PAAD_speed[commonCell,])
  }
  if (feat=="v17_PAAD_methyl.gex.cnv"){
    load("data/pipeline/v17/bem/PAAD_methyl.ro")
    load("data/pipeline/v17/bem/PAAD_gex.ro")
    load("data/pipeline/v17/bem/PAAD_cnv.ro")
    commonCell <- intersect(rownames(PAAD_methyl), intersect(rownames(PAAD_gex), rownames(PAAD_cnv)))
    feature <- cbind(PAAD_methyl[commonCell,], PAAD_gex[commonCell,], PAAD_cnv[commonCell,])
  }
  if (feat=="v17_PAAD_methyl.gex.speed"){
    load("data/pipeline/v17/bem/PAAD_methyl.ro")
    load("data/pipeline/v17/bem/PAAD_gex.ro")
    load("data/pipeline/v17/bem/PAAD_speed.ro")
    commonCell <- intersect(rownames(PAAD_methyl), intersect(rownames(PAAD_gex), rownames(PAAD_speed)))
    feature <- cbind(PAAD_methyl[commonCell,], PAAD_gex[commonCell,], PAAD_speed[commonCell,])
  }
  if (feat=="v17_PAAD_methyl.cnv.speed"){
    load("data/pipeline/v17/bem/PAAD_methyl.ro")
    load("data/pipeline/v17/bem/PAAD_cnv.ro")
    load("data/pipeline/v17/bem/PAAD_speed.ro")
    commonCell <- intersect(rownames(PAAD_methyl), intersect(rownames(PAAD_cnv), rownames(PAAD_speed)))
    feature <- cbind(PAAD_methyl[commonCell,], PAAD_cnv[commonCell,], PAAD_speed[commonCell,])
  }
  if (feat=="v17_PAAD_gex.cnv.speed"){
    load("data/pipeline/v17/bem/PAAD_gex.ro")
    load("data/pipeline/v17/bem/PAAD_cnv.ro")
    load("data/pipeline/v17/bem/PAAD_speed.ro")
    commonCell <- intersect(rownames(PAAD_gex), intersect(rownames(PAAD_cnv), rownames(PAAD_speed)))
    feature <- cbind(PAAD_gex[commonCell,], PAAD_cnv[commonCell,], PAAD_speed[commonCell,])
  }
  if (feat=="v17_PAAD_mut.methyl.gex.cnv"){
    load("data/pipeline/v17/bem/PAAD_mut.ro")
    load("data/pipeline/v17/bem/PAAD_methyl.ro")
    load("data/pipeline/v17/bem/PAAD_gex.ro")
    load("data/pipeline/v17/bem/PAAD_cnv.ro")
    commonCell <- intersect(rownames(PAAD_mut), intersect(rownames(PAAD_methyl), intersect(rownames(PAAD_gex), rownames(PAAD_cnv))))
    feature <- cbind(PAAD_mut[commonCell,], PAAD_methyl[commonCell,], PAAD_gex[commonCell,], PAAD_cnv[commonCell,])
  }
  if (feat=="v17_PAAD_mut.methyl.gex.speed"){
    load("data/pipeline/v17/bem/PAAD_mut.ro")
    load("data/pipeline/v17/bem/PAAD_methyl.ro")
    load("data/pipeline/v17/bem/PAAD_gex.ro")
    load("data/pipeline/v17/bem/PAAD_speed.ro")
    commonCell <- intersect(rownames(PAAD_mut), intersect(rownames(PAAD_methyl), intersect(rownames(PAAD_gex), rownames(PAAD_speed))))
    feature <- cbind(PAAD_mut[commonCell,], PAAD_methyl[commonCell,], PAAD_gex[commonCell,], PAAD_speed[commonCell,])
  }
  if (feat=="v17_PAAD_mut.methyl.cnv.speed"){
    load("data/pipeline/v17/bem/PAAD_mut.ro")
    load("data/pipeline/v17/bem/PAAD_methyl.ro")
    load("data/pipeline/v17/bem/PAAD_cnv.ro")
    load("data/pipeline/v17/bem/PAAD_speed.ro")
    commonCell <- intersect(rownames(PAAD_mut), intersect(rownames(PAAD_methyl), intersect(rownames(PAAD_cnv), rownames(PAAD_speed))))
    feature <- cbind(PAAD_mut[commonCell,], PAAD_methyl[commonCell,], PAAD_cnv[commonCell,], PAAD_speed[commonCell,])
  }
  if (feat=="v17_PAAD_mut.gex.cnv.speed"){
    load("data/pipeline/v17/bem/PAAD_mut.ro")
    load("data/pipeline/v17/bem/PAAD_gex.ro")
    load("data/pipeline/v17/bem/PAAD_cnv.ro")
    load("data/pipeline/v17/bem/PAAD_speed.ro")
    commonCell <- intersect(rownames(PAAD_mut), intersect(rownames(PAAD_gex), intersect(rownames(PAAD_cnv), rownames(PAAD_speed))))
    feature <- cbind(PAAD_mut[commonCell,], PAAD_gex[commonCell,], PAAD_cnv[commonCell,], PAAD_speed[commonCell,])
  }
  if (feat=="v17_PAAD_methyl.gex.cnv.speed"){
    load("data/pipeline/v17/bem/PAAD_methyl.ro")
    load("data/pipeline/v17/bem/PAAD_gex.ro")
    load("data/pipeline/v17/bem/PAAD_cnv.ro")
    load("data/pipeline/v17/bem/PAAD_speed.ro")
    commonCell <- intersect(rownames(PAAD_methyl), intersect(rownames(PAAD_gex), intersect(rownames(PAAD_cnv), rownames(PAAD_speed))))
    feature <- cbind(PAAD_methyl[commonCell,], PAAD_gex[commonCell,], PAAD_cnv[commonCell,], PAAD_speed[commonCell,])
  }
  if (feat=="v17_PAAD_mut.methyl.gex.cnv.speed"){
    load("data/pipeline/v17/bem/PAAD_mut.ro")
    load("data/pipeline/v17/bem/PAAD_methyl.ro")
    load("data/pipeline/v17/bem/PAAD_gex.ro")
    load("data/pipeline/v17/bem/PAAD_cnv.ro")
    load("data/pipeline/v17/bem/PAAD_speed.ro")
    commonCell <- intersect(rownames(PAAD_mut), intersect(rownames(PAAD_methyl), intersect(rownames(PAAD_gex), intersect(rownames(PAAD_cnv), rownames(PAAD_speed)))))
    feature <- cbind(PAAD_mut[commonCell,], PAAD_methyl[commonCell,], PAAD_gex[commonCell,], PAAD_cnv[commonCell,], PAAD_speed[commonCell,])
  }
  if (feat=="v17_PRAD_mut"){
    load("data/pipeline/v17/bem/PRAD_mut.ro")
    commonCell <- rownames(PRAD_mut)
    feature <- cbind(PRAD_mut[commonCell,])
  }
  if (feat=="v17_PRAD_methyl"){
    load("data/pipeline/v17/bem/PRAD_methyl.ro")
    commonCell <- rownames(PRAD_methyl)
    feature <- cbind(PRAD_methyl[commonCell,])
  }
  if (feat=="v17_PRAD_gex"){
    load("data/pipeline/v17/bem/PRAD_gex.ro")
    commonCell <- rownames(PRAD_gex)
    feature <- cbind(PRAD_gex[commonCell,])
  }
  if (feat=="v17_PRAD_cnv"){
    load("data/pipeline/v17/bem/PRAD_cnv.ro")
    commonCell <- rownames(PRAD_cnv)
    feature <- cbind(PRAD_cnv[commonCell,])
  }
  if (feat=="v17_PRAD_speed"){
    load("data/pipeline/v17/bem/PRAD_speed.ro")
    commonCell <- rownames(PRAD_speed)
    feature <- cbind(PRAD_speed[commonCell,])
  }
  if (feat=="v17_PRAD_mut.methyl"){
    load("data/pipeline/v17/bem/PRAD_mut.ro")
    load("data/pipeline/v17/bem/PRAD_methyl.ro")
    commonCell <- intersect(rownames(PRAD_mut), rownames(PRAD_methyl))
    feature <- cbind(PRAD_mut[commonCell,], PRAD_methyl[commonCell,])
  }
  if (feat=="v17_PRAD_mut.gex"){
    load("data/pipeline/v17/bem/PRAD_mut.ro")
    load("data/pipeline/v17/bem/PRAD_gex.ro")
    commonCell <- intersect(rownames(PRAD_mut), rownames(PRAD_gex))
    feature <- cbind(PRAD_mut[commonCell,], PRAD_gex[commonCell,])
  }
  if (feat=="v17_PRAD_mut.cnv"){
    load("data/pipeline/v17/bem/PRAD_mut.ro")
    load("data/pipeline/v17/bem/PRAD_cnv.ro")
    commonCell <- intersect(rownames(PRAD_mut), rownames(PRAD_cnv))
    feature <- cbind(PRAD_mut[commonCell,], PRAD_cnv[commonCell,])
  }
  if (feat=="v17_PRAD_mut.speed"){
    load("data/pipeline/v17/bem/PRAD_mut.ro")
    load("data/pipeline/v17/bem/PRAD_speed.ro")
    commonCell <- intersect(rownames(PRAD_mut), rownames(PRAD_speed))
    feature <- cbind(PRAD_mut[commonCell,], PRAD_speed[commonCell,])
  }
  if (feat=="v17_PRAD_methyl.gex"){
    load("data/pipeline/v17/bem/PRAD_methyl.ro")
    load("data/pipeline/v17/bem/PRAD_gex.ro")
    commonCell <- intersect(rownames(PRAD_methyl), rownames(PRAD_gex))
    feature <- cbind(PRAD_methyl[commonCell,], PRAD_gex[commonCell,])
  }
  if (feat=="v17_PRAD_methyl.cnv"){
    load("data/pipeline/v17/bem/PRAD_methyl.ro")
    load("data/pipeline/v17/bem/PRAD_cnv.ro")
    commonCell <- intersect(rownames(PRAD_methyl), rownames(PRAD_cnv))
    feature <- cbind(PRAD_methyl[commonCell,], PRAD_cnv[commonCell,])
  }
  if (feat=="v17_PRAD_methyl.speed"){
    load("data/pipeline/v17/bem/PRAD_methyl.ro")
    load("data/pipeline/v17/bem/PRAD_speed.ro")
    commonCell <- intersect(rownames(PRAD_methyl), rownames(PRAD_speed))
    feature <- cbind(PRAD_methyl[commonCell,], PRAD_speed[commonCell,])
  }
  if (feat=="v17_PRAD_gex.cnv"){
    load("data/pipeline/v17/bem/PRAD_gex.ro")
    load("data/pipeline/v17/bem/PRAD_cnv.ro")
    commonCell <- intersect(rownames(PRAD_gex), rownames(PRAD_cnv))
    feature <- cbind(PRAD_gex[commonCell,], PRAD_cnv[commonCell,])
  }
  if (feat=="v17_PRAD_gex.speed"){
    load("data/pipeline/v17/bem/PRAD_gex.ro")
    load("data/pipeline/v17/bem/PRAD_speed.ro")
    commonCell <- intersect(rownames(PRAD_gex), rownames(PRAD_speed))
    feature <- cbind(PRAD_gex[commonCell,], PRAD_speed[commonCell,])
  }
  if (feat=="v17_PRAD_cnv.speed"){
    load("data/pipeline/v17/bem/PRAD_cnv.ro")
    load("data/pipeline/v17/bem/PRAD_speed.ro")
    commonCell <- intersect(rownames(PRAD_cnv), rownames(PRAD_speed))
    feature <- cbind(PRAD_cnv[commonCell,], PRAD_speed[commonCell,])
  }
  if (feat=="v17_PRAD_mut.methyl.gex"){
    load("data/pipeline/v17/bem/PRAD_mut.ro")
    load("data/pipeline/v17/bem/PRAD_methyl.ro")
    load("data/pipeline/v17/bem/PRAD_gex.ro")
    commonCell <- intersect(rownames(PRAD_mut), intersect(rownames(PRAD_methyl), rownames(PRAD_gex)))
    feature <- cbind(PRAD_mut[commonCell,], PRAD_methyl[commonCell,], PRAD_gex[commonCell,])
  }
  if (feat=="v17_PRAD_mut.methyl.cnv"){
    load("data/pipeline/v17/bem/PRAD_mut.ro")
    load("data/pipeline/v17/bem/PRAD_methyl.ro")
    load("data/pipeline/v17/bem/PRAD_cnv.ro")
    commonCell <- intersect(rownames(PRAD_mut), intersect(rownames(PRAD_methyl), rownames(PRAD_cnv)))
    feature <- cbind(PRAD_mut[commonCell,], PRAD_methyl[commonCell,], PRAD_cnv[commonCell,])
  }
  if (feat=="v17_PRAD_mut.methyl.speed"){
    load("data/pipeline/v17/bem/PRAD_mut.ro")
    load("data/pipeline/v17/bem/PRAD_methyl.ro")
    load("data/pipeline/v17/bem/PRAD_speed.ro")
    commonCell <- intersect(rownames(PRAD_mut), intersect(rownames(PRAD_methyl), rownames(PRAD_speed)))
    feature <- cbind(PRAD_mut[commonCell,], PRAD_methyl[commonCell,], PRAD_speed[commonCell,])
  }
  if (feat=="v17_PRAD_mut.gex.cnv"){
    load("data/pipeline/v17/bem/PRAD_mut.ro")
    load("data/pipeline/v17/bem/PRAD_gex.ro")
    load("data/pipeline/v17/bem/PRAD_cnv.ro")
    commonCell <- intersect(rownames(PRAD_mut), intersect(rownames(PRAD_gex), rownames(PRAD_cnv)))
    feature <- cbind(PRAD_mut[commonCell,], PRAD_gex[commonCell,], PRAD_cnv[commonCell,])
  }
  if (feat=="v17_PRAD_mut.gex.speed"){
    load("data/pipeline/v17/bem/PRAD_mut.ro")
    load("data/pipeline/v17/bem/PRAD_gex.ro")
    load("data/pipeline/v17/bem/PRAD_speed.ro")
    commonCell <- intersect(rownames(PRAD_mut), intersect(rownames(PRAD_gex), rownames(PRAD_speed)))
    feature <- cbind(PRAD_mut[commonCell,], PRAD_gex[commonCell,], PRAD_speed[commonCell,])
  }
  if (feat=="v17_PRAD_mut.cnv.speed"){
    load("data/pipeline/v17/bem/PRAD_mut.ro")
    load("data/pipeline/v17/bem/PRAD_cnv.ro")
    load("data/pipeline/v17/bem/PRAD_speed.ro")
    commonCell <- intersect(rownames(PRAD_mut), intersect(rownames(PRAD_cnv), rownames(PRAD_speed)))
    feature <- cbind(PRAD_mut[commonCell,], PRAD_cnv[commonCell,], PRAD_speed[commonCell,])
  }
  if (feat=="v17_PRAD_methyl.gex.cnv"){
    load("data/pipeline/v17/bem/PRAD_methyl.ro")
    load("data/pipeline/v17/bem/PRAD_gex.ro")
    load("data/pipeline/v17/bem/PRAD_cnv.ro")
    commonCell <- intersect(rownames(PRAD_methyl), intersect(rownames(PRAD_gex), rownames(PRAD_cnv)))
    feature <- cbind(PRAD_methyl[commonCell,], PRAD_gex[commonCell,], PRAD_cnv[commonCell,])
  }
  if (feat=="v17_PRAD_methyl.gex.speed"){
    load("data/pipeline/v17/bem/PRAD_methyl.ro")
    load("data/pipeline/v17/bem/PRAD_gex.ro")
    load("data/pipeline/v17/bem/PRAD_speed.ro")
    commonCell <- intersect(rownames(PRAD_methyl), intersect(rownames(PRAD_gex), rownames(PRAD_speed)))
    feature <- cbind(PRAD_methyl[commonCell,], PRAD_gex[commonCell,], PRAD_speed[commonCell,])
  }
  if (feat=="v17_PRAD_methyl.cnv.speed"){
    load("data/pipeline/v17/bem/PRAD_methyl.ro")
    load("data/pipeline/v17/bem/PRAD_cnv.ro")
    load("data/pipeline/v17/bem/PRAD_speed.ro")
    commonCell <- intersect(rownames(PRAD_methyl), intersect(rownames(PRAD_cnv), rownames(PRAD_speed)))
    feature <- cbind(PRAD_methyl[commonCell,], PRAD_cnv[commonCell,], PRAD_speed[commonCell,])
  }
  if (feat=="v17_PRAD_gex.cnv.speed"){
    load("data/pipeline/v17/bem/PRAD_gex.ro")
    load("data/pipeline/v17/bem/PRAD_cnv.ro")
    load("data/pipeline/v17/bem/PRAD_speed.ro")
    commonCell <- intersect(rownames(PRAD_gex), intersect(rownames(PRAD_cnv), rownames(PRAD_speed)))
    feature <- cbind(PRAD_gex[commonCell,], PRAD_cnv[commonCell,], PRAD_speed[commonCell,])
  }
  if (feat=="v17_PRAD_mut.methyl.gex.cnv"){
    load("data/pipeline/v17/bem/PRAD_mut.ro")
    load("data/pipeline/v17/bem/PRAD_methyl.ro")
    load("data/pipeline/v17/bem/PRAD_gex.ro")
    load("data/pipeline/v17/bem/PRAD_cnv.ro")
    commonCell <- intersect(rownames(PRAD_mut), intersect(rownames(PRAD_methyl), intersect(rownames(PRAD_gex), rownames(PRAD_cnv))))
    feature <- cbind(PRAD_mut[commonCell,], PRAD_methyl[commonCell,], PRAD_gex[commonCell,], PRAD_cnv[commonCell,])
  }
  if (feat=="v17_PRAD_mut.methyl.gex.speed"){
    load("data/pipeline/v17/bem/PRAD_mut.ro")
    load("data/pipeline/v17/bem/PRAD_methyl.ro")
    load("data/pipeline/v17/bem/PRAD_gex.ro")
    load("data/pipeline/v17/bem/PRAD_speed.ro")
    commonCell <- intersect(rownames(PRAD_mut), intersect(rownames(PRAD_methyl), intersect(rownames(PRAD_gex), rownames(PRAD_speed))))
    feature <- cbind(PRAD_mut[commonCell,], PRAD_methyl[commonCell,], PRAD_gex[commonCell,], PRAD_speed[commonCell,])
  }
  if (feat=="v17_PRAD_mut.methyl.cnv.speed"){
    load("data/pipeline/v17/bem/PRAD_mut.ro")
    load("data/pipeline/v17/bem/PRAD_methyl.ro")
    load("data/pipeline/v17/bem/PRAD_cnv.ro")
    load("data/pipeline/v17/bem/PRAD_speed.ro")
    commonCell <- intersect(rownames(PRAD_mut), intersect(rownames(PRAD_methyl), intersect(rownames(PRAD_cnv), rownames(PRAD_speed))))
    feature <- cbind(PRAD_mut[commonCell,], PRAD_methyl[commonCell,], PRAD_cnv[commonCell,], PRAD_speed[commonCell,])
  }
  if (feat=="v17_PRAD_mut.gex.cnv.speed"){
    load("data/pipeline/v17/bem/PRAD_mut.ro")
    load("data/pipeline/v17/bem/PRAD_gex.ro")
    load("data/pipeline/v17/bem/PRAD_cnv.ro")
    load("data/pipeline/v17/bem/PRAD_speed.ro")
    commonCell <- intersect(rownames(PRAD_mut), intersect(rownames(PRAD_gex), intersect(rownames(PRAD_cnv), rownames(PRAD_speed))))
    feature <- cbind(PRAD_mut[commonCell,], PRAD_gex[commonCell,], PRAD_cnv[commonCell,], PRAD_speed[commonCell,])
  }
  if (feat=="v17_PRAD_methyl.gex.cnv.speed"){
    load("data/pipeline/v17/bem/PRAD_methyl.ro")
    load("data/pipeline/v17/bem/PRAD_gex.ro")
    load("data/pipeline/v17/bem/PRAD_cnv.ro")
    load("data/pipeline/v17/bem/PRAD_speed.ro")
    commonCell <- intersect(rownames(PRAD_methyl), intersect(rownames(PRAD_gex), intersect(rownames(PRAD_cnv), rownames(PRAD_speed))))
    feature <- cbind(PRAD_methyl[commonCell,], PRAD_gex[commonCell,], PRAD_cnv[commonCell,], PRAD_speed[commonCell,])
  }
  if (feat=="v17_PRAD_mut.methyl.gex.cnv.speed"){
    load("data/pipeline/v17/bem/PRAD_mut.ro")
    load("data/pipeline/v17/bem/PRAD_methyl.ro")
    load("data/pipeline/v17/bem/PRAD_gex.ro")
    load("data/pipeline/v17/bem/PRAD_cnv.ro")
    load("data/pipeline/v17/bem/PRAD_speed.ro")
    commonCell <- intersect(rownames(PRAD_mut), intersect(rownames(PRAD_methyl), intersect(rownames(PRAD_gex), intersect(rownames(PRAD_cnv), rownames(PRAD_speed)))))
    feature <- cbind(PRAD_mut[commonCell,], PRAD_methyl[commonCell,], PRAD_gex[commonCell,], PRAD_cnv[commonCell,], PRAD_speed[commonCell,])
  }
  if (feat=="v17_SKCM_mut"){
    load("data/pipeline/v17/bem/SKCM_mut.ro")
    commonCell <- rownames(SKCM_mut)
    feature <- cbind(SKCM_mut[commonCell,])
  }
  if (feat=="v17_SKCM_methyl"){
    load("data/pipeline/v17/bem/SKCM_methyl.ro")
    commonCell <- rownames(SKCM_methyl)
    feature <- cbind(SKCM_methyl[commonCell,])
  }
  if (feat=="v17_SKCM_gex"){
    load("data/pipeline/v17/bem/SKCM_gex.ro")
    commonCell <- rownames(SKCM_gex)
    feature <- cbind(SKCM_gex[commonCell,])
  }
  if (feat=="v17_SKCM_cnv"){
    load("data/pipeline/v17/bem/SKCM_cnv.ro")
    commonCell <- rownames(SKCM_cnv)
    feature <- cbind(SKCM_cnv[commonCell,])
  }
  if (feat=="v17_SKCM_speed"){
    load("data/pipeline/v17/bem/SKCM_speed.ro")
    commonCell <- rownames(SKCM_speed)
    feature <- cbind(SKCM_speed[commonCell,])
  }
  if (feat=="v17_SKCM_mut.methyl"){
    load("data/pipeline/v17/bem/SKCM_mut.ro")
    load("data/pipeline/v17/bem/SKCM_methyl.ro")
    commonCell <- intersect(rownames(SKCM_mut), rownames(SKCM_methyl))
    feature <- cbind(SKCM_mut[commonCell,], SKCM_methyl[commonCell,])
  }
  if (feat=="v17_SKCM_mut.gex"){
    load("data/pipeline/v17/bem/SKCM_mut.ro")
    load("data/pipeline/v17/bem/SKCM_gex.ro")
    commonCell <- intersect(rownames(SKCM_mut), rownames(SKCM_gex))
    feature <- cbind(SKCM_mut[commonCell,], SKCM_gex[commonCell,])
  }
  if (feat=="v17_SKCM_mut.cnv"){
    load("data/pipeline/v17/bem/SKCM_mut.ro")
    load("data/pipeline/v17/bem/SKCM_cnv.ro")
    commonCell <- intersect(rownames(SKCM_mut), rownames(SKCM_cnv))
    feature <- cbind(SKCM_mut[commonCell,], SKCM_cnv[commonCell,])
  }
  if (feat=="v17_SKCM_mut.speed"){
    load("data/pipeline/v17/bem/SKCM_mut.ro")
    load("data/pipeline/v17/bem/SKCM_speed.ro")
    commonCell <- intersect(rownames(SKCM_mut), rownames(SKCM_speed))
    feature <- cbind(SKCM_mut[commonCell,], SKCM_speed[commonCell,])
  }
  if (feat=="v17_SKCM_methyl.gex"){
    load("data/pipeline/v17/bem/SKCM_methyl.ro")
    load("data/pipeline/v17/bem/SKCM_gex.ro")
    commonCell <- intersect(rownames(SKCM_methyl), rownames(SKCM_gex))
    feature <- cbind(SKCM_methyl[commonCell,], SKCM_gex[commonCell,])
  }
  if (feat=="v17_SKCM_methyl.cnv"){
    load("data/pipeline/v17/bem/SKCM_methyl.ro")
    load("data/pipeline/v17/bem/SKCM_cnv.ro")
    commonCell <- intersect(rownames(SKCM_methyl), rownames(SKCM_cnv))
    feature <- cbind(SKCM_methyl[commonCell,], SKCM_cnv[commonCell,])
  }
  if (feat=="v17_SKCM_methyl.speed"){
    load("data/pipeline/v17/bem/SKCM_methyl.ro")
    load("data/pipeline/v17/bem/SKCM_speed.ro")
    commonCell <- intersect(rownames(SKCM_methyl), rownames(SKCM_speed))
    feature <- cbind(SKCM_methyl[commonCell,], SKCM_speed[commonCell,])
  }
  if (feat=="v17_SKCM_gex.cnv"){
    load("data/pipeline/v17/bem/SKCM_gex.ro")
    load("data/pipeline/v17/bem/SKCM_cnv.ro")
    commonCell <- intersect(rownames(SKCM_gex), rownames(SKCM_cnv))
    feature <- cbind(SKCM_gex[commonCell,], SKCM_cnv[commonCell,])
  }
  if (feat=="v17_SKCM_gex.speed"){
    load("data/pipeline/v17/bem/SKCM_gex.ro")
    load("data/pipeline/v17/bem/SKCM_speed.ro")
    commonCell <- intersect(rownames(SKCM_gex), rownames(SKCM_speed))
    feature <- cbind(SKCM_gex[commonCell,], SKCM_speed[commonCell,])
  }
  if (feat=="v17_SKCM_cnv.speed"){
    load("data/pipeline/v17/bem/SKCM_cnv.ro")
    load("data/pipeline/v17/bem/SKCM_speed.ro")
    commonCell <- intersect(rownames(SKCM_cnv), rownames(SKCM_speed))
    feature <- cbind(SKCM_cnv[commonCell,], SKCM_speed[commonCell,])
  }
  if (feat=="v17_SKCM_mut.methyl.gex"){
    load("data/pipeline/v17/bem/SKCM_mut.ro")
    load("data/pipeline/v17/bem/SKCM_methyl.ro")
    load("data/pipeline/v17/bem/SKCM_gex.ro")
    commonCell <- intersect(rownames(SKCM_mut), intersect(rownames(SKCM_methyl), rownames(SKCM_gex)))
    feature <- cbind(SKCM_mut[commonCell,], SKCM_methyl[commonCell,], SKCM_gex[commonCell,])
  }
  if (feat=="v17_SKCM_mut.methyl.cnv"){
    load("data/pipeline/v17/bem/SKCM_mut.ro")
    load("data/pipeline/v17/bem/SKCM_methyl.ro")
    load("data/pipeline/v17/bem/SKCM_cnv.ro")
    commonCell <- intersect(rownames(SKCM_mut), intersect(rownames(SKCM_methyl), rownames(SKCM_cnv)))
    feature <- cbind(SKCM_mut[commonCell,], SKCM_methyl[commonCell,], SKCM_cnv[commonCell,])
  }
  if (feat=="v17_SKCM_mut.methyl.speed"){
    load("data/pipeline/v17/bem/SKCM_mut.ro")
    load("data/pipeline/v17/bem/SKCM_methyl.ro")
    load("data/pipeline/v17/bem/SKCM_speed.ro")
    commonCell <- intersect(rownames(SKCM_mut), intersect(rownames(SKCM_methyl), rownames(SKCM_speed)))
    feature <- cbind(SKCM_mut[commonCell,], SKCM_methyl[commonCell,], SKCM_speed[commonCell,])
  }
  if (feat=="v17_SKCM_mut.gex.cnv"){
    load("data/pipeline/v17/bem/SKCM_mut.ro")
    load("data/pipeline/v17/bem/SKCM_gex.ro")
    load("data/pipeline/v17/bem/SKCM_cnv.ro")
    commonCell <- intersect(rownames(SKCM_mut), intersect(rownames(SKCM_gex), rownames(SKCM_cnv)))
    feature <- cbind(SKCM_mut[commonCell,], SKCM_gex[commonCell,], SKCM_cnv[commonCell,])
  }
  if (feat=="v17_SKCM_mut.gex.speed"){
    load("data/pipeline/v17/bem/SKCM_mut.ro")
    load("data/pipeline/v17/bem/SKCM_gex.ro")
    load("data/pipeline/v17/bem/SKCM_speed.ro")
    commonCell <- intersect(rownames(SKCM_mut), intersect(rownames(SKCM_gex), rownames(SKCM_speed)))
    feature <- cbind(SKCM_mut[commonCell,], SKCM_gex[commonCell,], SKCM_speed[commonCell,])
  }
  if (feat=="v17_SKCM_mut.cnv.speed"){
    load("data/pipeline/v17/bem/SKCM_mut.ro")
    load("data/pipeline/v17/bem/SKCM_cnv.ro")
    load("data/pipeline/v17/bem/SKCM_speed.ro")
    commonCell <- intersect(rownames(SKCM_mut), intersect(rownames(SKCM_cnv), rownames(SKCM_speed)))
    feature <- cbind(SKCM_mut[commonCell,], SKCM_cnv[commonCell,], SKCM_speed[commonCell,])
  }
  if (feat=="v17_SKCM_methyl.gex.cnv"){
    load("data/pipeline/v17/bem/SKCM_methyl.ro")
    load("data/pipeline/v17/bem/SKCM_gex.ro")
    load("data/pipeline/v17/bem/SKCM_cnv.ro")
    commonCell <- intersect(rownames(SKCM_methyl), intersect(rownames(SKCM_gex), rownames(SKCM_cnv)))
    feature <- cbind(SKCM_methyl[commonCell,], SKCM_gex[commonCell,], SKCM_cnv[commonCell,])
  }
  if (feat=="v17_SKCM_methyl.gex.speed"){
    load("data/pipeline/v17/bem/SKCM_methyl.ro")
    load("data/pipeline/v17/bem/SKCM_gex.ro")
    load("data/pipeline/v17/bem/SKCM_speed.ro")
    commonCell <- intersect(rownames(SKCM_methyl), intersect(rownames(SKCM_gex), rownames(SKCM_speed)))
    feature <- cbind(SKCM_methyl[commonCell,], SKCM_gex[commonCell,], SKCM_speed[commonCell,])
  }
  if (feat=="v17_SKCM_methyl.cnv.speed"){
    load("data/pipeline/v17/bem/SKCM_methyl.ro")
    load("data/pipeline/v17/bem/SKCM_cnv.ro")
    load("data/pipeline/v17/bem/SKCM_speed.ro")
    commonCell <- intersect(rownames(SKCM_methyl), intersect(rownames(SKCM_cnv), rownames(SKCM_speed)))
    feature <- cbind(SKCM_methyl[commonCell,], SKCM_cnv[commonCell,], SKCM_speed[commonCell,])
  }
  if (feat=="v17_SKCM_gex.cnv.speed"){
    load("data/pipeline/v17/bem/SKCM_gex.ro")
    load("data/pipeline/v17/bem/SKCM_cnv.ro")
    load("data/pipeline/v17/bem/SKCM_speed.ro")
    commonCell <- intersect(rownames(SKCM_gex), intersect(rownames(SKCM_cnv), rownames(SKCM_speed)))
    feature <- cbind(SKCM_gex[commonCell,], SKCM_cnv[commonCell,], SKCM_speed[commonCell,])
  }
  if (feat=="v17_SKCM_mut.methyl.gex.cnv"){
    load("data/pipeline/v17/bem/SKCM_mut.ro")
    load("data/pipeline/v17/bem/SKCM_methyl.ro")
    load("data/pipeline/v17/bem/SKCM_gex.ro")
    load("data/pipeline/v17/bem/SKCM_cnv.ro")
    commonCell <- intersect(rownames(SKCM_mut), intersect(rownames(SKCM_methyl), intersect(rownames(SKCM_gex), rownames(SKCM_cnv))))
    feature <- cbind(SKCM_mut[commonCell,], SKCM_methyl[commonCell,], SKCM_gex[commonCell,], SKCM_cnv[commonCell,])
  }
  if (feat=="v17_SKCM_mut.methyl.gex.speed"){
    load("data/pipeline/v17/bem/SKCM_mut.ro")
    load("data/pipeline/v17/bem/SKCM_methyl.ro")
    load("data/pipeline/v17/bem/SKCM_gex.ro")
    load("data/pipeline/v17/bem/SKCM_speed.ro")
    commonCell <- intersect(rownames(SKCM_mut), intersect(rownames(SKCM_methyl), intersect(rownames(SKCM_gex), rownames(SKCM_speed))))
    feature <- cbind(SKCM_mut[commonCell,], SKCM_methyl[commonCell,], SKCM_gex[commonCell,], SKCM_speed[commonCell,])
  }
  if (feat=="v17_SKCM_mut.methyl.cnv.speed"){
    load("data/pipeline/v17/bem/SKCM_mut.ro")
    load("data/pipeline/v17/bem/SKCM_methyl.ro")
    load("data/pipeline/v17/bem/SKCM_cnv.ro")
    load("data/pipeline/v17/bem/SKCM_speed.ro")
    commonCell <- intersect(rownames(SKCM_mut), intersect(rownames(SKCM_methyl), intersect(rownames(SKCM_cnv), rownames(SKCM_speed))))
    feature <- cbind(SKCM_mut[commonCell,], SKCM_methyl[commonCell,], SKCM_cnv[commonCell,], SKCM_speed[commonCell,])
  }
  if (feat=="v17_SKCM_mut.gex.cnv.speed"){
    load("data/pipeline/v17/bem/SKCM_mut.ro")
    load("data/pipeline/v17/bem/SKCM_gex.ro")
    load("data/pipeline/v17/bem/SKCM_cnv.ro")
    load("data/pipeline/v17/bem/SKCM_speed.ro")
    commonCell <- intersect(rownames(SKCM_mut), intersect(rownames(SKCM_gex), intersect(rownames(SKCM_cnv), rownames(SKCM_speed))))
    feature <- cbind(SKCM_mut[commonCell,], SKCM_gex[commonCell,], SKCM_cnv[commonCell,], SKCM_speed[commonCell,])
  }
  if (feat=="v17_SKCM_methyl.gex.cnv.speed"){
    load("data/pipeline/v17/bem/SKCM_methyl.ro")
    load("data/pipeline/v17/bem/SKCM_gex.ro")
    load("data/pipeline/v17/bem/SKCM_cnv.ro")
    load("data/pipeline/v17/bem/SKCM_speed.ro")
    commonCell <- intersect(rownames(SKCM_methyl), intersect(rownames(SKCM_gex), intersect(rownames(SKCM_cnv), rownames(SKCM_speed))))
    feature <- cbind(SKCM_methyl[commonCell,], SKCM_gex[commonCell,], SKCM_cnv[commonCell,], SKCM_speed[commonCell,])
  }
  if (feat=="v17_SKCM_mut.methyl.gex.cnv.speed"){
    load("data/pipeline/v17/bem/SKCM_mut.ro")
    load("data/pipeline/v17/bem/SKCM_methyl.ro")
    load("data/pipeline/v17/bem/SKCM_gex.ro")
    load("data/pipeline/v17/bem/SKCM_cnv.ro")
    load("data/pipeline/v17/bem/SKCM_speed.ro")
    commonCell <- intersect(rownames(SKCM_mut), intersect(rownames(SKCM_methyl), intersect(rownames(SKCM_gex), intersect(rownames(SKCM_cnv), rownames(SKCM_speed)))))
    feature <- cbind(SKCM_mut[commonCell,], SKCM_methyl[commonCell,], SKCM_gex[commonCell,], SKCM_cnv[commonCell,], SKCM_speed[commonCell,])
  }
  if (feat=="v17_STAD_mut"){
    load("data/pipeline/v17/bem/STAD_mut.ro")
    commonCell <- rownames(STAD_mut)
    feature <- cbind(STAD_mut[commonCell,])
  }
  if (feat=="v17_STAD_methyl"){
    load("data/pipeline/v17/bem/STAD_methyl.ro")
    commonCell <- rownames(STAD_methyl)
    feature <- cbind(STAD_methyl[commonCell,])
  }
  if (feat=="v17_STAD_gex"){
    load("data/pipeline/v17/bem/STAD_gex.ro")
    commonCell <- rownames(STAD_gex)
    feature <- cbind(STAD_gex[commonCell,])
  }
  if (feat=="v17_STAD_cnv"){
    load("data/pipeline/v17/bem/STAD_cnv.ro")
    commonCell <- rownames(STAD_cnv)
    feature <- cbind(STAD_cnv[commonCell,])
  }
  if (feat=="v17_STAD_speed"){
    load("data/pipeline/v17/bem/STAD_speed.ro")
    commonCell <- rownames(STAD_speed)
    feature <- cbind(STAD_speed[commonCell,])
  }
  if (feat=="v17_STAD_mut.methyl"){
    load("data/pipeline/v17/bem/STAD_mut.ro")
    load("data/pipeline/v17/bem/STAD_methyl.ro")
    commonCell <- intersect(rownames(STAD_mut), rownames(STAD_methyl))
    feature <- cbind(STAD_mut[commonCell,], STAD_methyl[commonCell,])
  }
  if (feat=="v17_STAD_mut.gex"){
    load("data/pipeline/v17/bem/STAD_mut.ro")
    load("data/pipeline/v17/bem/STAD_gex.ro")
    commonCell <- intersect(rownames(STAD_mut), rownames(STAD_gex))
    feature <- cbind(STAD_mut[commonCell,], STAD_gex[commonCell,])
  }
  if (feat=="v17_STAD_mut.cnv"){
    load("data/pipeline/v17/bem/STAD_mut.ro")
    load("data/pipeline/v17/bem/STAD_cnv.ro")
    commonCell <- intersect(rownames(STAD_mut), rownames(STAD_cnv))
    feature <- cbind(STAD_mut[commonCell,], STAD_cnv[commonCell,])
  }
  if (feat=="v17_STAD_mut.speed"){
    load("data/pipeline/v17/bem/STAD_mut.ro")
    load("data/pipeline/v17/bem/STAD_speed.ro")
    commonCell <- intersect(rownames(STAD_mut), rownames(STAD_speed))
    feature <- cbind(STAD_mut[commonCell,], STAD_speed[commonCell,])
  }
  if (feat=="v17_STAD_methyl.gex"){
    load("data/pipeline/v17/bem/STAD_methyl.ro")
    load("data/pipeline/v17/bem/STAD_gex.ro")
    commonCell <- intersect(rownames(STAD_methyl), rownames(STAD_gex))
    feature <- cbind(STAD_methyl[commonCell,], STAD_gex[commonCell,])
  }
  if (feat=="v17_STAD_methyl.cnv"){
    load("data/pipeline/v17/bem/STAD_methyl.ro")
    load("data/pipeline/v17/bem/STAD_cnv.ro")
    commonCell <- intersect(rownames(STAD_methyl), rownames(STAD_cnv))
    feature <- cbind(STAD_methyl[commonCell,], STAD_cnv[commonCell,])
  }
  if (feat=="v17_STAD_methyl.speed"){
    load("data/pipeline/v17/bem/STAD_methyl.ro")
    load("data/pipeline/v17/bem/STAD_speed.ro")
    commonCell <- intersect(rownames(STAD_methyl), rownames(STAD_speed))
    feature <- cbind(STAD_methyl[commonCell,], STAD_speed[commonCell,])
  }
  if (feat=="v17_STAD_gex.cnv"){
    load("data/pipeline/v17/bem/STAD_gex.ro")
    load("data/pipeline/v17/bem/STAD_cnv.ro")
    commonCell <- intersect(rownames(STAD_gex), rownames(STAD_cnv))
    feature <- cbind(STAD_gex[commonCell,], STAD_cnv[commonCell,])
  }
  if (feat=="v17_STAD_gex.speed"){
    load("data/pipeline/v17/bem/STAD_gex.ro")
    load("data/pipeline/v17/bem/STAD_speed.ro")
    commonCell <- intersect(rownames(STAD_gex), rownames(STAD_speed))
    feature <- cbind(STAD_gex[commonCell,], STAD_speed[commonCell,])
  }
  if (feat=="v17_STAD_cnv.speed"){
    load("data/pipeline/v17/bem/STAD_cnv.ro")
    load("data/pipeline/v17/bem/STAD_speed.ro")
    commonCell <- intersect(rownames(STAD_cnv), rownames(STAD_speed))
    feature <- cbind(STAD_cnv[commonCell,], STAD_speed[commonCell,])
  }
  if (feat=="v17_STAD_mut.methyl.gex"){
    load("data/pipeline/v17/bem/STAD_mut.ro")
    load("data/pipeline/v17/bem/STAD_methyl.ro")
    load("data/pipeline/v17/bem/STAD_gex.ro")
    commonCell <- intersect(rownames(STAD_mut), intersect(rownames(STAD_methyl), rownames(STAD_gex)))
    feature <- cbind(STAD_mut[commonCell,], STAD_methyl[commonCell,], STAD_gex[commonCell,])
  }
  if (feat=="v17_STAD_mut.methyl.cnv"){
    load("data/pipeline/v17/bem/STAD_mut.ro")
    load("data/pipeline/v17/bem/STAD_methyl.ro")
    load("data/pipeline/v17/bem/STAD_cnv.ro")
    commonCell <- intersect(rownames(STAD_mut), intersect(rownames(STAD_methyl), rownames(STAD_cnv)))
    feature <- cbind(STAD_mut[commonCell,], STAD_methyl[commonCell,], STAD_cnv[commonCell,])
  }
  if (feat=="v17_STAD_mut.methyl.speed"){
    load("data/pipeline/v17/bem/STAD_mut.ro")
    load("data/pipeline/v17/bem/STAD_methyl.ro")
    load("data/pipeline/v17/bem/STAD_speed.ro")
    commonCell <- intersect(rownames(STAD_mut), intersect(rownames(STAD_methyl), rownames(STAD_speed)))
    feature <- cbind(STAD_mut[commonCell,], STAD_methyl[commonCell,], STAD_speed[commonCell,])
  }
  if (feat=="v17_STAD_mut.gex.cnv"){
    load("data/pipeline/v17/bem/STAD_mut.ro")
    load("data/pipeline/v17/bem/STAD_gex.ro")
    load("data/pipeline/v17/bem/STAD_cnv.ro")
    commonCell <- intersect(rownames(STAD_mut), intersect(rownames(STAD_gex), rownames(STAD_cnv)))
    feature <- cbind(STAD_mut[commonCell,], STAD_gex[commonCell,], STAD_cnv[commonCell,])
  }
  if (feat=="v17_STAD_mut.gex.speed"){
    load("data/pipeline/v17/bem/STAD_mut.ro")
    load("data/pipeline/v17/bem/STAD_gex.ro")
    load("data/pipeline/v17/bem/STAD_speed.ro")
    commonCell <- intersect(rownames(STAD_mut), intersect(rownames(STAD_gex), rownames(STAD_speed)))
    feature <- cbind(STAD_mut[commonCell,], STAD_gex[commonCell,], STAD_speed[commonCell,])
  }
  if (feat=="v17_STAD_mut.cnv.speed"){
    load("data/pipeline/v17/bem/STAD_mut.ro")
    load("data/pipeline/v17/bem/STAD_cnv.ro")
    load("data/pipeline/v17/bem/STAD_speed.ro")
    commonCell <- intersect(rownames(STAD_mut), intersect(rownames(STAD_cnv), rownames(STAD_speed)))
    feature <- cbind(STAD_mut[commonCell,], STAD_cnv[commonCell,], STAD_speed[commonCell,])
  }
  if (feat=="v17_STAD_methyl.gex.cnv"){
    load("data/pipeline/v17/bem/STAD_methyl.ro")
    load("data/pipeline/v17/bem/STAD_gex.ro")
    load("data/pipeline/v17/bem/STAD_cnv.ro")
    commonCell <- intersect(rownames(STAD_methyl), intersect(rownames(STAD_gex), rownames(STAD_cnv)))
    feature <- cbind(STAD_methyl[commonCell,], STAD_gex[commonCell,], STAD_cnv[commonCell,])
  }
  if (feat=="v17_STAD_methyl.gex.speed"){
    load("data/pipeline/v17/bem/STAD_methyl.ro")
    load("data/pipeline/v17/bem/STAD_gex.ro")
    load("data/pipeline/v17/bem/STAD_speed.ro")
    commonCell <- intersect(rownames(STAD_methyl), intersect(rownames(STAD_gex), rownames(STAD_speed)))
    feature <- cbind(STAD_methyl[commonCell,], STAD_gex[commonCell,], STAD_speed[commonCell,])
  }
  if (feat=="v17_STAD_methyl.cnv.speed"){
    load("data/pipeline/v17/bem/STAD_methyl.ro")
    load("data/pipeline/v17/bem/STAD_cnv.ro")
    load("data/pipeline/v17/bem/STAD_speed.ro")
    commonCell <- intersect(rownames(STAD_methyl), intersect(rownames(STAD_cnv), rownames(STAD_speed)))
    feature <- cbind(STAD_methyl[commonCell,], STAD_cnv[commonCell,], STAD_speed[commonCell,])
  }
  if (feat=="v17_STAD_gex.cnv.speed"){
    load("data/pipeline/v17/bem/STAD_gex.ro")
    load("data/pipeline/v17/bem/STAD_cnv.ro")
    load("data/pipeline/v17/bem/STAD_speed.ro")
    commonCell <- intersect(rownames(STAD_gex), intersect(rownames(STAD_cnv), rownames(STAD_speed)))
    feature <- cbind(STAD_gex[commonCell,], STAD_cnv[commonCell,], STAD_speed[commonCell,])
  }
  if (feat=="v17_STAD_mut.methyl.gex.cnv"){
    load("data/pipeline/v17/bem/STAD_mut.ro")
    load("data/pipeline/v17/bem/STAD_methyl.ro")
    load("data/pipeline/v17/bem/STAD_gex.ro")
    load("data/pipeline/v17/bem/STAD_cnv.ro")
    commonCell <- intersect(rownames(STAD_mut), intersect(rownames(STAD_methyl), intersect(rownames(STAD_gex), rownames(STAD_cnv))))
    feature <- cbind(STAD_mut[commonCell,], STAD_methyl[commonCell,], STAD_gex[commonCell,], STAD_cnv[commonCell,])
  }
  if (feat=="v17_STAD_mut.methyl.gex.speed"){
    load("data/pipeline/v17/bem/STAD_mut.ro")
    load("data/pipeline/v17/bem/STAD_methyl.ro")
    load("data/pipeline/v17/bem/STAD_gex.ro")
    load("data/pipeline/v17/bem/STAD_speed.ro")
    commonCell <- intersect(rownames(STAD_mut), intersect(rownames(STAD_methyl), intersect(rownames(STAD_gex), rownames(STAD_speed))))
    feature <- cbind(STAD_mut[commonCell,], STAD_methyl[commonCell,], STAD_gex[commonCell,], STAD_speed[commonCell,])
  }
  if (feat=="v17_STAD_mut.methyl.cnv.speed"){
    load("data/pipeline/v17/bem/STAD_mut.ro")
    load("data/pipeline/v17/bem/STAD_methyl.ro")
    load("data/pipeline/v17/bem/STAD_cnv.ro")
    load("data/pipeline/v17/bem/STAD_speed.ro")
    commonCell <- intersect(rownames(STAD_mut), intersect(rownames(STAD_methyl), intersect(rownames(STAD_cnv), rownames(STAD_speed))))
    feature <- cbind(STAD_mut[commonCell,], STAD_methyl[commonCell,], STAD_cnv[commonCell,], STAD_speed[commonCell,])
  }
  if (feat=="v17_STAD_mut.gex.cnv.speed"){
    load("data/pipeline/v17/bem/STAD_mut.ro")
    load("data/pipeline/v17/bem/STAD_gex.ro")
    load("data/pipeline/v17/bem/STAD_cnv.ro")
    load("data/pipeline/v17/bem/STAD_speed.ro")
    commonCell <- intersect(rownames(STAD_mut), intersect(rownames(STAD_gex), intersect(rownames(STAD_cnv), rownames(STAD_speed))))
    feature <- cbind(STAD_mut[commonCell,], STAD_gex[commonCell,], STAD_cnv[commonCell,], STAD_speed[commonCell,])
  }
  if (feat=="v17_STAD_methyl.gex.cnv.speed"){
    load("data/pipeline/v17/bem/STAD_methyl.ro")
    load("data/pipeline/v17/bem/STAD_gex.ro")
    load("data/pipeline/v17/bem/STAD_cnv.ro")
    load("data/pipeline/v17/bem/STAD_speed.ro")
    commonCell <- intersect(rownames(STAD_methyl), intersect(rownames(STAD_gex), intersect(rownames(STAD_cnv), rownames(STAD_speed))))
    feature <- cbind(STAD_methyl[commonCell,], STAD_gex[commonCell,], STAD_cnv[commonCell,], STAD_speed[commonCell,])
  }
  if (feat=="v17_STAD_mut.methyl.gex.cnv.speed"){
    load("data/pipeline/v17/bem/STAD_mut.ro")
    load("data/pipeline/v17/bem/STAD_methyl.ro")
    load("data/pipeline/v17/bem/STAD_gex.ro")
    load("data/pipeline/v17/bem/STAD_cnv.ro")
    load("data/pipeline/v17/bem/STAD_speed.ro")
    commonCell <- intersect(rownames(STAD_mut), intersect(rownames(STAD_methyl), intersect(rownames(STAD_gex), intersect(rownames(STAD_cnv), rownames(STAD_speed)))))
    feature <- cbind(STAD_mut[commonCell,], STAD_methyl[commonCell,], STAD_gex[commonCell,], STAD_cnv[commonCell,], STAD_speed[commonCell,])
  }
  if (feat=="v17_THCA_mut"){
    load("data/pipeline/v17/bem/THCA_mut.ro")
    commonCell <- rownames(THCA_mut)
    feature <- cbind(THCA_mut[commonCell,])
  }
  if (feat=="v17_THCA_methyl"){
    load("data/pipeline/v17/bem/THCA_methyl.ro")
    commonCell <- rownames(THCA_methyl)
    feature <- cbind(THCA_methyl[commonCell,])
  }
  if (feat=="v17_THCA_gex"){
    load("data/pipeline/v17/bem/THCA_gex.ro")
    commonCell <- rownames(THCA_gex)
    feature <- cbind(THCA_gex[commonCell,])
  }
  if (feat=="v17_THCA_cnv"){
    load("data/pipeline/v17/bem/THCA_cnv.ro")
    commonCell <- rownames(THCA_cnv)
    feature <- cbind(THCA_cnv[commonCell,])
  }
  if (feat=="v17_THCA_speed"){
    load("data/pipeline/v17/bem/THCA_speed.ro")
    commonCell <- rownames(THCA_speed)
    feature <- cbind(THCA_speed[commonCell,])
  }
  if (feat=="v17_THCA_mut.methyl"){
    load("data/pipeline/v17/bem/THCA_mut.ro")
    load("data/pipeline/v17/bem/THCA_methyl.ro")
    commonCell <- intersect(rownames(THCA_mut), rownames(THCA_methyl))
    feature <- cbind(THCA_mut[commonCell,], THCA_methyl[commonCell,])
  }
  if (feat=="v17_THCA_mut.gex"){
    load("data/pipeline/v17/bem/THCA_mut.ro")
    load("data/pipeline/v17/bem/THCA_gex.ro")
    commonCell <- intersect(rownames(THCA_mut), rownames(THCA_gex))
    feature <- cbind(THCA_mut[commonCell,], THCA_gex[commonCell,])
  }
  if (feat=="v17_THCA_mut.cnv"){
    load("data/pipeline/v17/bem/THCA_mut.ro")
    load("data/pipeline/v17/bem/THCA_cnv.ro")
    commonCell <- intersect(rownames(THCA_mut), rownames(THCA_cnv))
    feature <- cbind(THCA_mut[commonCell,], THCA_cnv[commonCell,])
  }
  if (feat=="v17_THCA_mut.speed"){
    load("data/pipeline/v17/bem/THCA_mut.ro")
    load("data/pipeline/v17/bem/THCA_speed.ro")
    commonCell <- intersect(rownames(THCA_mut), rownames(THCA_speed))
    feature <- cbind(THCA_mut[commonCell,], THCA_speed[commonCell,])
  }
  if (feat=="v17_THCA_methyl.gex"){
    load("data/pipeline/v17/bem/THCA_methyl.ro")
    load("data/pipeline/v17/bem/THCA_gex.ro")
    commonCell <- intersect(rownames(THCA_methyl), rownames(THCA_gex))
    feature <- cbind(THCA_methyl[commonCell,], THCA_gex[commonCell,])
  }
  if (feat=="v17_THCA_methyl.cnv"){
    load("data/pipeline/v17/bem/THCA_methyl.ro")
    load("data/pipeline/v17/bem/THCA_cnv.ro")
    commonCell <- intersect(rownames(THCA_methyl), rownames(THCA_cnv))
    feature <- cbind(THCA_methyl[commonCell,], THCA_cnv[commonCell,])
  }
  if (feat=="v17_THCA_methyl.speed"){
    load("data/pipeline/v17/bem/THCA_methyl.ro")
    load("data/pipeline/v17/bem/THCA_speed.ro")
    commonCell <- intersect(rownames(THCA_methyl), rownames(THCA_speed))
    feature <- cbind(THCA_methyl[commonCell,], THCA_speed[commonCell,])
  }
  if (feat=="v17_THCA_gex.cnv"){
    load("data/pipeline/v17/bem/THCA_gex.ro")
    load("data/pipeline/v17/bem/THCA_cnv.ro")
    commonCell <- intersect(rownames(THCA_gex), rownames(THCA_cnv))
    feature <- cbind(THCA_gex[commonCell,], THCA_cnv[commonCell,])
  }
  if (feat=="v17_THCA_gex.speed"){
    load("data/pipeline/v17/bem/THCA_gex.ro")
    load("data/pipeline/v17/bem/THCA_speed.ro")
    commonCell <- intersect(rownames(THCA_gex), rownames(THCA_speed))
    feature <- cbind(THCA_gex[commonCell,], THCA_speed[commonCell,])
  }
  if (feat=="v17_THCA_cnv.speed"){
    load("data/pipeline/v17/bem/THCA_cnv.ro")
    load("data/pipeline/v17/bem/THCA_speed.ro")
    commonCell <- intersect(rownames(THCA_cnv), rownames(THCA_speed))
    feature <- cbind(THCA_cnv[commonCell,], THCA_speed[commonCell,])
  }
  if (feat=="v17_THCA_mut.methyl.gex"){
    load("data/pipeline/v17/bem/THCA_mut.ro")
    load("data/pipeline/v17/bem/THCA_methyl.ro")
    load("data/pipeline/v17/bem/THCA_gex.ro")
    commonCell <- intersect(rownames(THCA_mut), intersect(rownames(THCA_methyl), rownames(THCA_gex)))
    feature <- cbind(THCA_mut[commonCell,], THCA_methyl[commonCell,], THCA_gex[commonCell,])
  }
  if (feat=="v17_THCA_mut.methyl.cnv"){
    load("data/pipeline/v17/bem/THCA_mut.ro")
    load("data/pipeline/v17/bem/THCA_methyl.ro")
    load("data/pipeline/v17/bem/THCA_cnv.ro")
    commonCell <- intersect(rownames(THCA_mut), intersect(rownames(THCA_methyl), rownames(THCA_cnv)))
    feature <- cbind(THCA_mut[commonCell,], THCA_methyl[commonCell,], THCA_cnv[commonCell,])
  }
  if (feat=="v17_THCA_mut.methyl.speed"){
    load("data/pipeline/v17/bem/THCA_mut.ro")
    load("data/pipeline/v17/bem/THCA_methyl.ro")
    load("data/pipeline/v17/bem/THCA_speed.ro")
    commonCell <- intersect(rownames(THCA_mut), intersect(rownames(THCA_methyl), rownames(THCA_speed)))
    feature <- cbind(THCA_mut[commonCell,], THCA_methyl[commonCell,], THCA_speed[commonCell,])
  }
  if (feat=="v17_THCA_mut.gex.cnv"){
    load("data/pipeline/v17/bem/THCA_mut.ro")
    load("data/pipeline/v17/bem/THCA_gex.ro")
    load("data/pipeline/v17/bem/THCA_cnv.ro")
    commonCell <- intersect(rownames(THCA_mut), intersect(rownames(THCA_gex), rownames(THCA_cnv)))
    feature <- cbind(THCA_mut[commonCell,], THCA_gex[commonCell,], THCA_cnv[commonCell,])
  }
  if (feat=="v17_THCA_mut.gex.speed"){
    load("data/pipeline/v17/bem/THCA_mut.ro")
    load("data/pipeline/v17/bem/THCA_gex.ro")
    load("data/pipeline/v17/bem/THCA_speed.ro")
    commonCell <- intersect(rownames(THCA_mut), intersect(rownames(THCA_gex), rownames(THCA_speed)))
    feature <- cbind(THCA_mut[commonCell,], THCA_gex[commonCell,], THCA_speed[commonCell,])
  }
  if (feat=="v17_THCA_mut.cnv.speed"){
    load("data/pipeline/v17/bem/THCA_mut.ro")
    load("data/pipeline/v17/bem/THCA_cnv.ro")
    load("data/pipeline/v17/bem/THCA_speed.ro")
    commonCell <- intersect(rownames(THCA_mut), intersect(rownames(THCA_cnv), rownames(THCA_speed)))
    feature <- cbind(THCA_mut[commonCell,], THCA_cnv[commonCell,], THCA_speed[commonCell,])
  }
  if (feat=="v17_THCA_methyl.gex.cnv"){
    load("data/pipeline/v17/bem/THCA_methyl.ro")
    load("data/pipeline/v17/bem/THCA_gex.ro")
    load("data/pipeline/v17/bem/THCA_cnv.ro")
    commonCell <- intersect(rownames(THCA_methyl), intersect(rownames(THCA_gex), rownames(THCA_cnv)))
    feature <- cbind(THCA_methyl[commonCell,], THCA_gex[commonCell,], THCA_cnv[commonCell,])
  }
  if (feat=="v17_THCA_methyl.gex.speed"){
    load("data/pipeline/v17/bem/THCA_methyl.ro")
    load("data/pipeline/v17/bem/THCA_gex.ro")
    load("data/pipeline/v17/bem/THCA_speed.ro")
    commonCell <- intersect(rownames(THCA_methyl), intersect(rownames(THCA_gex), rownames(THCA_speed)))
    feature <- cbind(THCA_methyl[commonCell,], THCA_gex[commonCell,], THCA_speed[commonCell,])
  }
  if (feat=="v17_THCA_methyl.cnv.speed"){
    load("data/pipeline/v17/bem/THCA_methyl.ro")
    load("data/pipeline/v17/bem/THCA_cnv.ro")
    load("data/pipeline/v17/bem/THCA_speed.ro")
    commonCell <- intersect(rownames(THCA_methyl), intersect(rownames(THCA_cnv), rownames(THCA_speed)))
    feature <- cbind(THCA_methyl[commonCell,], THCA_cnv[commonCell,], THCA_speed[commonCell,])
  }
  if (feat=="v17_THCA_gex.cnv.speed"){
    load("data/pipeline/v17/bem/THCA_gex.ro")
    load("data/pipeline/v17/bem/THCA_cnv.ro")
    load("data/pipeline/v17/bem/THCA_speed.ro")
    commonCell <- intersect(rownames(THCA_gex), intersect(rownames(THCA_cnv), rownames(THCA_speed)))
    feature <- cbind(THCA_gex[commonCell,], THCA_cnv[commonCell,], THCA_speed[commonCell,])
  }
  if (feat=="v17_THCA_mut.methyl.gex.cnv"){
    load("data/pipeline/v17/bem/THCA_mut.ro")
    load("data/pipeline/v17/bem/THCA_methyl.ro")
    load("data/pipeline/v17/bem/THCA_gex.ro")
    load("data/pipeline/v17/bem/THCA_cnv.ro")
    commonCell <- intersect(rownames(THCA_mut), intersect(rownames(THCA_methyl), intersect(rownames(THCA_gex), rownames(THCA_cnv))))
    feature <- cbind(THCA_mut[commonCell,], THCA_methyl[commonCell,], THCA_gex[commonCell,], THCA_cnv[commonCell,])
  }
  if (feat=="v17_THCA_mut.methyl.gex.speed"){
    load("data/pipeline/v17/bem/THCA_mut.ro")
    load("data/pipeline/v17/bem/THCA_methyl.ro")
    load("data/pipeline/v17/bem/THCA_gex.ro")
    load("data/pipeline/v17/bem/THCA_speed.ro")
    commonCell <- intersect(rownames(THCA_mut), intersect(rownames(THCA_methyl), intersect(rownames(THCA_gex), rownames(THCA_speed))))
    feature <- cbind(THCA_mut[commonCell,], THCA_methyl[commonCell,], THCA_gex[commonCell,], THCA_speed[commonCell,])
  }
  if (feat=="v17_THCA_mut.methyl.cnv.speed"){
    load("data/pipeline/v17/bem/THCA_mut.ro")
    load("data/pipeline/v17/bem/THCA_methyl.ro")
    load("data/pipeline/v17/bem/THCA_cnv.ro")
    load("data/pipeline/v17/bem/THCA_speed.ro")
    commonCell <- intersect(rownames(THCA_mut), intersect(rownames(THCA_methyl), intersect(rownames(THCA_cnv), rownames(THCA_speed))))
    feature <- cbind(THCA_mut[commonCell,], THCA_methyl[commonCell,], THCA_cnv[commonCell,], THCA_speed[commonCell,])
  }
  if (feat=="v17_THCA_mut.gex.cnv.speed"){
    load("data/pipeline/v17/bem/THCA_mut.ro")
    load("data/pipeline/v17/bem/THCA_gex.ro")
    load("data/pipeline/v17/bem/THCA_cnv.ro")
    load("data/pipeline/v17/bem/THCA_speed.ro")
    commonCell <- intersect(rownames(THCA_mut), intersect(rownames(THCA_gex), intersect(rownames(THCA_cnv), rownames(THCA_speed))))
    feature <- cbind(THCA_mut[commonCell,], THCA_gex[commonCell,], THCA_cnv[commonCell,], THCA_speed[commonCell,])
  }
  if (feat=="v17_THCA_methyl.gex.cnv.speed"){
    load("data/pipeline/v17/bem/THCA_methyl.ro")
    load("data/pipeline/v17/bem/THCA_gex.ro")
    load("data/pipeline/v17/bem/THCA_cnv.ro")
    load("data/pipeline/v17/bem/THCA_speed.ro")
    commonCell <- intersect(rownames(THCA_methyl), intersect(rownames(THCA_gex), intersect(rownames(THCA_cnv), rownames(THCA_speed))))
    feature <- cbind(THCA_methyl[commonCell,], THCA_gex[commonCell,], THCA_cnv[commonCell,], THCA_speed[commonCell,])
  }
  if (feat=="v17_THCA_mut.methyl.gex.cnv.speed"){
    load("data/pipeline/v17/bem/THCA_mut.ro")
    load("data/pipeline/v17/bem/THCA_methyl.ro")
    load("data/pipeline/v17/bem/THCA_gex.ro")
    load("data/pipeline/v17/bem/THCA_cnv.ro")
    load("data/pipeline/v17/bem/THCA_speed.ro")
    commonCell <- intersect(rownames(THCA_mut), intersect(rownames(THCA_methyl), intersect(rownames(THCA_gex), intersect(rownames(THCA_cnv), rownames(THCA_speed)))))
    feature <- cbind(THCA_mut[commonCell,], THCA_methyl[commonCell,], THCA_gex[commonCell,], THCA_cnv[commonCell,], THCA_speed[commonCell,])
  }
  if (feat=="v17_PANCAN_mut"){
    load("data/pipeline/v17/bem/PANCAN_mut.ro")
    commonCell <- rownames(PANCAN_mut)
    feature <- cbind(PANCAN_mut[commonCell,])
  }
  if (feat=="v17_PANCAN_methyl"){
    load("data/pipeline/v17/bem/PANCAN_methyl.ro")
    commonCell <- rownames(PANCAN_methyl)
    feature <- cbind(PANCAN_methyl[commonCell,])
  }
  if (feat=="v17_PANCAN_gex"){
    load("data/pipeline/v17/bem/PANCAN_gex.ro")
    commonCell <- rownames(PANCAN_gex)
    feature <- cbind(PANCAN_gex[commonCell,])
  }
  if (feat=="v17_PANCAN_cnv"){
    load("data/pipeline/v17/bem/PANCAN_cnv.ro")
    commonCell <- rownames(PANCAN_cnv)
    feature <- cbind(PANCAN_cnv[commonCell,])
  }
  if (feat=="v17_PANCAN_speed"){
    load("data/pipeline/v17/bem/PANCAN_speed.ro")
    commonCell <- rownames(PANCAN_speed)
    feature <- cbind(PANCAN_speed[commonCell,])
  }
  if (feat=="v17_PANCAN_tissue"){
    load("data/pipeline/v17/tissue.ro")
    commonCell <- rownames(tissue)
    feature <- cbind(tissue[commonCell,])
  }
  if (feat=="v17_PANCAN_mut.methyl"){
    load("data/pipeline/v17/bem/PANCAN_mut.ro")
    load("data/pipeline/v17/bem/PANCAN_methyl.ro")
    commonCell <- intersect(rownames(PANCAN_mut), rownames(PANCAN_methyl))
    feature <- cbind(PANCAN_mut[commonCell,], PANCAN_methyl[commonCell,])
  }
  if (feat=="v17_PANCAN_mut.gex"){
    load("data/pipeline/v17/bem/PANCAN_mut.ro")
    load("data/pipeline/v17/bem/PANCAN_gex.ro")
    commonCell <- intersect(rownames(PANCAN_mut), rownames(PANCAN_gex))
    feature <- cbind(PANCAN_mut[commonCell,], PANCAN_gex[commonCell,])
  }
  if (feat=="v17_PANCAN_mut.cnv"){
    load("data/pipeline/v17/bem/PANCAN_mut.ro")
    load("data/pipeline/v17/bem/PANCAN_cnv.ro")
    commonCell <- intersect(rownames(PANCAN_mut), rownames(PANCAN_cnv))
    feature <- cbind(PANCAN_mut[commonCell,], PANCAN_cnv[commonCell,])
  }
  if (feat=="v17_PANCAN_mut.speed"){
    load("data/pipeline/v17/bem/PANCAN_mut.ro")
    load("data/pipeline/v17/bem/PANCAN_speed.ro")
    commonCell <- intersect(rownames(PANCAN_mut), rownames(PANCAN_speed))
    feature <- cbind(PANCAN_mut[commonCell,], PANCAN_speed[commonCell,])
  }
  if (feat=="v17_PANCAN_mut.tissue"){
    load("data/pipeline/v17/bem/PANCAN_mut.ro")
    load("data/pipeline/v17/tissue.ro")
    commonCell <- intersect(rownames(PANCAN_mut), rownames(tissue))
    feature <- cbind(PANCAN_mut[commonCell,], tissue[commonCell,])
  }
  if (feat=="v17_PANCAN_methyl.gex"){
    load("data/pipeline/v17/bem/PANCAN_methyl.ro")
    load("data/pipeline/v17/bem/PANCAN_gex.ro")
    commonCell <- intersect(rownames(PANCAN_methyl), rownames(PANCAN_gex))
    feature <- cbind(PANCAN_methyl[commonCell,], PANCAN_gex[commonCell,])
  }
  if (feat=="v17_PANCAN_methyl.cnv"){
    load("data/pipeline/v17/bem/PANCAN_methyl.ro")
    load("data/pipeline/v17/bem/PANCAN_cnv.ro")
    commonCell <- intersect(rownames(PANCAN_methyl), rownames(PANCAN_cnv))
    feature <- cbind(PANCAN_methyl[commonCell,], PANCAN_cnv[commonCell,])
  }
  if (feat=="v17_PANCAN_methyl.speed"){
    load("data/pipeline/v17/bem/PANCAN_methyl.ro")
    load("data/pipeline/v17/bem/PANCAN_speed.ro")
    commonCell <- intersect(rownames(PANCAN_methyl), rownames(PANCAN_speed))
    feature <- cbind(PANCAN_methyl[commonCell,], PANCAN_speed[commonCell,])
  }
  if (feat=="v17_PANCAN_methyl.tissue"){
    load("data/pipeline/v17/bem/PANCAN_methyl.ro")
    load("data/pipeline/v17/tissue.ro")
    commonCell <- intersect(rownames(PANCAN_methyl), rownames(tissue))
    feature <- cbind(PANCAN_methyl[commonCell,], tissue[commonCell,])
  }
  if (feat=="v17_PANCAN_gex.cnv"){
    load("data/pipeline/v17/bem/PANCAN_gex.ro")
    load("data/pipeline/v17/bem/PANCAN_cnv.ro")
    commonCell <- intersect(rownames(PANCAN_gex), rownames(PANCAN_cnv))
    feature <- cbind(PANCAN_gex[commonCell,], PANCAN_cnv[commonCell,])
  }
  if (feat=="v17_PANCAN_gex.speed"){
    load("data/pipeline/v17/bem/PANCAN_gex.ro")
    load("data/pipeline/v17/bem/PANCAN_speed.ro")
    commonCell <- intersect(rownames(PANCAN_gex), rownames(PANCAN_speed))
    feature <- cbind(PANCAN_gex[commonCell,], PANCAN_speed[commonCell,])
  }
  if (feat=="v17_PANCAN_gex.tissue"){
    load("data/pipeline/v17/bem/PANCAN_gex.ro")
    load("data/pipeline/v17/tissue.ro")
    commonCell <- intersect(rownames(PANCAN_gex), rownames(tissue))
    feature <- cbind(PANCAN_gex[commonCell,], tissue[commonCell,])
  }
  if (feat=="v17_PANCAN_cnv.speed"){
    load("data/pipeline/v17/bem/PANCAN_cnv.ro")
    load("data/pipeline/v17/bem/PANCAN_speed.ro")
    commonCell <- intersect(rownames(PANCAN_cnv), rownames(PANCAN_speed))
    feature <- cbind(PANCAN_cnv[commonCell,], PANCAN_speed[commonCell,])
  }
  if (feat=="v17_PANCAN_cnv.tissue"){
    load("data/pipeline/v17/bem/PANCAN_cnv.ro")
    load("data/pipeline/v17/tissue.ro")
    commonCell <- intersect(rownames(PANCAN_cnv), rownames(tissue))
    feature <- cbind(PANCAN_cnv[commonCell,], tissue[commonCell,])
  }
  if (feat=="v17_PANCAN_speed.tissue"){
    load("data/pipeline/v17/bem/PANCAN_speed.ro")
    load("data/pipeline/v17/tissue.ro")
    commonCell <- intersect(rownames(PANCAN_speed), rownames(tissue))
    feature <- cbind(PANCAN_speed[commonCell,], tissue[commonCell,])
  }
  if (feat=="v17_PANCAN_mut.methyl.gex"){
    load("data/pipeline/v17/bem/PANCAN_mut.ro")
    load("data/pipeline/v17/bem/PANCAN_methyl.ro")
    load("data/pipeline/v17/bem/PANCAN_gex.ro")
    commonCell <- intersect(rownames(PANCAN_mut), intersect(rownames(PANCAN_methyl), rownames(PANCAN_gex)))
    feature <- cbind(PANCAN_mut[commonCell,], PANCAN_methyl[commonCell,], PANCAN_gex[commonCell,])
  }
  if (feat=="v17_PANCAN_mut.methyl.cnv"){
    load("data/pipeline/v17/bem/PANCAN_mut.ro")
    load("data/pipeline/v17/bem/PANCAN_methyl.ro")
    load("data/pipeline/v17/bem/PANCAN_cnv.ro")
    commonCell <- intersect(rownames(PANCAN_mut), intersect(rownames(PANCAN_methyl), rownames(PANCAN_cnv)))
    feature <- cbind(PANCAN_mut[commonCell,], PANCAN_methyl[commonCell,], PANCAN_cnv[commonCell,])
  }
  if (feat=="v17_PANCAN_mut.methyl.speed"){
    load("data/pipeline/v17/bem/PANCAN_mut.ro")
    load("data/pipeline/v17/bem/PANCAN_methyl.ro")
    load("data/pipeline/v17/bem/PANCAN_speed.ro")
    commonCell <- intersect(rownames(PANCAN_mut), intersect(rownames(PANCAN_methyl), rownames(PANCAN_speed)))
    feature <- cbind(PANCAN_mut[commonCell,], PANCAN_methyl[commonCell,], PANCAN_speed[commonCell,])
  }
  if (feat=="v17_PANCAN_mut.methyl.tissue"){
    load("data/pipeline/v17/bem/PANCAN_mut.ro")
    load("data/pipeline/v17/bem/PANCAN_methyl.ro")
    load("data/pipeline/v17/tissue.ro")
    commonCell <- intersect(rownames(PANCAN_mut), intersect(rownames(PANCAN_methyl), rownames(tissue)))
    feature <- cbind(PANCAN_mut[commonCell,], PANCAN_methyl[commonCell,], tissue[commonCell,])
  }
  if (feat=="v17_PANCAN_mut.gex.cnv"){
    load("data/pipeline/v17/bem/PANCAN_mut.ro")
    load("data/pipeline/v17/bem/PANCAN_gex.ro")
    load("data/pipeline/v17/bem/PANCAN_cnv.ro")
    commonCell <- intersect(rownames(PANCAN_mut), intersect(rownames(PANCAN_gex), rownames(PANCAN_cnv)))
    feature <- cbind(PANCAN_mut[commonCell,], PANCAN_gex[commonCell,], PANCAN_cnv[commonCell,])
  }
  if (feat=="v17_PANCAN_mut.gex.speed"){
    load("data/pipeline/v17/bem/PANCAN_mut.ro")
    load("data/pipeline/v17/bem/PANCAN_gex.ro")
    load("data/pipeline/v17/bem/PANCAN_speed.ro")
    commonCell <- intersect(rownames(PANCAN_mut), intersect(rownames(PANCAN_gex), rownames(PANCAN_speed)))
    feature <- cbind(PANCAN_mut[commonCell,], PANCAN_gex[commonCell,], PANCAN_speed[commonCell,])
  }
  if (feat=="v17_PANCAN_mut.gex.tissue"){
    load("data/pipeline/v17/bem/PANCAN_mut.ro")
    load("data/pipeline/v17/bem/PANCAN_gex.ro")
    load("data/pipeline/v17/tissue.ro")
    commonCell <- intersect(rownames(PANCAN_mut), intersect(rownames(PANCAN_gex), rownames(tissue)))
    feature <- cbind(PANCAN_mut[commonCell,], PANCAN_gex[commonCell,], tissue[commonCell,])
  }
  if (feat=="v17_PANCAN_mut.cnv.speed"){
    load("data/pipeline/v17/bem/PANCAN_mut.ro")
    load("data/pipeline/v17/bem/PANCAN_cnv.ro")
    load("data/pipeline/v17/bem/PANCAN_speed.ro")
    commonCell <- intersect(rownames(PANCAN_mut), intersect(rownames(PANCAN_cnv), rownames(PANCAN_speed)))
    feature <- cbind(PANCAN_mut[commonCell,], PANCAN_cnv[commonCell,], PANCAN_speed[commonCell,])
  }
  if (feat=="v17_PANCAN_mut.cnv.tissue"){
    load("data/pipeline/v17/bem/PANCAN_mut.ro")
    load("data/pipeline/v17/bem/PANCAN_cnv.ro")
    load("data/pipeline/v17/tissue.ro")
    commonCell <- intersect(rownames(PANCAN_mut), intersect(rownames(PANCAN_cnv), rownames(tissue)))
    feature <- cbind(PANCAN_mut[commonCell,], PANCAN_cnv[commonCell,], tissue[commonCell,])
  }
  if (feat=="v17_PANCAN_mut.speed.tissue"){
    load("data/pipeline/v17/bem/PANCAN_mut.ro")
    load("data/pipeline/v17/bem/PANCAN_speed.ro")
    load("data/pipeline/v17/tissue.ro")
    commonCell <- intersect(rownames(PANCAN_mut), intersect(rownames(PANCAN_speed), rownames(tissue)))
    feature <- cbind(PANCAN_mut[commonCell,], PANCAN_speed[commonCell,], tissue[commonCell,])
  }
  if (feat=="v17_PANCAN_methyl.gex.cnv"){
    load("data/pipeline/v17/bem/PANCAN_methyl.ro")
    load("data/pipeline/v17/bem/PANCAN_gex.ro")
    load("data/pipeline/v17/bem/PANCAN_cnv.ro")
    commonCell <- intersect(rownames(PANCAN_methyl), intersect(rownames(PANCAN_gex), rownames(PANCAN_cnv)))
    feature <- cbind(PANCAN_methyl[commonCell,], PANCAN_gex[commonCell,], PANCAN_cnv[commonCell,])
  }
  if (feat=="v17_PANCAN_methyl.gex.speed"){
    load("data/pipeline/v17/bem/PANCAN_methyl.ro")
    load("data/pipeline/v17/bem/PANCAN_gex.ro")
    load("data/pipeline/v17/bem/PANCAN_speed.ro")
    commonCell <- intersect(rownames(PANCAN_methyl), intersect(rownames(PANCAN_gex), rownames(PANCAN_speed)))
    feature <- cbind(PANCAN_methyl[commonCell,], PANCAN_gex[commonCell,], PANCAN_speed[commonCell,])
  }
  if (feat=="v17_PANCAN_methyl.gex.tissue"){
    load("data/pipeline/v17/bem/PANCAN_methyl.ro")
    load("data/pipeline/v17/bem/PANCAN_gex.ro")
    load("data/pipeline/v17/tissue.ro")
    commonCell <- intersect(rownames(PANCAN_methyl), intersect(rownames(PANCAN_gex), rownames(tissue)))
    feature <- cbind(PANCAN_methyl[commonCell,], PANCAN_gex[commonCell,], tissue[commonCell,])
  }
  if (feat=="v17_PANCAN_methyl.cnv.speed"){
    load("data/pipeline/v17/bem/PANCAN_methyl.ro")
    load("data/pipeline/v17/bem/PANCAN_cnv.ro")
    load("data/pipeline/v17/bem/PANCAN_speed.ro")
    commonCell <- intersect(rownames(PANCAN_methyl), intersect(rownames(PANCAN_cnv), rownames(PANCAN_speed)))
    feature <- cbind(PANCAN_methyl[commonCell,], PANCAN_cnv[commonCell,], PANCAN_speed[commonCell,])
  }
  if (feat=="v17_PANCAN_methyl.cnv.tissue"){
    load("data/pipeline/v17/bem/PANCAN_methyl.ro")
    load("data/pipeline/v17/bem/PANCAN_cnv.ro")
    load("data/pipeline/v17/tissue.ro")
    commonCell <- intersect(rownames(PANCAN_methyl), intersect(rownames(PANCAN_cnv), rownames(tissue)))
    feature <- cbind(PANCAN_methyl[commonCell,], PANCAN_cnv[commonCell,], tissue[commonCell,])
  }
  if (feat=="v17_PANCAN_methyl.speed.tissue"){
    load("data/pipeline/v17/bem/PANCAN_methyl.ro")
    load("data/pipeline/v17/bem/PANCAN_speed.ro")
    load("data/pipeline/v17/tissue.ro")
    commonCell <- intersect(rownames(PANCAN_methyl), intersect(rownames(PANCAN_speed), rownames(tissue)))
    feature <- cbind(PANCAN_methyl[commonCell,], PANCAN_speed[commonCell,], tissue[commonCell,])
  }
  if (feat=="v17_PANCAN_gex.cnv.speed"){
    load("data/pipeline/v17/bem/PANCAN_gex.ro")
    load("data/pipeline/v17/bem/PANCAN_cnv.ro")
    load("data/pipeline/v17/bem/PANCAN_speed.ro")
    commonCell <- intersect(rownames(PANCAN_gex), intersect(rownames(PANCAN_cnv), rownames(PANCAN_speed)))
    feature <- cbind(PANCAN_gex[commonCell,], PANCAN_cnv[commonCell,], PANCAN_speed[commonCell,])
  }
  if (feat=="v17_PANCAN_gex.cnv.tissue"){
    load("data/pipeline/v17/bem/PANCAN_gex.ro")
    load("data/pipeline/v17/bem/PANCAN_cnv.ro")
    load("data/pipeline/v17/tissue.ro")
    commonCell <- intersect(rownames(PANCAN_gex), intersect(rownames(PANCAN_cnv), rownames(tissue)))
    feature <- cbind(PANCAN_gex[commonCell,], PANCAN_cnv[commonCell,], tissue[commonCell,])
  }
  if (feat=="v17_PANCAN_gex.speed.tissue"){
    load("data/pipeline/v17/bem/PANCAN_gex.ro")
    load("data/pipeline/v17/bem/PANCAN_speed.ro")
    load("data/pipeline/v17/tissue.ro")
    commonCell <- intersect(rownames(PANCAN_gex), intersect(rownames(PANCAN_speed), rownames(tissue)))
    feature <- cbind(PANCAN_gex[commonCell,], PANCAN_speed[commonCell,], tissue[commonCell,])
  }
  if (feat=="v17_PANCAN_cnv.speed.tissue"){
    load("data/pipeline/v17/bem/PANCAN_cnv.ro")
    load("data/pipeline/v17/bem/PANCAN_speed.ro")
    load("data/pipeline/v17/tissue.ro")
    commonCell <- intersect(rownames(PANCAN_cnv), intersect(rownames(PANCAN_speed), rownames(tissue)))
    feature <- cbind(PANCAN_cnv[commonCell,], PANCAN_speed[commonCell,], tissue[commonCell,])
  }
  if (feat=="v17_PANCAN_mut.methyl.gex.cnv"){
    load("data/pipeline/v17/bem/PANCAN_mut.ro")
    load("data/pipeline/v17/bem/PANCAN_methyl.ro")
    load("data/pipeline/v17/bem/PANCAN_gex.ro")
    load("data/pipeline/v17/bem/PANCAN_cnv.ro")
    commonCell <- intersect(rownames(PANCAN_mut), intersect(rownames(PANCAN_methyl), intersect(rownames(PANCAN_gex), rownames(PANCAN_cnv))))
    feature <- cbind(PANCAN_mut[commonCell,], PANCAN_methyl[commonCell,], PANCAN_gex[commonCell,], PANCAN_cnv[commonCell,])
  }
  if (feat=="v17_PANCAN_mut.methyl.gex.speed"){
    load("data/pipeline/v17/bem/PANCAN_mut.ro")
    load("data/pipeline/v17/bem/PANCAN_methyl.ro")
    load("data/pipeline/v17/bem/PANCAN_gex.ro")
    load("data/pipeline/v17/bem/PANCAN_speed.ro")
    commonCell <- intersect(rownames(PANCAN_mut), intersect(rownames(PANCAN_methyl), intersect(rownames(PANCAN_gex), rownames(PANCAN_speed))))
    feature <- cbind(PANCAN_mut[commonCell,], PANCAN_methyl[commonCell,], PANCAN_gex[commonCell,], PANCAN_speed[commonCell,])
  }
  if (feat=="v17_PANCAN_mut.methyl.gex.tissue"){
    load("data/pipeline/v17/bem/PANCAN_mut.ro")
    load("data/pipeline/v17/bem/PANCAN_methyl.ro")
    load("data/pipeline/v17/bem/PANCAN_gex.ro")
    load("data/pipeline/v17/tissue.ro")
    commonCell <- intersect(rownames(PANCAN_mut), intersect(rownames(PANCAN_methyl), intersect(rownames(PANCAN_gex), rownames(tissue))))
    feature <- cbind(PANCAN_mut[commonCell,], PANCAN_methyl[commonCell,], PANCAN_gex[commonCell,], tissue[commonCell,])
  }
  if (feat=="v17_PANCAN_mut.methyl.cnv.speed"){
    load("data/pipeline/v17/bem/PANCAN_mut.ro")
    load("data/pipeline/v17/bem/PANCAN_methyl.ro")
    load("data/pipeline/v17/bem/PANCAN_cnv.ro")
    load("data/pipeline/v17/bem/PANCAN_speed.ro")
    commonCell <- intersect(rownames(PANCAN_mut), intersect(rownames(PANCAN_methyl), intersect(rownames(PANCAN_cnv), rownames(PANCAN_speed))))
    feature <- cbind(PANCAN_mut[commonCell,], PANCAN_methyl[commonCell,], PANCAN_cnv[commonCell,], PANCAN_speed[commonCell,])
  }
  if (feat=="v17_PANCAN_mut.methyl.cnv.tissue"){
    load("data/pipeline/v17/bem/PANCAN_mut.ro")
    load("data/pipeline/v17/bem/PANCAN_methyl.ro")
    load("data/pipeline/v17/bem/PANCAN_cnv.ro")
    load("data/pipeline/v17/tissue.ro")
    commonCell <- intersect(rownames(PANCAN_mut), intersect(rownames(PANCAN_methyl), intersect(rownames(PANCAN_cnv), rownames(tissue))))
    feature <- cbind(PANCAN_mut[commonCell,], PANCAN_methyl[commonCell,], PANCAN_cnv[commonCell,], tissue[commonCell,])
  }
  if (feat=="v17_PANCAN_mut.methyl.speed.tissue"){
    load("data/pipeline/v17/bem/PANCAN_mut.ro")
    load("data/pipeline/v17/bem/PANCAN_methyl.ro")
    load("data/pipeline/v17/bem/PANCAN_speed.ro")
    load("data/pipeline/v17/tissue.ro")
    commonCell <- intersect(rownames(PANCAN_mut), intersect(rownames(PANCAN_methyl), intersect(rownames(PANCAN_speed), rownames(tissue))))
    feature <- cbind(PANCAN_mut[commonCell,], PANCAN_methyl[commonCell,], PANCAN_speed[commonCell,], tissue[commonCell,])
  }
  if (feat=="v17_PANCAN_mut.gex.cnv.speed"){
    load("data/pipeline/v17/bem/PANCAN_mut.ro")
    load("data/pipeline/v17/bem/PANCAN_gex.ro")
    load("data/pipeline/v17/bem/PANCAN_cnv.ro")
    load("data/pipeline/v17/bem/PANCAN_speed.ro")
    commonCell <- intersect(rownames(PANCAN_mut), intersect(rownames(PANCAN_gex), intersect(rownames(PANCAN_cnv), rownames(PANCAN_speed))))
    feature <- cbind(PANCAN_mut[commonCell,], PANCAN_gex[commonCell,], PANCAN_cnv[commonCell,], PANCAN_speed[commonCell,])
  }
  if (feat=="v17_PANCAN_mut.gex.cnv.tissue"){
    load("data/pipeline/v17/bem/PANCAN_mut.ro")
    load("data/pipeline/v17/bem/PANCAN_gex.ro")
    load("data/pipeline/v17/bem/PANCAN_cnv.ro")
    load("data/pipeline/v17/tissue.ro")
    commonCell <- intersect(rownames(PANCAN_mut), intersect(rownames(PANCAN_gex), intersect(rownames(PANCAN_cnv), rownames(tissue))))
    feature <- cbind(PANCAN_mut[commonCell,], PANCAN_gex[commonCell,], PANCAN_cnv[commonCell,], tissue[commonCell,])
  }
  if (feat=="v17_PANCAN_mut.gex.speed.tissue"){
    load("data/pipeline/v17/bem/PANCAN_mut.ro")
    load("data/pipeline/v17/bem/PANCAN_gex.ro")
    load("data/pipeline/v17/bem/PANCAN_speed.ro")
    load("data/pipeline/v17/tissue.ro")
    commonCell <- intersect(rownames(PANCAN_mut), intersect(rownames(PANCAN_gex), intersect(rownames(PANCAN_speed), rownames(tissue))))
    feature <- cbind(PANCAN_mut[commonCell,], PANCAN_gex[commonCell,], PANCAN_speed[commonCell,], tissue[commonCell,])
  }
  if (feat=="v17_PANCAN_mut.cnv.speed.tissue"){
    load("data/pipeline/v17/bem/PANCAN_mut.ro")
    load("data/pipeline/v17/bem/PANCAN_cnv.ro")
    load("data/pipeline/v17/bem/PANCAN_speed.ro")
    load("data/pipeline/v17/tissue.ro")
    commonCell <- intersect(rownames(PANCAN_mut), intersect(rownames(PANCAN_cnv), intersect(rownames(PANCAN_speed), rownames(tissue))))
    feature <- cbind(PANCAN_mut[commonCell,], PANCAN_cnv[commonCell,], PANCAN_speed[commonCell,], tissue[commonCell,])
  }
  if (feat=="v17_PANCAN_methyl.gex.cnv.speed"){
    load("data/pipeline/v17/bem/PANCAN_methyl.ro")
    load("data/pipeline/v17/bem/PANCAN_gex.ro")
    load("data/pipeline/v17/bem/PANCAN_cnv.ro")
    load("data/pipeline/v17/bem/PANCAN_speed.ro")
    commonCell <- intersect(rownames(PANCAN_methyl), intersect(rownames(PANCAN_gex), intersect(rownames(PANCAN_cnv), rownames(PANCAN_speed))))
    feature <- cbind(PANCAN_methyl[commonCell,], PANCAN_gex[commonCell,], PANCAN_cnv[commonCell,], PANCAN_speed[commonCell,])
  }
  if (feat=="v17_PANCAN_methyl.gex.cnv.tissue"){
    load("data/pipeline/v17/bem/PANCAN_methyl.ro")
    load("data/pipeline/v17/bem/PANCAN_gex.ro")
    load("data/pipeline/v17/bem/PANCAN_cnv.ro")
    load("data/pipeline/v17/tissue.ro")
    commonCell <- intersect(rownames(PANCAN_methyl), intersect(rownames(PANCAN_gex), intersect(rownames(PANCAN_cnv), rownames(tissue))))
    feature <- cbind(PANCAN_methyl[commonCell,], PANCAN_gex[commonCell,], PANCAN_cnv[commonCell,], tissue[commonCell,])
  }
  if (feat=="v17_PANCAN_methyl.gex.speed.tissue"){
    load("data/pipeline/v17/bem/PANCAN_methyl.ro")
    load("data/pipeline/v17/bem/PANCAN_gex.ro")
    load("data/pipeline/v17/bem/PANCAN_speed.ro")
    load("data/pipeline/v17/tissue.ro")
    commonCell <- intersect(rownames(PANCAN_methyl), intersect(rownames(PANCAN_gex), intersect(rownames(PANCAN_speed), rownames(tissue))))
    feature <- cbind(PANCAN_methyl[commonCell,], PANCAN_gex[commonCell,], PANCAN_speed[commonCell,], tissue[commonCell,])
  }
  if (feat=="v17_PANCAN_methyl.cnv.speed.tissue"){
    load("data/pipeline/v17/bem/PANCAN_methyl.ro")
    load("data/pipeline/v17/bem/PANCAN_cnv.ro")
    load("data/pipeline/v17/bem/PANCAN_speed.ro")
    load("data/pipeline/v17/tissue.ro")
    commonCell <- intersect(rownames(PANCAN_methyl), intersect(rownames(PANCAN_cnv), intersect(rownames(PANCAN_speed), rownames(tissue))))
    feature <- cbind(PANCAN_methyl[commonCell,], PANCAN_cnv[commonCell,], PANCAN_speed[commonCell,], tissue[commonCell,])
  }
  if (feat=="v17_PANCAN_gex.cnv.speed.tissue"){
    load("data/pipeline/v17/bem/PANCAN_gex.ro")
    load("data/pipeline/v17/bem/PANCAN_cnv.ro")
    load("data/pipeline/v17/bem/PANCAN_speed.ro")
    load("data/pipeline/v17/tissue.ro")
    commonCell <- intersect(rownames(PANCAN_gex), intersect(rownames(PANCAN_cnv), intersect(rownames(PANCAN_speed), rownames(tissue))))
    feature <- cbind(PANCAN_gex[commonCell,], PANCAN_cnv[commonCell,], PANCAN_speed[commonCell,], tissue[commonCell,])
  }
  if (feat=="v17_PANCAN_mut.methyl.gex.cnv.speed"){
    load("data/pipeline/v17/bem/PANCAN_mut.ro")
    load("data/pipeline/v17/bem/PANCAN_methyl.ro")
    load("data/pipeline/v17/bem/PANCAN_gex.ro")
    load("data/pipeline/v17/bem/PANCAN_cnv.ro")
    load("data/pipeline/v17/bem/PANCAN_speed.ro")
    commonCell <- intersect(rownames(PANCAN_mut), intersect(rownames(PANCAN_methyl), intersect(rownames(PANCAN_gex), intersect(rownames(PANCAN_cnv), rownames(PANCAN_speed)))))
    feature <- cbind(PANCAN_mut[commonCell,], PANCAN_methyl[commonCell,], PANCAN_gex[commonCell,], PANCAN_cnv[commonCell,], PANCAN_speed[commonCell,])
  }
  if (feat=="v17_PANCAN_mut.methyl.gex.cnv.tissue"){
    load("data/pipeline/v17/bem/PANCAN_mut.ro")
    load("data/pipeline/v17/bem/PANCAN_methyl.ro")
    load("data/pipeline/v17/bem/PANCAN_gex.ro")
    load("data/pipeline/v17/bem/PANCAN_cnv.ro")
    load("data/pipeline/v17/tissue.ro")
    commonCell <- intersect(rownames(PANCAN_mut), intersect(rownames(PANCAN_methyl), intersect(rownames(PANCAN_gex), intersect(rownames(PANCAN_cnv), rownames(tissue)))))
    feature <- cbind(PANCAN_mut[commonCell,], PANCAN_methyl[commonCell,], PANCAN_gex[commonCell,], PANCAN_cnv[commonCell,], tissue[commonCell,])
  }
  if (feat=="v17_PANCAN_mut.methyl.gex.speed.tissue"){
    load("data/pipeline/v17/bem/PANCAN_mut.ro")
    load("data/pipeline/v17/bem/PANCAN_methyl.ro")
    load("data/pipeline/v17/bem/PANCAN_gex.ro")
    load("data/pipeline/v17/bem/PANCAN_speed.ro")
    load("data/pipeline/v17/tissue.ro")
    commonCell <- intersect(rownames(PANCAN_mut), intersect(rownames(PANCAN_methyl), intersect(rownames(PANCAN_gex), intersect(rownames(PANCAN_speed), rownames(tissue)))))
    feature <- cbind(PANCAN_mut[commonCell,], PANCAN_methyl[commonCell,], PANCAN_gex[commonCell,], PANCAN_speed[commonCell,], tissue[commonCell,])
  }
  if (feat=="v17_PANCAN_mut.methyl.cnv.speed.tissue"){
    load("data/pipeline/v17/bem/PANCAN_mut.ro")
    load("data/pipeline/v17/bem/PANCAN_methyl.ro")
    load("data/pipeline/v17/bem/PANCAN_cnv.ro")
    load("data/pipeline/v17/bem/PANCAN_speed.ro")
    load("data/pipeline/v17/tissue.ro")
    commonCell <- intersect(rownames(PANCAN_mut), intersect(rownames(PANCAN_methyl), intersect(rownames(PANCAN_cnv), intersect(rownames(PANCAN_speed), rownames(tissue)))))
    feature <- cbind(PANCAN_mut[commonCell,], PANCAN_methyl[commonCell,], PANCAN_cnv[commonCell,], PANCAN_speed[commonCell,], tissue[commonCell,])
  }
  if (feat=="v17_PANCAN_mut.gex.cnv.speed.tissue"){
    load("data/pipeline/v17/bem/PANCAN_mut.ro")
    load("data/pipeline/v17/bem/PANCAN_gex.ro")
    load("data/pipeline/v17/bem/PANCAN_cnv.ro")
    load("data/pipeline/v17/bem/PANCAN_speed.ro")
    load("data/pipeline/v17/tissue.ro")
    commonCell <- intersect(rownames(PANCAN_mut), intersect(rownames(PANCAN_gex), intersect(rownames(PANCAN_cnv), intersect(rownames(PANCAN_speed), rownames(tissue)))))
    feature <- cbind(PANCAN_mut[commonCell,], PANCAN_gex[commonCell,], PANCAN_cnv[commonCell,], PANCAN_speed[commonCell,], tissue[commonCell,])
  }
  if (feat=="v17_PANCAN_methyl.gex.cnv.speed.tissue"){
    load("data/pipeline/v17/bem/PANCAN_methyl.ro")
    load("data/pipeline/v17/bem/PANCAN_gex.ro")
    load("data/pipeline/v17/bem/PANCAN_cnv.ro")
    load("data/pipeline/v17/bem/PANCAN_speed.ro")
    load("data/pipeline/v17/tissue.ro")
    commonCell <- intersect(rownames(PANCAN_methyl), intersect(rownames(PANCAN_gex), intersect(rownames(PANCAN_cnv), intersect(rownames(PANCAN_speed), rownames(tissue)))))
    feature <- cbind(PANCAN_methyl[commonCell,], PANCAN_gex[commonCell,], PANCAN_cnv[commonCell,], PANCAN_speed[commonCell,], tissue[commonCell,])
  }
  if (feat=="v17_PANCAN_mut.methyl.gex.cnv.speed.tissue"){
    load("data/pipeline/v17/bem/PANCAN_mut.ro")
    load("data/pipeline/v17/bem/PANCAN_methyl.ro")
    load("data/pipeline/v17/bem/PANCAN_gex.ro")
    load("data/pipeline/v17/bem/PANCAN_cnv.ro")
    load("data/pipeline/v17/bem/PANCAN_speed.ro")
    load("data/pipeline/v17/tissue.ro")
    commonCell <- intersect(rownames(PANCAN_mut), intersect(rownames(PANCAN_methyl), intersect(rownames(PANCAN_gex), intersect(rownames(PANCAN_cnv), intersect(rownames(PANCAN_speed), rownames(tissue))))))
    feature <- cbind(PANCAN_mut[commonCell,], PANCAN_methyl[commonCell,], PANCAN_gex[commonCell,], PANCAN_cnv[commonCell,], PANCAN_speed[commonCell,], tissue[commonCell,])
  }

  #############################################################################################
  # load NCI60_v2_01_combinations
  #############################################################################################
  if (feat=="NCI60_TISSUE"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    commonCell <- rownames(TISSUE)
    feature <- cbind(TISSUE[commonCell,])
  }
  if (feat=="NCI60_RPPA"){
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    commonCell <- rownames(RPPA)
    feature <- cbind(RPPA[commonCell,])
  }
  if (feat=="NCI60_mRNA"){
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    commonCell <- rownames(mRNA)
    feature <- cbind(mRNA[commonCell,])
  }
  if (feat=="NCI60_SHOTGUN"){
    load("data/pipeline/NCI60_final_v2/SHOTGUN.ro")
    commonCell <- rownames(SHOTGUN)
    feature <- cbind(SHOTGUN[commonCell,])
  }
  if (feat=="NCI60_genomics"){
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    commonCell <- rownames(genomics)
    feature <- cbind(genomics[commonCell,])
  }
  if (feat=="NCI60_SWATH"){
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    commonCell <- rownames(SWATH)
    feature <- cbind(SWATH[commonCell,])
  }
  if (feat=="NCI60_SHOTGUN_500"){
    load("data/pipeline/NCI60_final_v2/SHOTGUN_500.ro")
    commonCell <- rownames(SHOTGUN_500)
    feature <- cbind(SHOTGUN_500[commonCell,])
  }
  if (feat=="NCI60_TISSUE.RPPA"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    commonCell <- intersect(rownames(TISSUE), rownames(RPPA))
    feature <- cbind(TISSUE[commonCell,], RPPA[commonCell,])
  }
  if (feat=="NCI60_TISSUE.mRNA"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    commonCell <- intersect(rownames(TISSUE), rownames(mRNA))
    feature <- cbind(TISSUE[commonCell,], mRNA[commonCell,])
  }
  if (feat=="NCI60_TISSUE.SHOTGUN"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN.ro")
    commonCell <- intersect(rownames(TISSUE), rownames(SHOTGUN))
    feature <- cbind(TISSUE[commonCell,], SHOTGUN[commonCell,])
  }
  if (feat=="NCI60_TISSUE.genomics"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    commonCell <- intersect(rownames(TISSUE), rownames(genomics))
    feature <- cbind(TISSUE[commonCell,], genomics[commonCell,])
  }
  if (feat=="NCI60_TISSUE.SWATH"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    commonCell <- intersect(rownames(TISSUE), rownames(SWATH))
    feature <- cbind(TISSUE[commonCell,], SWATH[commonCell,])
  }
  if (feat=="NCI60_TISSUE.SHOTGUN_500"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN_500.ro")
    commonCell <- intersect(rownames(TISSUE), rownames(SHOTGUN_500))
    feature <- cbind(TISSUE[commonCell,], SHOTGUN_500[commonCell,])
  }
  if (feat=="NCI60_RPPA.mRNA"){
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    commonCell <- intersect(rownames(RPPA), rownames(mRNA))
    feature <- cbind(RPPA[commonCell,], mRNA[commonCell,])
  }
  if (feat=="NCI60_RPPA.SHOTGUN"){
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN.ro")
    commonCell <- intersect(rownames(RPPA), rownames(SHOTGUN))
    feature <- cbind(RPPA[commonCell,], SHOTGUN[commonCell,])
  }
  if (feat=="NCI60_RPPA.genomics"){
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    commonCell <- intersect(rownames(RPPA), rownames(genomics))
    feature <- cbind(RPPA[commonCell,], genomics[commonCell,])
  }
  if (feat=="NCI60_RPPA.SWATH"){
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    commonCell <- intersect(rownames(RPPA), rownames(SWATH))
    feature <- cbind(RPPA[commonCell,], SWATH[commonCell,])
  }
  if (feat=="NCI60_RPPA.SHOTGUN_500"){
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN_500.ro")
    commonCell <- intersect(rownames(RPPA), rownames(SHOTGUN_500))
    feature <- cbind(RPPA[commonCell,], SHOTGUN_500[commonCell,])
  }
  if (feat=="NCI60_mRNA.SHOTGUN"){
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN.ro")
    commonCell <- intersect(rownames(mRNA), rownames(SHOTGUN))
    feature <- cbind(mRNA[commonCell,], SHOTGUN[commonCell,])
  }
  if (feat=="NCI60_mRNA.genomics"){
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    commonCell <- intersect(rownames(mRNA), rownames(genomics))
    feature <- cbind(mRNA[commonCell,], genomics[commonCell,])
  }
  if (feat=="NCI60_mRNA.SWATH"){
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    commonCell <- intersect(rownames(mRNA), rownames(SWATH))
    feature <- cbind(mRNA[commonCell,], SWATH[commonCell,])
  }
  if (feat=="NCI60_mRNA.SHOTGUN_500"){
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN_500.ro")
    commonCell <- intersect(rownames(mRNA), rownames(SHOTGUN_500))
    feature <- cbind(mRNA[commonCell,], SHOTGUN_500[commonCell,])
  }
  if (feat=="NCI60_SHOTGUN.genomics"){
    load("data/pipeline/NCI60_final_v2/SHOTGUN.ro")
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    commonCell <- intersect(rownames(SHOTGUN), rownames(genomics))
    feature <- cbind(SHOTGUN[commonCell,], genomics[commonCell,])
  }
  if (feat=="NCI60_SHOTGUN.SWATH"){
    load("data/pipeline/NCI60_final_v2/SHOTGUN.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    commonCell <- intersect(rownames(SHOTGUN), rownames(SWATH))
    feature <- cbind(SHOTGUN[commonCell,], SWATH[commonCell,])
  }
  if (feat=="NCI60_genomics.SWATH"){
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    commonCell <- intersect(rownames(genomics), rownames(SWATH))
    feature <- cbind(genomics[commonCell,], SWATH[commonCell,])
  }
  if (feat=="NCI60_genomics.SHOTGUN_500"){
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN_500.ro")
    commonCell <- intersect(rownames(genomics), rownames(SHOTGUN_500))
    feature <- cbind(genomics[commonCell,], SHOTGUN_500[commonCell,])
  }
  if (feat=="NCI60_SWATH.SHOTGUN_500"){
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN_500.ro")
    commonCell <- intersect(rownames(SWATH), rownames(SHOTGUN_500))
    feature <- cbind(SWATH[commonCell,], SHOTGUN_500[commonCell,])
  }
  if (feat=="NCI60_TISSUE.RPPA.mRNA"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    commonCell <- intersect(rownames(TISSUE), intersect(rownames(RPPA), rownames(mRNA)))
    feature <- cbind(TISSUE[commonCell,], RPPA[commonCell,], mRNA[commonCell,])
  }
  if (feat=="NCI60_TISSUE.RPPA.SHOTGUN"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN.ro")
    commonCell <- intersect(rownames(TISSUE), intersect(rownames(RPPA), rownames(SHOTGUN)))
    feature <- cbind(TISSUE[commonCell,], RPPA[commonCell,], SHOTGUN[commonCell,])
  }
  if (feat=="NCI60_TISSUE.RPPA.genomics"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    commonCell <- intersect(rownames(TISSUE), intersect(rownames(RPPA), rownames(genomics)))
    feature <- cbind(TISSUE[commonCell,], RPPA[commonCell,], genomics[commonCell,])
  }
  if (feat=="NCI60_TISSUE.RPPA.SWATH"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    commonCell <- intersect(rownames(TISSUE), intersect(rownames(RPPA), rownames(SWATH)))
    feature <- cbind(TISSUE[commonCell,], RPPA[commonCell,], SWATH[commonCell,])
  }
  if (feat=="NCI60_TISSUE.RPPA.SHOTGUN_500"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN_500.ro")
    commonCell <- intersect(rownames(TISSUE), intersect(rownames(RPPA), rownames(SHOTGUN_500)))
    feature <- cbind(TISSUE[commonCell,], RPPA[commonCell,], SHOTGUN_500[commonCell,])
  }
  if (feat=="NCI60_TISSUE.mRNA.SHOTGUN"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN.ro")
    commonCell <- intersect(rownames(TISSUE), intersect(rownames(mRNA), rownames(SHOTGUN)))
    feature <- cbind(TISSUE[commonCell,], mRNA[commonCell,], SHOTGUN[commonCell,])
  }
  if (feat=="NCI60_TISSUE.mRNA.genomics"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    commonCell <- intersect(rownames(TISSUE), intersect(rownames(mRNA), rownames(genomics)))
    feature <- cbind(TISSUE[commonCell,], mRNA[commonCell,], genomics[commonCell,])
  }
  if (feat=="NCI60_TISSUE.mRNA.SWATH"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    commonCell <- intersect(rownames(TISSUE), intersect(rownames(mRNA), rownames(SWATH)))
    feature <- cbind(TISSUE[commonCell,], mRNA[commonCell,], SWATH[commonCell,])
  }
  if (feat=="NCI60_TISSUE.mRNA.SHOTGUN_500"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN_500.ro")
    commonCell <- intersect(rownames(TISSUE), intersect(rownames(mRNA), rownames(SHOTGUN_500)))
    feature <- cbind(TISSUE[commonCell,], mRNA[commonCell,], SHOTGUN_500[commonCell,])
  }
  if (feat=="NCI60_TISSUE.SHOTGUN.genomics"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN.ro")
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    commonCell <- intersect(rownames(TISSUE), intersect(rownames(SHOTGUN), rownames(genomics)))
    feature <- cbind(TISSUE[commonCell,], SHOTGUN[commonCell,], genomics[commonCell,])
  }
  if (feat=="NCI60_TISSUE.SHOTGUN.SWATH"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    commonCell <- intersect(rownames(TISSUE), intersect(rownames(SHOTGUN), rownames(SWATH)))
    feature <- cbind(TISSUE[commonCell,], SHOTGUN[commonCell,], SWATH[commonCell,])
  }
  if (feat=="NCI60_TISSUE.genomics.SWATH"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    commonCell <- intersect(rownames(TISSUE), intersect(rownames(genomics), rownames(SWATH)))
    feature <- cbind(TISSUE[commonCell,], genomics[commonCell,], SWATH[commonCell,])
  }
  if (feat=="NCI60_TISSUE.genomics.SHOTGUN_500"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN_500.ro")
    commonCell <- intersect(rownames(TISSUE), intersect(rownames(genomics), rownames(SHOTGUN_500)))
    feature <- cbind(TISSUE[commonCell,], genomics[commonCell,], SHOTGUN_500[commonCell,])
  }
  if (feat=="NCI60_TISSUE.SWATH.SHOTGUN_500"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN_500.ro")
    commonCell <- intersect(rownames(TISSUE), intersect(rownames(SWATH), rownames(SHOTGUN_500)))
    feature <- cbind(TISSUE[commonCell,], SWATH[commonCell,], SHOTGUN_500[commonCell,])
  }
  if (feat=="NCI60_RPPA.mRNA.SHOTGUN"){
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN.ro")
    commonCell <- intersect(rownames(RPPA), intersect(rownames(mRNA), rownames(SHOTGUN)))
    feature <- cbind(RPPA[commonCell,], mRNA[commonCell,], SHOTGUN[commonCell,])
  }
  if (feat=="NCI60_RPPA.mRNA.genomics"){
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    commonCell <- intersect(rownames(RPPA), intersect(rownames(mRNA), rownames(genomics)))
    feature <- cbind(RPPA[commonCell,], mRNA[commonCell,], genomics[commonCell,])
  }
  if (feat=="NCI60_RPPA.mRNA.SWATH"){
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    commonCell <- intersect(rownames(RPPA), intersect(rownames(mRNA), rownames(SWATH)))
    feature <- cbind(RPPA[commonCell,], mRNA[commonCell,], SWATH[commonCell,])
  }
  if (feat=="NCI60_RPPA.mRNA.SHOTGUN_500"){
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN_500.ro")
    commonCell <- intersect(rownames(RPPA), intersect(rownames(mRNA), rownames(SHOTGUN_500)))
    feature <- cbind(RPPA[commonCell,], mRNA[commonCell,], SHOTGUN_500[commonCell,])
  }
  if (feat=="NCI60_RPPA.SHOTGUN.genomics"){
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN.ro")
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    commonCell <- intersect(rownames(RPPA), intersect(rownames(SHOTGUN), rownames(genomics)))
    feature <- cbind(RPPA[commonCell,], SHOTGUN[commonCell,], genomics[commonCell,])
  }
  if (feat=="NCI60_RPPA.SHOTGUN.SWATH"){
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    commonCell <- intersect(rownames(RPPA), intersect(rownames(SHOTGUN), rownames(SWATH)))
    feature <- cbind(RPPA[commonCell,], SHOTGUN[commonCell,], SWATH[commonCell,])
  }
  if (feat=="NCI60_RPPA.genomics.SWATH"){
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    commonCell <- intersect(rownames(RPPA), intersect(rownames(genomics), rownames(SWATH)))
    feature <- cbind(RPPA[commonCell,], genomics[commonCell,], SWATH[commonCell,])
  }
  if (feat=="NCI60_RPPA.genomics.SHOTGUN_500"){
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN_500.ro")
    commonCell <- intersect(rownames(RPPA), intersect(rownames(genomics), rownames(SHOTGUN_500)))
    feature <- cbind(RPPA[commonCell,], genomics[commonCell,], SHOTGUN_500[commonCell,])
  }
  if (feat=="NCI60_RPPA.SWATH.SHOTGUN_500"){
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN_500.ro")
    commonCell <- intersect(rownames(RPPA), intersect(rownames(SWATH), rownames(SHOTGUN_500)))
    feature <- cbind(RPPA[commonCell,], SWATH[commonCell,], SHOTGUN_500[commonCell,])
  }
  if (feat=="NCI60_mRNA.SHOTGUN.genomics"){
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN.ro")
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    commonCell <- intersect(rownames(mRNA), intersect(rownames(SHOTGUN), rownames(genomics)))
    feature <- cbind(mRNA[commonCell,], SHOTGUN[commonCell,], genomics[commonCell,])
  }
  if (feat=="NCI60_mRNA.SHOTGUN.SWATH"){
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    commonCell <- intersect(rownames(mRNA), intersect(rownames(SHOTGUN), rownames(SWATH)))
    feature <- cbind(mRNA[commonCell,], SHOTGUN[commonCell,], SWATH[commonCell,])
  }
  if (feat=="NCI60_mRNA.genomics.SWATH"){
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    commonCell <- intersect(rownames(mRNA), intersect(rownames(genomics), rownames(SWATH)))
    feature <- cbind(mRNA[commonCell,], genomics[commonCell,], SWATH[commonCell,])
  }
  if (feat=="NCI60_mRNA.genomics.SHOTGUN_500"){
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN_500.ro")
    commonCell <- intersect(rownames(mRNA), intersect(rownames(genomics), rownames(SHOTGUN_500)))
    feature <- cbind(mRNA[commonCell,], genomics[commonCell,], SHOTGUN_500[commonCell,])
  }
  if (feat=="NCI60_mRNA.SWATH.SHOTGUN_500"){
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN_500.ro")
    commonCell <- intersect(rownames(mRNA), intersect(rownames(SWATH), rownames(SHOTGUN_500)))
    feature <- cbind(mRNA[commonCell,], SWATH[commonCell,], SHOTGUN_500[commonCell,])
  }
  if (feat=="NCI60_SHOTGUN.genomics.SWATH"){
    load("data/pipeline/NCI60_final_v2/SHOTGUN.ro")
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    commonCell <- intersect(rownames(SHOTGUN), intersect(rownames(genomics), rownames(SWATH)))
    feature <- cbind(SHOTGUN[commonCell,], genomics[commonCell,], SWATH[commonCell,])
  }
  if (feat=="NCI60_genomics.SWATH.SHOTGUN_500"){
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN_500.ro")
    commonCell <- intersect(rownames(genomics), intersect(rownames(SWATH), rownames(SHOTGUN_500)))
    feature <- cbind(genomics[commonCell,], SWATH[commonCell,], SHOTGUN_500[commonCell,])
  }
  if (feat=="NCI60_TISSUE.RPPA.mRNA.SHOTGUN"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN.ro")
    commonCell <- intersect(rownames(TISSUE), intersect(rownames(RPPA), intersect(rownames(mRNA), rownames(SHOTGUN))))
    feature <- cbind(TISSUE[commonCell,], RPPA[commonCell,], mRNA[commonCell,], SHOTGUN[commonCell,])
  }
  if (feat=="NCI60_TISSUE.RPPA.mRNA.genomics"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    commonCell <- intersect(rownames(TISSUE), intersect(rownames(RPPA), intersect(rownames(mRNA), rownames(genomics))))
    feature <- cbind(TISSUE[commonCell,], RPPA[commonCell,], mRNA[commonCell,], genomics[commonCell,])
  }
  if (feat=="NCI60_TISSUE.RPPA.mRNA.SWATH"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    commonCell <- intersect(rownames(TISSUE), intersect(rownames(RPPA), intersect(rownames(mRNA), rownames(SWATH))))
    feature <- cbind(TISSUE[commonCell,], RPPA[commonCell,], mRNA[commonCell,], SWATH[commonCell,])
  }
  if (feat=="NCI60_TISSUE.RPPA.mRNA.SHOTGUN_500"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN_500.ro")
    commonCell <- intersect(rownames(TISSUE), intersect(rownames(RPPA), intersect(rownames(mRNA), rownames(SHOTGUN_500))))
    feature <- cbind(TISSUE[commonCell,], RPPA[commonCell,], mRNA[commonCell,], SHOTGUN_500[commonCell,])
  }
  if (feat=="NCI60_TISSUE.RPPA.SHOTGUN.genomics"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN.ro")
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    commonCell <- intersect(rownames(TISSUE), intersect(rownames(RPPA), intersect(rownames(SHOTGUN), rownames(genomics))))
    feature <- cbind(TISSUE[commonCell,], RPPA[commonCell,], SHOTGUN[commonCell,], genomics[commonCell,])
  }
  if (feat=="NCI60_TISSUE.RPPA.SHOTGUN.SWATH"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    commonCell <- intersect(rownames(TISSUE), intersect(rownames(RPPA), intersect(rownames(SHOTGUN), rownames(SWATH))))
    feature <- cbind(TISSUE[commonCell,], RPPA[commonCell,], SHOTGUN[commonCell,], SWATH[commonCell,])
  }
  if (feat=="NCI60_TISSUE.RPPA.genomics.SWATH"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    commonCell <- intersect(rownames(TISSUE), intersect(rownames(RPPA), intersect(rownames(genomics), rownames(SWATH))))
    feature <- cbind(TISSUE[commonCell,], RPPA[commonCell,], genomics[commonCell,], SWATH[commonCell,])
  }
  if (feat=="NCI60_TISSUE.RPPA.genomics.SHOTGUN_500"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN_500.ro")
    commonCell <- intersect(rownames(TISSUE), intersect(rownames(RPPA), intersect(rownames(genomics), rownames(SHOTGUN_500))))
    feature <- cbind(TISSUE[commonCell,], RPPA[commonCell,], genomics[commonCell,], SHOTGUN_500[commonCell,])
  }
  if (feat=="NCI60_TISSUE.RPPA.SWATH.SHOTGUN_500"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN_500.ro")
    commonCell <- intersect(rownames(TISSUE), intersect(rownames(RPPA), intersect(rownames(SWATH), rownames(SHOTGUN_500))))
    feature <- cbind(TISSUE[commonCell,], RPPA[commonCell,], SWATH[commonCell,], SHOTGUN_500[commonCell,])
  }
  if (feat=="NCI60_TISSUE.mRNA.SHOTGUN.genomics"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN.ro")
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    commonCell <- intersect(rownames(TISSUE), intersect(rownames(mRNA), intersect(rownames(SHOTGUN), rownames(genomics))))
    feature <- cbind(TISSUE[commonCell,], mRNA[commonCell,], SHOTGUN[commonCell,], genomics[commonCell,])
  }
  if (feat=="NCI60_TISSUE.mRNA.SHOTGUN.SWATH"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    commonCell <- intersect(rownames(TISSUE), intersect(rownames(mRNA), intersect(rownames(SHOTGUN), rownames(SWATH))))
    feature <- cbind(TISSUE[commonCell,], mRNA[commonCell,], SHOTGUN[commonCell,], SWATH[commonCell,])
  }
  if (feat=="NCI60_TISSUE.mRNA.genomics.SWATH"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    commonCell <- intersect(rownames(TISSUE), intersect(rownames(mRNA), intersect(rownames(genomics), rownames(SWATH))))
    feature <- cbind(TISSUE[commonCell,], mRNA[commonCell,], genomics[commonCell,], SWATH[commonCell,])
  }
  if (feat=="NCI60_TISSUE.mRNA.genomics.SHOTGUN_500"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN_500.ro")
    commonCell <- intersect(rownames(TISSUE), intersect(rownames(mRNA), intersect(rownames(genomics), rownames(SHOTGUN_500))))
    feature <- cbind(TISSUE[commonCell,], mRNA[commonCell,], genomics[commonCell,], SHOTGUN_500[commonCell,])
  }
  if (feat=="NCI60_TISSUE.mRNA.SWATH.SHOTGUN_500"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN_500.ro")
    commonCell <- intersect(rownames(TISSUE), intersect(rownames(mRNA), intersect(rownames(SWATH), rownames(SHOTGUN_500))))
    feature <- cbind(TISSUE[commonCell,], mRNA[commonCell,], SWATH[commonCell,], SHOTGUN_500[commonCell,])
  }
  if (feat=="NCI60_TISSUE.SHOTGUN.genomics.SWATH"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN.ro")
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    commonCell <- intersect(rownames(TISSUE), intersect(rownames(SHOTGUN), intersect(rownames(genomics), rownames(SWATH))))
    feature <- cbind(TISSUE[commonCell,], SHOTGUN[commonCell,], genomics[commonCell,], SWATH[commonCell,])
  }
  if (feat=="NCI60_TISSUE.genomics.SWATH.SHOTGUN_500"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN_500.ro")
    commonCell <- intersect(rownames(TISSUE), intersect(rownames(genomics), intersect(rownames(SWATH), rownames(SHOTGUN_500))))
    feature <- cbind(TISSUE[commonCell,], genomics[commonCell,], SWATH[commonCell,], SHOTGUN_500[commonCell,])
  }
  if (feat=="NCI60_RPPA.mRNA.SHOTGUN.genomics"){
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN.ro")
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    commonCell <- intersect(rownames(RPPA), intersect(rownames(mRNA), intersect(rownames(SHOTGUN), rownames(genomics))))
    feature <- cbind(RPPA[commonCell,], mRNA[commonCell,], SHOTGUN[commonCell,], genomics[commonCell,])
  }
  if (feat=="NCI60_RPPA.mRNA.SHOTGUN.SWATH"){
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    commonCell <- intersect(rownames(RPPA), intersect(rownames(mRNA), intersect(rownames(SHOTGUN), rownames(SWATH))))
    feature <- cbind(RPPA[commonCell,], mRNA[commonCell,], SHOTGUN[commonCell,], SWATH[commonCell,])
  }
  if (feat=="NCI60_RPPA.mRNA.genomics.SWATH"){
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    commonCell <- intersect(rownames(RPPA), intersect(rownames(mRNA), intersect(rownames(genomics), rownames(SWATH))))
    feature <- cbind(RPPA[commonCell,], mRNA[commonCell,], genomics[commonCell,], SWATH[commonCell,])
  }
  if (feat=="NCI60_RPPA.mRNA.genomics.SHOTGUN_500"){
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN_500.ro")
    commonCell <- intersect(rownames(RPPA), intersect(rownames(mRNA), intersect(rownames(genomics), rownames(SHOTGUN_500))))
    feature <- cbind(RPPA[commonCell,], mRNA[commonCell,], genomics[commonCell,], SHOTGUN_500[commonCell,])
  }
  if (feat=="NCI60_RPPA.mRNA.SWATH.SHOTGUN_500"){
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN_500.ro")
    commonCell <- intersect(rownames(RPPA), intersect(rownames(mRNA), intersect(rownames(SWATH), rownames(SHOTGUN_500))))
    feature <- cbind(RPPA[commonCell,], mRNA[commonCell,], SWATH[commonCell,], SHOTGUN_500[commonCell,])
  }
  if (feat=="NCI60_RPPA.SHOTGUN.genomics.SWATH"){
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN.ro")
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    commonCell <- intersect(rownames(RPPA), intersect(rownames(SHOTGUN), intersect(rownames(genomics), rownames(SWATH))))
    feature <- cbind(RPPA[commonCell,], SHOTGUN[commonCell,], genomics[commonCell,], SWATH[commonCell,])
  }
  if (feat=="NCI60_RPPA.genomics.SWATH.SHOTGUN_500"){
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN_500.ro")
    commonCell <- intersect(rownames(RPPA), intersect(rownames(genomics), intersect(rownames(SWATH), rownames(SHOTGUN_500))))
    feature <- cbind(RPPA[commonCell,], genomics[commonCell,], SWATH[commonCell,], SHOTGUN_500[commonCell,])
  }
  if (feat=="NCI60_mRNA.SHOTGUN.genomics.SWATH"){
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN.ro")
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    commonCell <- intersect(rownames(mRNA), intersect(rownames(SHOTGUN), intersect(rownames(genomics), rownames(SWATH))))
    feature <- cbind(mRNA[commonCell,], SHOTGUN[commonCell,], genomics[commonCell,], SWATH[commonCell,])
  }
  if (feat=="NCI60_mRNA.genomics.SWATH.SHOTGUN_500"){
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN_500.ro")
    commonCell <- intersect(rownames(mRNA), intersect(rownames(genomics), intersect(rownames(SWATH), rownames(SHOTGUN_500))))
    feature <- cbind(mRNA[commonCell,], genomics[commonCell,], SWATH[commonCell,], SHOTGUN_500[commonCell,])
  }
  if (feat=="NCI60_TISSUE.RPPA.mRNA.SHOTGUN.genomics"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN.ro")
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    commonCell <- intersect(rownames(TISSUE), intersect(rownames(RPPA), intersect(rownames(mRNA), intersect(rownames(SHOTGUN), rownames(genomics)))))
    feature <- cbind(TISSUE[commonCell,], RPPA[commonCell,], mRNA[commonCell,], SHOTGUN[commonCell,], genomics[commonCell,])
  }
  if (feat=="NCI60_TISSUE.RPPA.mRNA.SHOTGUN.SWATH"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    commonCell <- intersect(rownames(TISSUE), intersect(rownames(RPPA), intersect(rownames(mRNA), intersect(rownames(SHOTGUN), rownames(SWATH)))))
    feature <- cbind(TISSUE[commonCell,], RPPA[commonCell,], mRNA[commonCell,], SHOTGUN[commonCell,], SWATH[commonCell,])
  }
  if (feat=="NCI60_TISSUE.RPPA.mRNA.genomics.SWATH"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    commonCell <- intersect(rownames(TISSUE), intersect(rownames(RPPA), intersect(rownames(mRNA), intersect(rownames(genomics), rownames(SWATH)))))
    feature <- cbind(TISSUE[commonCell,], RPPA[commonCell,], mRNA[commonCell,], genomics[commonCell,], SWATH[commonCell,])
  }
  if (feat=="NCI60_TISSUE.RPPA.mRNA.genomics.SHOTGUN_500"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN_500.ro")
    commonCell <- intersect(rownames(TISSUE), intersect(rownames(RPPA), intersect(rownames(mRNA), intersect(rownames(genomics), rownames(SHOTGUN_500)))))
    feature <- cbind(TISSUE[commonCell,], RPPA[commonCell,], mRNA[commonCell,], genomics[commonCell,], SHOTGUN_500[commonCell,])
  }
  if (feat=="NCI60_TISSUE.RPPA.mRNA.SWATH.SHOTGUN_500"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN_500.ro")
    commonCell <- intersect(rownames(TISSUE), intersect(rownames(RPPA), intersect(rownames(mRNA), intersect(rownames(SWATH), rownames(SHOTGUN_500)))))
    feature <- cbind(TISSUE[commonCell,], RPPA[commonCell,], mRNA[commonCell,], SWATH[commonCell,], SHOTGUN_500[commonCell,])
  }
  if (feat=="NCI60_TISSUE.RPPA.SHOTGUN.genomics.SWATH"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN.ro")
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    commonCell <- intersect(rownames(TISSUE), intersect(rownames(RPPA), intersect(rownames(SHOTGUN), intersect(rownames(genomics), rownames(SWATH)))))
    feature <- cbind(TISSUE[commonCell,], RPPA[commonCell,], SHOTGUN[commonCell,], genomics[commonCell,], SWATH[commonCell,])
  }
  if (feat=="NCI60_TISSUE.RPPA.genomics.SWATH.SHOTGUN_500"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN_500.ro")
    commonCell <- intersect(rownames(TISSUE), intersect(rownames(RPPA), intersect(rownames(genomics), intersect(rownames(SWATH), rownames(SHOTGUN_500)))))
    feature <- cbind(TISSUE[commonCell,], RPPA[commonCell,], genomics[commonCell,], SWATH[commonCell,], SHOTGUN_500[commonCell,])
  }
  if (feat=="NCI60_TISSUE.mRNA.SHOTGUN.genomics.SWATH"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN.ro")
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    commonCell <- intersect(rownames(TISSUE), intersect(rownames(mRNA), intersect(rownames(SHOTGUN), intersect(rownames(genomics), rownames(SWATH)))))
    feature <- cbind(TISSUE[commonCell,], mRNA[commonCell,], SHOTGUN[commonCell,], genomics[commonCell,], SWATH[commonCell,])
  }
  if (feat=="NCI60_TISSUE.mRNA.genomics.SWATH.SHOTGUN_500"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN_500.ro")
    commonCell <- intersect(rownames(TISSUE), intersect(rownames(mRNA), intersect(rownames(genomics), intersect(rownames(SWATH), rownames(SHOTGUN_500)))))
    feature <- cbind(TISSUE[commonCell,], mRNA[commonCell,], genomics[commonCell,], SWATH[commonCell,], SHOTGUN_500[commonCell,])
  }
  if (feat=="NCI60_RPPA.mRNA.SHOTGUN.genomics.SWATH"){
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN.ro")
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    commonCell <- intersect(rownames(RPPA), intersect(rownames(mRNA), intersect(rownames(SHOTGUN), intersect(rownames(genomics), rownames(SWATH)))))
    feature <- cbind(RPPA[commonCell,], mRNA[commonCell,], SHOTGUN[commonCell,], genomics[commonCell,], SWATH[commonCell,])
  }
  if (feat=="NCI60_RPPA.mRNA.genomics.SWATH.SHOTGUN_500"){
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN_500.ro")
    commonCell <- intersect(rownames(RPPA), intersect(rownames(mRNA), intersect(rownames(genomics), intersect(rownames(SWATH), rownames(SHOTGUN_500)))))
    feature <- cbind(RPPA[commonCell,], mRNA[commonCell,], genomics[commonCell,], SWATH[commonCell,], SHOTGUN_500[commonCell,])
  }
  if (feat=="NCI60_TISSUE.RPPA.mRNA.SHOTGUN.genomics.SWATH"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN.ro")
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    commonCell <- intersect(rownames(TISSUE), intersect(rownames(RPPA), intersect(rownames(mRNA), intersect(rownames(SHOTGUN), intersect(rownames(genomics), rownames(SWATH))))))
    feature <- cbind(TISSUE[commonCell,], RPPA[commonCell,], mRNA[commonCell,], SHOTGUN[commonCell,], genomics[commonCell,], SWATH[commonCell,])
  }
  if (feat=="NCI60_TISSUE.RPPA.mRNA.genomics.SWATH.SHOTGUN_500"){
    load("data/pipeline/NCI60_final_v2/TISSUE.ro")
    load("data/pipeline/NCI60_final_v2/RPPA.ro")
    load("data/pipeline/NCI60_final_v2/mRNA.ro")
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    load("data/pipeline/NCI60_final_v2/SHOTGUN_500.ro")
    commonCell <- intersect(rownames(TISSUE), intersect(rownames(RPPA), intersect(rownames(mRNA), intersect(rownames(genomics), intersect(rownames(SWATH), rownames(SHOTGUN_500))))))
    feature <- cbind(TISSUE[commonCell,], RPPA[commonCell,], mRNA[commonCell,], genomics[commonCell,], SWATH[commonCell,], SHOTGUN_500[commonCell,])
  }
  
  
  
  #############################################################################################
  # Test
  #############################################################################################
  if (feat=="SWATH_paper"){
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    feature <- SWATH
  }

  if (feat=="SWATH_log10"){
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    feature <- 10^SWATH
  }
  
  if (feat=="SWATH_ln"){
    load("data/pipeline/NCI60_final_v2/SWATH.ro")
    feature <- exp(1)^SWATH
  }
  
  if (feat=="SWATH_PROT_QUANT_R1"){
    load("data/pipeline/NCI60_final_v2/SWATH_PROT_QUANT_R1.ro")
    feature <- PROT_QUANT_R1
  }
  
  if (feat=="SWATH_PROT_QUANT_R2"){
    load("data/pipeline/NCI60_final_v2/SWATH_PROT_QUANT_R2.ro")
    feature <- PROT_QUANT_R2
  }
  
  if (feat=="gGDSC"){
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    feature <- genomics[,1:70]
  }
  
  if (feat=="gNCI"){
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    feature <- genomics[,71:ncol(genomics)]
  }
  
  if (feat=="gNCI_gGDSC"){
    load("data/pipeline/NCI60_final_v2/genomics.ro")
    feature <- genomics
  }
  
  #############################################################################################
  # DREAM pilot
  #############################################################################################
  
  # as In Sock + tissue
  if (feat=="DREAM_BEM_2src"){
    load("data/pipeline/DREAM/BEM_2src.ro")
    feature <- BEM_2src
  }
  if (feat=="DREAM_gex"){
    load("data/pipeline/DREAM/gex.ro")
    feature <- gex
  }
  if (feat=="DREAM_methyl"){
    load("data/pipeline/DREAM/methyl.ro")
    feature <- methyl
  }
  if (feat=="DREAM_BEM_2src.gex"){
    load("data/pipeline/DREAM/BEM_2src.ro")
    load("data/pipeline/DREAM/gex.ro")
    feature <- cbind(BEM_2src, gex)
  }
  if (feat=="DREAM_BEM_2src.methyl"){
    load("data/pipeline/DREAM/BEM_2src.ro")
    load("data/pipeline/DREAM/methyl.ro")
    feature <- cbind(BEM_2src, methyl)
  }
  if (feat=="DREAM_methyl.gex"){
    load("data/pipeline/DREAM/methyl.ro")
    load("data/pipeline/DREAM/gex.ro")
    feature <- cbind(methyl, gex)
  }
  if (feat=="DREAM_BEM_2src.methyl.gex"){
    load("data/pipeline/DREAM/BEM_2src.ro")
    load("data/pipeline/DREAM/methyl.ro")
    load("data/pipeline/DREAM/gex.ro")
    feature <- cbind(BEM_2src, methyl, gex)
  }
  # add tissue
  if (feat=="DREAM_tissue"){
    load("data/pipeline/DREAM/tissue.ro")
    feature <- tissue
  }
  if (feat=="DREAM_tissue.BEM_2src"){
    load("data/pipeline/DREAM/tissue.ro")
    load("data/pipeline/DREAM/BEM_2src.ro")
    feature <- cbind(tissue, BEM_2src)
  }
  if (feat=="DREAM_tissue.gex"){
    load("data/pipeline/DREAM/tissue.ro")
    load("data/pipeline/DREAM/gex.ro")
    feature <- cbind(tissue, gex)
  }
  if (feat=="DREAM_tissue.methyl"){
    load("data/pipeline/DREAM/tissue.ro")
    load("data/pipeline/DREAM/methyl.ro")
    feature <- cbind(tissue, methyl)
  }
  if (feat=="DREAM_tissue.BEM_2src.gex"){
    load("data/pipeline/DREAM/tissue.ro")
    load("data/pipeline/DREAM/BEM_2src.ro")
    load("data/pipeline/DREAM/gex.ro")
    feature <- cbind(tissue, BEM_2src, gex)
  }
  if (feat=="DREAM_tissue.BEM_2src.methyl"){
    load("data/pipeline/DREAM/tissue.ro")
    load("data/pipeline/DREAM/BEM_2src.ro")
    load("data/pipeline/DREAM/methyl.ro")
    feature <- cbind(tissue, BEM_2src, methyl)
  }
  if (feat=="DREAM_tissue.methyl.gex"){
    load("data/pipeline/DREAM/tissue.ro")
    load("data/pipeline/DREAM/methyl.ro")
    load("data/pipeline/DREAM/gex.ro")
    feature <- cbind(tissue, methyl, gex)
  }
  if (feat=="DREAM_tissue.BEM_2src.methyl.gex"){
    load("data/pipeline/DREAM/tissue.ro")
    load("data/pipeline/DREAM/BEM_2src.ro")
    load("data/pipeline/DREAM/methyl.ro")
    load("data/pipeline/DREAM/gex.ro")
    feature <- cbind(tissue, BEM_2src, methyl, gex)
  }
  # as 2nd BEM implementation
  if (feat=="DREAM_BEM_union"){
    load("data/pipeline/DREAM/BEM_union.ro")
    feature <- BEM_union
  }
  if (feat=="DREAM_BEM_union.gex"){
    load("data/pipeline/DREAM/BEM_union.ro")
    load("data/pipeline/DREAM/gex.ro")
    feature <- cbind(BEM_union, gex)
  }
  if (feat=="DREAM_BEM_union.methyl"){
    load("data/pipeline/DREAM/BEM_union.ro")
    load("data/pipeline/DREAM/methyl.ro")
    feature <- cbind(BEM_union, methyl)
  }
  if (feat=="DREAM_BEM_union.methyl.gex"){
    load("data/pipeline/DREAM/BEM_union.ro")
    load("data/pipeline/DREAM/methyl.ro")
    load("data/pipeline/DREAM/gex.ro")
    feature <- cbind(BEM_union, methyl, gex)
  }
  if (feat=="DREAM_tissue.BEM_union"){
    load("data/pipeline/DREAM/tissue.ro")
    load("data/pipeline/DREAM/BEM_union.ro")
    feature <- cbind(tissue, BEM_union)
  }
  if (feat=="DREAM_tissue.BEM_union.gex"){
    load("data/pipeline/DREAM/tissue.ro")
    load("data/pipeline/DREAM/BEM_union.ro")
    load("data/pipeline/DREAM/gex.ro")
    feature <- cbind(tissue, BEM_union, gex)
  }
  if (feat=="DREAM_tissue.BEM_union.methyl"){
    load("data/pipeline/DREAM/tissue.ro")
    load("data/pipeline/DREAM/BEM_union.ro")
    load("data/pipeline/DREAM/methyl.ro")
    feature <- cbind(tissue, BEM_union, methyl)
  }
  if (feat=="DREAM_tissue.BEM_union.methyl.gex"){
    load("data/pipeline/DREAM/tissue.ro")
    load("data/pipeline/DREAM/BEM_union.ro")
    load("data/pipeline/DREAM/methyl.ro")
    load("data/pipeline/DREAM/gex.ro")
    feature <- cbind(tissue, BEM_union, methyl, gex)
  }
  
  
  #############################################################################################
  # feature which contains an NA will be removed.
  #############################################################################################
  feature <- feature[apply(feature, 1, function(x) sum(is.na(x))!=length(x)),]
  feature <- feature[,apply(feature, 2, function(x) !any(is.na(x)))]
  return(feature)
}
