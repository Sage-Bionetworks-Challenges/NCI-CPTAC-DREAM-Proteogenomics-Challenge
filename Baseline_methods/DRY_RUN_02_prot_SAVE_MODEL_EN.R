library("doParallel")
cl <- makeCluster(72)
registerDoParallel(cl)

path <- "/home/my871390/MI_YANG/RWTH_Aachen"
# path <- "~/Documents/RWTH_Aachen"
source(paste(path,"/ML_CODE/crossvalidate_data/CrossVal_LIBRARY.R", sep=""))
source(paste(path,"/ML_CODE/ML/LIB_METRICS.R", sep=""))
source(paste(path,"/FUNCTIONS/EN_functions.R", sep=""))
source(paste(path,"/FUNCTIONS/general_functions.R", sep=""))
result_folder <- paste(path,"/DREAM_CPTAC/PROT_PREDICTION/DATA_RESULT_STORAGE/", sep="")

################################################ HGSC ################################################
######################################################################################################

####### load data
protein_ovarian <- read.csv( paste(path,"/DREAM_CPTAC/PROT_PREDICTION/DATA/protein_ovarian", sep=""), row.names = 1, check.names = F )
features <- read.csv(paste(path,"/DREAM_CPTAC/PROT_PREDICTION/DATA/rna_ovarian", sep=""),check.names = F,row.names=1)

####### model
RES <- as.matrix(protein_ovarian)
features <- as.matrix(features)
common_sample <- intersect(rownames(features),colnames(RES)) ; RES <- RES[ ,common_sample] ; features <- features[common_sample, ]

fold=5 ; iteration=10 ; alpha = 0.5 ; N <- length(rownames(RES))

model <- foreach(i = 1:N ) %dopar% {

  library('methods')
  path <- "/home/my871390/MI_YANG/RWTH_Aachen"
  # path <- "~/Documents/RWTH_Aachen"
  source(paste(path,"/ML_CODE/crossvalidate_data/CrossVal_LIBRARY.R", sep=""))
  source(paste(path,"/ML_CODE/ML/LIB_METRICS.R", sep=""))
  source(paste(path,"/FUNCTIONS/EN_functions.R", sep=""))
  source(paste(path,"/FUNCTIONS/general_functions.R", sep=""))

  ####### load data
  protein_ovarian <- read.csv( paste(path,"/DREAM_CPTAC/PROT_PREDICTION/DATA/protein_ovarian", sep=""), row.names = 1, check.names = F )
  features <- read.csv(paste(path,"/DREAM_CPTAC/PROT_PREDICTION/DATA/rna_ovarian", sep=""),check.names = F,row.names=1)

  ####### model
  RES <- as.matrix(protein_ovarian)
  features <- as.matrix(features)
  common_sample <- intersect(rownames(features),colnames(RES)) ; RES <- RES[ ,common_sample] ; features <- features[common_sample, ]

  model <- SAVE_MODEL ( i , fold=fold, iteration=iteration , alpha = alpha )

}
### save models
names(model) <- rownames(RES)[1:N] #
save(model, file=paste0(result_folder,"SAVE_MODEL_ElasticNet_HGSC_ALL_PROT_microarray_","fold",as.character(fold),"_ite",as.character(iteration),".Rdata"))

stopCluster( cl )



############################################### breast ################################################
#######################################################################################################

####### load data
protein_breast <- read.csv( paste(path,"/DREAM_CPTAC/PROT_PREDICTION/DATA/protein_breast_knn", sep=""), row.names = 1, check.names = F )
features <- read.csv(paste(path,"/DREAM_CPTAC/PROT_PREDICTION/DATA/rna_breast_knn", sep=""),check.names = F,row.names=1)

####### model 
RES <- as.matrix(protein_breast)
features <- t(as.matrix(features))
common_sample <- intersect(rownames(features),colnames(RES)) ; RES <- RES[ ,common_sample] ; features <- features[common_sample, ]

fold=5 ; iteration=10 ; alpha = 0.5 ; N <- length(rownames(RES))

model <- foreach(i = 1:N ) %dopar% { 
  
  library('methods')
  path <- "/home/my871390/MI_YANG/RWTH_Aachen"
  # path <- "~/Documents/RWTH_Aachen"
  source(paste(path,"/ML_CODE/crossvalidate_data/CrossVal_LIBRARY.R", sep=""))
  source(paste(path,"/ML_CODE/ML/LIB_METRICS.R", sep=""))
  source(paste(path,"/FUNCTIONS/EN_functions.R", sep=""))
  source(paste(path,"/FUNCTIONS/general_functions.R", sep=""))
  
  ####### load data
  protein_breast <- read.csv( paste(path,"/DREAM_CPTAC/PROT_PREDICTION/DATA/protein_breast_knn", sep=""), row.names = 1, check.names = F )
  features <- read.csv(paste(path,"/DREAM_CPTAC/PROT_PREDICTION/DATA/rna_breast_knn", sep=""),check.names = F,row.names=1)
  
  ####### model 
  RES <- as.matrix(protein_breast)
  features <- t(as.matrix(features))
  common_sample <- intersect(rownames(features),colnames(RES)) ; RES <- RES[ ,common_sample] ; features <- features[common_sample, ]
  
#   number_obs <- apply(RES, MARGIN = 1, FUN = function(x) length(x[!is.na(x)]) ) 
#   RES <- RES[which(number_obs < 4), ]  
  
  model <- SAVE_MODEL ( index=i , fold=fold, iteration=iteration , alpha =alpha )
  
}
names(model) <- rownames(RES) # 
save(model, file=paste0(result_folder,"SAVE_MODEL_ElasticNet_breast_PROT_RNAseq_","fold",as.character(fold),"_ite",as.character(iteration),".Rdata"))

stopCluster( cl )


