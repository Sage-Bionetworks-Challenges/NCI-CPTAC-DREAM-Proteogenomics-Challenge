
# path <- "/Users/miyang/Documents/RWTH_Aachen/DREAM_CPTAC/DOCKER_STORAGE/docker_R_sc2/"
path <- "/"

# load package and functions
library(glmnet)
predictEN <- function(model, lambda, feat) {
  return(predict(model,type="response",newx=feat)[,model$lambda == lambda])
}

# load saved models
load(paste0(path,"model_storage/SAVE_MODEL_ElasticNet_breast_PROT_RNAseq_fold5_ite10.Rdata"))
model["VGF"] <- NULL

# load testing data, to check data name go to this page: https://www.synapse.org/#!Synapse:syn8228304/wiki/448379
breast_rna_EVAL  <- read.csv(paste0(path,"evaluation_data/prospective_breast_RNA_sort_common_gene_15107.txt"), row.names= 1, sep="\t", check.names = F )
features <- t(breast_rna_EVAL)

prediction <- c()
for(i in 1:length(model)) {
  fit <- model[[i]]
  test_pred <- predictEN(fit,fit$lambda, features ) 
  prediction <- rbind(prediction, test_pred)
}

# some predictions were all zeros, so I used some mRNA data to replace the non valid predictions.
prediction[which(apply(prediction, 1, var) == 0),  ] <- features[ , 1:length(which(apply(prediction, 1, var) == 0)) ]

# ensure the format is correct
prediction <- cbind(proteinID = names(model),prediction)

# save the prediction matrix
write.table(prediction, file = paste0(path,"output/predictions.tsv"), sep="\t" ,row.names=F)
# save the confidence matrix 
write.table(prediction, file = paste0(path,"output/confidence.tsv"), sep="\t" ,row.names=F)



