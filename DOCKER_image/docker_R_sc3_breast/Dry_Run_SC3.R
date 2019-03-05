
# path <- "/Users/miyang/Documents/RWTH_Aachen/DREAM_CPTAC/DOCKER_STORAGE/docker_R_sc3_breast/"
path <- "/"

# load package and functions
library(randomForest)
predictRF <- function(model,newx) {
  return(predict(model,newdata=newx,type="response"))
}

############################################### BREAST ###############################################

# load saved models
out.new<-vector("list", 32313)
phosphoID<-NULL
for (j in 1:324) {
  load(paste0(path,"model_storage/Breast/Model_",j,".rda"))
  phosphoID<-c(phosphoID,phospho.reg)
  for (s in 1:length(out))  out.new[[(j-1)*100+s]]<-out[[s]]
  print(j)
}


# load testing data, to check data name go to this page: https://www.synapse.org/#!Synapse:syn8228300/wiki/413417
features  <- read.csv(paste0(path,"evaluation_data/prospective_breast_proteome_sort_common_gene_10005.txt"), row.names= 1, sep="\t", check.names = F )

# select features used in trained model
load(paste0(path,"model_storage/Breast/protein_subset.rda")) # --- load proteins with complete measurement

features<-features[match(protein.subset,rownames(features)), ]
rownames(features) <- NULL

prediction <- c()
for(i in 1:length(out.new)) {
  fit <- out.new[[i]]
  test_pred <- predictRF(fit, t(features) )
  prediction <- rbind(prediction, test_pred)
}

# ensure the format is correct
prediction <- cbind(phosphoID = phosphoID, prediction)

# save the prediction matrix
write.table(prediction, file = paste0(path,"output/predictions.tsv"), sep="\t",row.names=F )
# save the confidence matrix 
write.table(prediction, file = paste0(path,"output/confidence.tsv"), sep="\t" ,row.names=F)


