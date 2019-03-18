#--------------------------------------------------------------------------------------------------------------------------------
#DESC:  This function "create_train_plot" creates a plot containing two lines against nTrees. One line is represanting the 
#       xTrain_cIdxs and the other
#       one the train_cIdxs.
#IN:    store_path      ==> a path where the plot will be stored 
#       plot_name       ==> the name of the plot 
#       rf_model        ==> the RandomForest model with, which the cindexes have been calculated
#       cv_fold_number  ==> a numeric shwon in the name of the plot
#                           represanting for which fold of the crossvalidation the plot is created
#       drug_name       ==> the name of the drug from which the data has been generated
#OUT:   returns a png named after the value of "plot_name"
#--------------------------------------------------------------------------------------------------------------------------------
create_train_plot <-function(store_path, plot_name, rf_model, cv_fold_number, drug_name){
  
#save TrainPlots
#define the save_path for the xTrainPlots
save_xTrain_path<-file.path(store_path, plot_name)

png(file=save_xTrain_path)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(rf_model$nTrees, 
     rf_model$performance$xTrain_cIdx,
     xlab="nTree",
     ylab="Cindex",
     ylim=c(min(rf_model$performance$xTrain_cIdx, rf_model$performance$train_cIdx), 
            max(rf_model$performance$xTrain_cIdx, rf_model$performance$train_cIdx)),
     col="3",
     type="b", 
     main=drug_name,
     sub=paste("nTrees=",rf_model$nTree_best_performer))
lines(rf_model$nTrees,
      rf_model$performance$train_cIdx, 
      type="b", 
      col="2")
legend("topright", 
       inset=c(-0.3,0), 
       legend=c("xTrain_cIdx","train_cIdx"), 
       cex=0.9, 
       fill=c("3","2"))
dev.off()
}