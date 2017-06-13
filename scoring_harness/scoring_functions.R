
################################## load pearson correlation function ##################################
correlation_by_row <- function(pred_path, truth_path)  {
  prediction <- read.csv( pred_path, row.names = 1 , check.names = F, sep="\t") 
  test_prot  <- read.csv( truth_path, row.names = 1 , check.names = F, sep="\t")
  common_protein <- intersect(rownames(prediction), rownames(test_prot))
  prediction <- prediction[common_protein ,]
  test_prot <- test_prot[common_protein ,]

  mat1 <- as.matrix(prediction)
  mat2 <- as.matrix(test_prot) 
  corr_vec <- c()
  for(i in 1:length(mat1[ ,1]) ) {
    temp <- cor.test(mat1[ i, ], mat2[ i , ])
    pcorr <- temp$estimate # pearson correlation
    corr_vec <- c(corr_vec , pcorr)
  }
  names(corr_vec) <- rownames(mat1)
  return(mean(corr_vec))
}

#result_corr <- correlation_by_row("predictions.tsv", "pros_ova_proteome_sort_common_gene_6577.txt")

########################################## load RMSE function #########################################
RMSE_by_row <- function(pred_path, truth_path)  {
  suppressPackageStartupMessages(library(hydroGOF))
  prediction <- read.csv( pred_path, row.names = 1 , check.names = F, sep="\t") 
  test_prot  <- read.csv( truth_path, row.names = 1 , check.names = F, sep="\t")
  common_protein <- intersect(rownames(prediction), rownames(test_prot))
  prediction <- prediction[common_protein ,]
  test_prot <- test_prot[common_protein ,]
  mat1 <- as.matrix(prediction)
  mat2 <- as.matrix(test_prot) 
  rmse_vec <- c()
  for(i in 1:length(mat1[ ,1]) ) {
    temp <- hydroGOF::rmse(mat1[ i, ], mat2[ i , ])
    rmse_vec <- c(rmse_vec , temp)
  }
  names(rmse_vec) <- rownames(mat1)
  return(mean(rmse_vec))
}

#result_rmse <- RMSE_by_row("predictions.tsv", "pros_ova_proteome_sort_common_gene_6577.txt")

