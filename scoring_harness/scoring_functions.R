
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

########################################## load NRMSE function #########################################
NRMSE_by_row <- function(pred_path, truth_path)  {
  suppressPackageStartupMessages(library(hydroGOF))
  prediction <- read.csv( pred_path, row.names = 1 , check.names = F, sep="\t") 
  test_prot  <- read.csv( truth_path, row.names = 1 , check.names = F, sep="\t")
  common_protein <- intersect(rownames(prediction), rownames(test_prot))
  prediction <- prediction[common_protein ,]
  test_prot <- test_prot[common_protein ,]
  mat1 <- as.matrix(prediction)
  mat2 <- as.matrix(test_prot) 
  nrmse_vec <- c()
  for(i in 1:length(mat1[ ,1]) ) {
    temp <- hydroGOF::rmse(mat1[ i, ], mat2[ i , ])
    nrmse_vec <- c(nrmse_vec , temp/(max(mat2,na.rm = T)-min(mat2,na.rm = T)))
  }
  names(nrmse_vec) <- rownames(mat1)
  return(mean(nrmse_vec))
}

#result_nrmse <- NRMSE_by_row("predictions.tsv", "pros_ova_proteome_sort_common_gene_6577.txt")

