library(hydroGOF)

################################## load pearson correlation function ##################################
correlation_by_column <- function(pred_path, test_path)  {
  prediction <- read.csv( pred_path, row.names = 1 ,sep="\t",check.names = F) 
  test_prot  <- read.csv( test_path, row.names = 1 ,sep="\t",check.names = F)
  common_protein <- intersect(rownames(prediction), rownames(test_prot))
  prediction <- prediction[common_protein,]
  test_prot <- test_prot[common_protein,]

  mat1 <- as.matrix(prediction)
  mat2 <- as.matrix(test_prot) 
  corr_vec <- c()
  for(i in 1:length(mat1[1,]) ) {
    temp <- cor.test(mat1[ , i], mat2[ , i ])
    pcorr <- temp$estimate # pearson correlation
    corr_vec <- c(corr_vec , pcorr)
  }
  names(corr_vec) <- colnames(mat1)
  return(corr_vec)
}
#result_corr <- correlation_by_column("scoring_harness/predictions.tsv", "scoring_harness/pros_ova_proteome_sort_common_gene_6577.txt")

########################################## load RMSE function #########################################
RMSE_by_column <- function(pred_path, test_path)  {
  prediction <- read.csv( pred_path, row.names = 1  ,sep="\t",check.names = F) 
  test_prot  <- read.csv( test_path, row.names = 1 ,sep="\t",check.names = F)
  common_protein <- intersect(rownames(prediction), rownames(test_prot))
  prediction <- prediction[ common_protein, ]
  test_prot <- test_prot[ common_protein, ]

  mat1 <- as.matrix(prediction)
  mat2 <- as.matrix(test_prot) 
  rmse_vec <- c()
  for(i in 1:length(mat1[1,]) ) {
    temp <- rmse(mat1[ , i ], mat2[ , i ])
    rmse_vec <- c(rmse_vec , temp)
  }
  names(rmse_vec) <- colnames(mat1)
  return(rmse_vec)
}

#result_rmse <- RMSE_by_column("scoring_harness/predictions.tsv", "scoring_harness/pros_ova_proteome_sort_common_gene_6577.txt")


