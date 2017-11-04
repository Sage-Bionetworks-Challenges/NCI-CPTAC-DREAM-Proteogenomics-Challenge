##SC1

score.nrmsd = function(pred_path, observed_path, truth_path)  
{
  d.predict = as.matrix(read.csv( pred_path, sep="\t", row.names=1));
  d.true = as.matrix(read.csv( truth_path, sep="\t", row.names=1));
  d.obs = as.matrix(read.csv( observed_path, sep="\t",row.names = 1));
  d.obs[d.obs==''] = NA; 
  missing.ind = is.na(d.obs);
  
  diff.true = apply(d.true, 1, function(x){diff(range(x[x>0]))});
  d.predict[!missing.ind]= NA;
  nrmsd = sqrt(apply((d.predict-d.true)^2,1,mean,na.rm=T))/diff.true;
  return(mean(nrmsd));
}

score.cor = function(pred_path, observed_path, truth_path)  
{
  d.predict = as.matrix(read.csv( pred_path, sep="\t", row.names=1));
  d.true = as.matrix(read.csv( truth_path, sep="\t", row.names=1));
  d.obs = as.matrix(read.csv( observed_path, sep="\t",row.names = 1));
  d.obs[d.obs==''] = NA; 
  missing.ind = is.na(d.obs);
    
  L = dim(d.true)[1]
  d.predict[!missing.ind]= NA;
  d.predict[d.true==0]= NA;
  cor.p = sapply(1:L,function(l){cor(d.predict[l,],d.true[l,],use = 'pairwise.complete.obs')});
  cor.p[is.na(cor.p)]=0;
  return(mean(cor.p));
}
#get.score.sc1 isn't used in the scoring harness
get.score.sc1 = function(path_pred='/',path_obs='/',path_true='/')
{
  nrmsd.all=c()
  cor.all=c()
  
  for(i in 1:100)
  {
    pred_file = paste0(path_pred,sprintf('predictions_%s.tsv', i));
    observed_file = paste0(path_obs,sprintf('data_test_obs_%s.txt', i));
    truth_file = paste0(path_obs,sprintf('data_test_true_%s.txt', i));
    nrmsd.all = c(nrmsd.all,score.nrmsd(pred_file, observed_file, truth_file));  
    cor.all = c(cor.all,score.cor(pred_file, observed_file, truth_file));  
  }  
  
  return(list(nrmsd = mean(nrmsd.all),  cor = mean(cor.all) )); 
}


###  SC2 
################################## load pearson correlation function ##################################
correlation_by_row <- function(pred_path, truth_path) {
  prediction <- read.csv( pred_path, row.names = 1 , check.names = F, sep="\t") 
  test_prot  <- read.csv( truth_path, row.names = 1 , check.names = F, sep="\t")
  #common_protein <- intersect(rownames(prediction), rownames(test_prot))
  prediction <- prediction[rownames(test_prot), colnames(test_prot)]
  test_prot <- test_prot[rownames(test_prot), colnames(test_prot)]

  mat1 <- as.matrix(prediction)
  mat2 <- as.matrix(test_prot) 
  
  corr_vec <- c()
  for(i in 1:length(mat1[ ,1]) ) {
    c <- rbind(mat1[i, ], mat2[i, ]) ; c <- c[ ,complete.cases(t(c))]
    temp <- cor.test(mat1[ i, ], mat2[ i , ])
    pcorr <- temp$estimate # pearson correlation
    if (is.na(pcorr)) {pcorr<-0}
    corr_vec <- c(corr_vec , pcorr)
  }
  names(corr_vec) <- rownames(mat1)
  return(mean(corr_vec))
}
#result_corr <- correlation_by_row("predictions.tsv", "prospective_ova_proteome_sort_common_gene_7061.txt")

########################################## load NRMSE function #########################################
NRMSE_by_row <- function(pred_path, truth_path)  {
  suppressPackageStartupMessages(library(hydroGOF))
  prediction <- read.csv( pred_path, row.names = 1 , check.names = F, sep="\t") 
  test_prot  <- read.csv( truth_path, row.names = 1 , check.names = F, sep="\t")
  common_protein <- intersect(rownames(prediction), rownames(test_prot))
  prediction <- prediction[common_protein , colnames(test_prot) ]
  test_prot <- test_prot[common_protein , colnames(test_prot)]
  mat1 <- as.matrix(prediction)
  mat2 <- as.matrix(test_prot) 

  nrmse_vec <- c()
  for(i in 1:length(mat1[ ,1]) ) {
    temp <- hydroGOF::rmse(mat1[i,], mat2[i,],na.rm=T)
    nrmse_vec <- c(nrmse_vec , temp/(max(mat2[i,],na.rm=T)-min(mat2[i,],na.rm=T)))
  }
  names(nrmse_vec) <- rownames(mat1)
  return(mean(nrmse_vec))
}
#result_nrmse <- NRMSE_by_row("predictions.tsv", "prospective_ova_proteome_sort_common_gene_7061.txt")



###  SC3
################################## load pearson correlation function ##################################
correlation_by_row_ALL_OBSERVED <- function(pred_path, truth_path) {
  prediction <- read.csv( pred_path, row.names = 1 , check.names = F, sep="\t") 
  test_prot  <- read.csv( truth_path, row.names = 1 , check.names = F, sep="\t")
  #common_protein <- intersect(rownames(prediction), rownames(test_prot))
  prediction <- prediction[rownames(test_prot), colnames(test_prot)]
  test_prot <- test_prot[rownames(test_prot), colnames(test_prot)]
  
  mat1 <- as.matrix(prediction)
  mat2 <- as.matrix(test_prot) 
  
  mat2 <- mat2[complete.cases(mat2), ]
  mat1 <- mat1[rownames(mat2), ]
  
  corr_vec <- c()
  for(i in 1:length(mat1[ ,1]) ) {
    c <- rbind(mat1[i, ], mat2[i, ]) ; c <- c[ ,complete.cases(t(c))]
    temp <- cor.test(mat1[ i, ], mat2[ i , ])
    pcorr <- temp$estimate # pearson correlation
    if (is.na(pcorr)) {pcorr<-0}
    corr_vec <- c(corr_vec , pcorr)
  }
  names(corr_vec) <- rownames(mat1)
  return(mean(corr_vec))
}

########################################## load RMSE function #########################################
NRMSE_by_row_ALL_OBSERVED <- function(pred_path, truth_path)  {
  suppressPackageStartupMessages(library(hydroGOF))
  prediction <- read.csv( pred_path, row.names = 1 , check.names = F, sep="\t") 
  test_prot  <- read.csv( truth_path, row.names = 1 , check.names = F, sep="\t")
  common_protein <- intersect(rownames(prediction), rownames(test_prot))
  prediction <- prediction[common_protein , colnames(test_prot) ]
  test_prot <- test_prot[common_protein , colnames(test_prot)]
  mat1 <- as.matrix(prediction)
  mat2 <- as.matrix(test_prot) 
  
  mat2 <- mat2[complete.cases(mat2), ]
  mat1 <- mat1[rownames(mat2), ]
  
  nrmse_vec <- c()
  for(i in 1:length(mat1[ ,1]) ) {
    temp <- hydroGOF::rmse(mat1[i,], mat2[i,],na.rm=T)
    nrmse_vec <- c(nrmse_vec , temp/(max(mat2[i,],na.rm=T)-min(mat2[i,],na.rm=T)))
  }
  names(nrmse_vec) <- rownames(mat1)
  return(mean(nrmse_vec))
}


################################## load pearson correlation function ##################################
correlation_by_row_less30percMissing <- function(pred_path, truth_path) {
  prediction <- read.csv( pred_path, row.names = 1 , check.names = F, sep="\t") 
  test_prot  <- read.csv( truth_path, row.names = 1 , check.names = F, sep="\t")
  #common_protein <- intersect(rownames(prediction), rownames(test_prot))
  prediction <- prediction[rownames(test_prot), colnames(test_prot)]
  test_prot <- test_prot[rownames(test_prot), colnames(test_prot)]
  
  mat1 <- as.matrix(prediction)
  mat2 <- as.matrix(test_prot) 
  
  mat2 <- mat2[which(rowMeans(is.na(mat2)) < 0.3), ]
  mat1 <- mat1[rownames(mat2), ]
  
  corr_vec <- c()
  for(i in 1:length(mat1[ ,1]) ) { 
    c <- rbind(mat1[i, ], mat2[i, ]) ; c <- c[ ,complete.cases(t(c))]
    temp <- cor.test(mat1[ i, ], mat2[ i , ])
    pcorr <- temp$estimate # pearson correlation
    if (is.na(pcorr)) {pcorr<-0}
    corr_vec <- c(corr_vec , pcorr)
  }
  names(corr_vec) <- rownames(mat1)
  return(mean(corr_vec))
}

########################################## load RMSE function #########################################
NRMSE_by_row_ALL_less30percMissing <- function(pred_path, truth_path)  {
  suppressPackageStartupMessages(library(hydroGOF))
  prediction <- read.csv( pred_path, row.names = 1 , check.names = F, sep="\t") 
  test_prot  <- read.csv( truth_path, row.names = 1 , check.names = F, sep="\t")
  common_protein <- intersect(rownames(prediction), rownames(test_prot))
  prediction <- prediction[common_protein , colnames(test_prot) ]
  test_prot <- test_prot[common_protein , colnames(test_prot)]
  mat1 <- as.matrix(prediction)
  mat2 <- as.matrix(test_prot) 
  
  mat2 <- mat2[which(rowMeans(is.na(mat2)) < 0.3), ]
  mat1 <- mat1[rownames(mat2), ]
  
  nrmse_vec <- c()
  for(i in 1:length(mat1[ ,1]) ) {
    temp <- hydroGOF::rmse(mat1[i,], mat2[i,],na.rm=T)
    nrmse_vec <- c(nrmse_vec , temp/(max(mat2[i,],na.rm=T)-min(mat2[i,],na.rm=T)))
  }
  names(nrmse_vec) <- rownames(mat1)
  return(mean(nrmse_vec))
}




