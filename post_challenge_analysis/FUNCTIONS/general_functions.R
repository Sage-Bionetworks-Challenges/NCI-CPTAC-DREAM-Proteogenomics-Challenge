
########################## usual functiuons #########################
`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))
movetolast <- function(data, move) {  data[c(setdiff(names(data), move), move)] }
normalization_01 <- function(x) { x <- (x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)) ; return(x) }
## vapply(strsplit( vector , split="[.]"), "[", "", 1) take first element of a split 
## read_excel:  library(readxl) library(openxlsx)
########################  decompose tissues  ########################

retrieve_tissue <- function(name, tissue_label, features, mutation, response, tissue) { 

    subDir <- name
    dir.create(file.path(result_folder, subDir), showWarnings = FALSE)
    setwd(file.path(result_folder, subDir))
    t <- tissue_label[,1] [   which(tissue_label$tissue %in% tissue )  ] 
    t <- as.character(t)
    response_t <- response[rownames(response) %in%  t , ]
    
#     CRISPR <- read.csv("~/Documents/RWTH_Aachen/CRISPR/DATA/CRISPR", row.names=1, check.names = F)
#     common_cell <- intersect(rownames(CRISPR) , rownames(response_t) )
#     CRISPR <- CRISPR[common_cell, ]
#     response_t <- response_t[common_cell, ]   ;   response_t <- t(response_t)
#     CRISPR <- CRISPR[ ,colnames(CRISPR) %in% relevant_gene]
#     write.csv( CRISPR , paste( getwd() ,"/","CRISPR", sep = "") ) 
#     write.csv( response_t , paste( getwd() ,"/","RES_CRISPR", sep = "") ) 
    
    
    GEX_t <- features[ rownames(features) %in% t  , ]
    common_cell <- intersect(rownames(GEX_t) , rownames(response_t) )
    GEX_t <- GEX_t[common_cell, ]
    response_t <- response_t[common_cell, ]   ;   response_t <- t(response_t)
    write.csv( GEX_t , paste( getwd() ,"/","GEX", sep = "") ) 
    write.csv( response_t , paste( getwd() ,"/","RES", sep = "") ) 
    
    
#     SLC_ABC <- read.delim("~/Documents/RWTH_Aachen/SLC_TRANSPORTER/ADRIAN_SCRIPT/transporters_full_list.txt") ; SLC_ABC <- as.character(SLC_ABC$hgnc_symbol)
#     GEX_SLC_ABC <- GEX_t[ , colnames(GEX_t) %in% SLC_ABC ]
#     write.csv( GEX_SLC_ABC , paste( getwd() ,"/","GEX_SLC_ABC", sep = "") ) 
    
    # ECFP4 <- read.csv("/Users/miyang/Documents/RWTH_Aachen/SANGER_DATA/CHEMOINFORMATICS/Drugs_MACAU_ECFP4_names_edited.csv", sep=";",row.names = 1)
    # ECFP4 <- ECFP4[order(rownames(ECFP4)) , ]
    # common_drug <- intersect(rownames(ECFP4),rownames(response_t))
    # response_t_ECFP4 <- response_t[rownames(response_t) %in% common_drug,  ]
    # response_t_ECFP4 <- response_t_ECFP4[order(rownames(response_t_ECFP4)) , ]
    # ECFP4 <- ECFP4[common_drug , ]
    # response_t_ECFP4 <- response_t_ECFP4[unique(rownames(response_t_ECFP4)) , ]
    # write.csv( response_t_ECFP4 , paste( getwd() ,"/","RES_ECFP4", sep = "") ) 
    # write.csv( ECFP4 , paste( getwd() ,"/","ECFP4", sep = "") ) 
    
    
    # IC50_target_Leiden <- read.csv("/Users/miyang/Documents/RWTH_Aachen/macau_work_dir/macau_test_sanger/DATA/IC50_target_Leiden", row.names = 1, check.names = F)
    # IC50_target_Leiden <- IC50_target_Leiden[ ,common_cell]
    # write.csv( IC50_target_Leiden , paste( getwd() ,"/","RES_target_Leiden", sep = "") ) 
    
    model <- read.csv("/Users/miyang/Documents/RWTH_Aachen/PROGENY/model_14PW.csv", row.names = 1)
    common_gene <- intersect(colnames(GEX_t) , rownames(model))
    GEX <- as.matrix(GEX_t[ ,common_gene ]) ; model <- as.matrix(model[ common_gene, ])
    progeny14 <- GEX %*% model
    progeny14 <- scale(progeny14) ## scale every column
    write.csv( progeny14 , paste( getwd() ,"/","progeny14", sep = "") ) 
    
    load("~/Documents/RWTH_Aachen/PROGENY/model.RData")
    common_gene <- intersect(colnames(GEX_t) , rownames(model))
    G <- as.matrix(GEX_t[ ,common_gene ]) ; model <- as.matrix(model[ common_gene, ])
    progeny11 <- G %*% model
    progeny11 <- scale(progeny11) ## scale every column
    write.csv( progeny11 , paste( getwd() ,"/","progeny11", sep = "") ) 
    # 
    # progeny_old11_new3 <- cbind(progeny11, progeny14[ ,c("Androgen","Estrogen","WNT")])
    # write.csv(progeny_old11_new3, paste( getwd() ,"/","progeny_old11_new3", sep = "") ) 
    
    # GEX_by_gene <- c()
    # for(i in 1:length(colnames(GEX_t))) {  GEX_by_gene <- cbind(GEX_by_gene,  (GEX_t[ ,i]-min(GEX_t[ ,i]))/(max(GEX_t[ ,i])-min(GEX_t[ ,i]))  )  }
    # rownames(GEX_by_gene) <- rownames(GEX_t) ; colnames(GEX_by_gene) <- colnames(GEX_t) 
    # write.csv( GEX_by_gene , paste( getwd() ,"/","GEX_by_gene", sep = "") ) 
    
    # x <- scale(GEX_t)
    # write.csv( x , paste( getwd() ,"/","GEX_scale_by_gene", sep = "") ) 
    # 
    # y <- t(scale(t(GEX_t)))
    # write.csv( y , paste( getwd() ,"/","GEX_scale_by_sample", sep = "") ) 
    # 
    # SNP_CNV <- mutation[ common_cell ,   ] ; SNP_CNV <- SNP_CNV[ , - which(colSums(SNP_CNV) == 0 ) ] 
    # write.csv( SNP_CNV , paste( getwd() ,"/","SNP_CNV", sep = "") ) 
    # SNP_CNV_to_Remove <- colnames(SNP_CNV)[which(colSums(SNP_CNV) == 1) ]
    # write.csv( SNP_CNV_to_Remove , paste( getwd() ,"/","SNP_CNV_to_Remove", sep = "") ) 
    # 
    # NES_GDSC_GSVA <- NES_GDSC_GSVA[common_cell, ] ; TF_GSVA <- NES_GDSC_GSVA
    # write.csv( TF_GSVA , paste( getwd() ,"/","TF_GSVA", sep = "") ) 
    # 
    # NES_GDSC_VIPER <- NES_GDSC_VIPER[common, ] ; TF_VIPER <- NES_GDSC_VIPER
    # write.csv( TF_VIPER , paste( getwd() ,"/","TF_VIPER", sep = "") ) 
    # 
    # progeny11 <- read.csv( paste( getwd() ,"/","progeny11", sep = ""), row.names = 1 )
    # colnames(progeny11) <- paste("Progeny_",colnames(progeny11),sep = "")
    # colnames(TF_GSVA) <- paste("TF_",colnames(TF_GSVA),sep = "")
    # GEX_progeny_TF_GSVA <- cbind(progeny11, TF_GSVA)
    # write.csv( GEX_progeny_TF_GSVA , paste( getwd() ,"/","progeny11_TF_GSVA", sep = "") ) 
    # 
    # colnames(progeny14) <- paste("Progeny_",colnames(progeny14),sep = "")
    # colnames(TF_GSVA) <- paste("TF_",colnames(TF_GSVA),sep = "")
    # progeny14_TF_GSVA <- cbind(progeny14, TF_GSVA)
    # write.csv( progeny14_TF_GSVA , paste( getwd() ,"/","progeny14_TF_GSVA", sep = "") ) 

}


retrieve_tissue_CTRP <- function(name, tissue_label, features, mutation, response, tissue) { 
  
  subDir <- name
  dir.create(file.path(result_folder, subDir), showWarnings = FALSE)
  setwd(file.path(result_folder, subDir))
  t <- tissue_label[,1] [   which(tissue_label$tissue %in% tissue )  ] 
  t <- as.character(t)
  response_t <- response[rownames(response) %in%  t , ]
  
  GEX_t <- features[ rownames(features) %in% t  , ]
  common_cell <- intersect(rownames(GEX_t) , rownames(response_t) )
  GEX_t <- GEX_t[common_cell, ]
  response_t <- response_t[common_cell, ]   ;   response_t <- t(response_t)
  write.csv( GEX_t , paste( getwd() ,"/","GEX", sep = "") ) 
  write.csv( response_t , paste( getwd() ,"/","RES", sep = "") ) 

  load("~/Documents/RWTH_Aachen/PROGENY/model.RData")
  common_gene <- intersect(colnames(GEX_t) , rownames(model))
  G <- as.matrix(GEX_t[ ,common_gene ]) ; model <- as.matrix(model[ common_gene, ])
  progeny11 <- G %*% model
  progeny11 <- scale(progeny11) ## scale every column
  write.csv( progeny11 , paste( getwd() ,"/","progeny11", sep = "") ) 
  
  CNV <- features[ rownames(features) %in% t  , ]
  common_cell <- intersect(rownames(CNV) , rownames(response_t) )
  CNV <- CNV[common_cell, ]
  response_t <- response_t[common_cell, ]   ;   response_t <- t(response_t)
  write.csv( CNV , paste( getwd() ,"/","CNV", sep = "") ) 
  write.csv( response_t , paste( getwd() ,"/","RES_CNV", sep = "") ) 

}


retrieve_tissue_CCLE <- function(name, tissue_label, features, mutation, response, tissue) { 
  
  subDir <- name
  dir.create(file.path(result_folder, subDir), showWarnings = FALSE)
  setwd(file.path(result_folder, subDir))
  t <- tissue_label[,1] [   which(tissue_label$tissue %in% tissue )  ] 
  t <- as.character(t)
  response_t <- response[rownames(response) %in%  t , ]
  
  GEX_t <- features[ rownames(features) %in% t  , ]
  common_cell <- intersect(rownames(GEX_t) , rownames(response_t) )
  GEX_t <- GEX_t[common_cell, ]
  response_t <- response_t[common_cell, ]   ;   response_t <- t(response_t)
  write.csv( GEX_t , paste( getwd() ,"/","GEX", sep = "") ) 
  write.csv( response_t , paste( getwd() ,"/","RES", sep = "") ) 

  load("~/Documents/RWTH_Aachen/PROGENY/model.RData")
  common_gene <- intersect(colnames(GEX_t) , rownames(model))
  G <- as.matrix(GEX_t[ ,common_gene ]) ; model <- as.matrix(model[ common_gene, ])
  progeny11 <- G %*% model
  progeny11 <- scale(progeny11) ## scale every column
  write.csv( progeny11 , paste( getwd() ,"/","progeny11", sep = "") ) 
  
#   CNV <- features[ rownames(features) %in% t  , ]
#   common_cell <- intersect(rownames(CNV) , rownames(response_t) )
#   CNV <- CNV[common_cell, ]
#   response_t <- response_t[common_cell, ]   ;   response_t <- t(response_t)
#   write.csv( CNV , paste( getwd() ,"/","CNV", sep = "") ) 
#   write.csv( response_t , paste( getwd() ,"/","RES_CNV", sep = "") ) 
  
}

retrieve_tissue_NOVARTIS_PDX <- function(name, tissue_label, features, mutation, response, tissue) { 
  
  subDir <- name
  dir.create(file.path(result_folder, subDir), showWarnings = FALSE)
  setwd(file.path(result_folder, subDir))
  t <- tissue_label[,1] [   which(tissue_label$HISTOLOGY %in% tissue )  ] 
  t <- as.character(t)
  response_t <- response[rownames(response) %in%  t , ]
  
  GEX_t <- features[ rownames(features) %in% t  , ]
  common_cell <- intersect(rownames(GEX_t) , rownames(response_t) )
  GEX_t <- GEX_t[common_cell, ]
  response_t <- response_t[common_cell, ]   ;   response_t <- t(response_t)
  write.csv( GEX_t , paste( getwd() ,"/","GEX", sep = "") ) 
  write.csv( response_t , paste( getwd() ,"/","RES", sep = "") ) 
  
  load("~/Documents/RWTH_Aachen/PROGENY/model.RData")
  common_gene <- intersect(colnames(GEX_t) , rownames(model))
  G <- as.matrix(GEX_t[ ,common_gene ]) ; model <- as.matrix(model[ common_gene, ])
  progeny11 <- G %*% model
  progeny11 <- scale(progeny11) ## scale every column
  write.csv( progeny11 , paste( getwd() ,"/","progeny11", sep = "") ) 
  
  #   CNV <- features[ rownames(features) %in% t  , ]
  #   common_cell <- intersect(rownames(CNV) , rownames(response_t) )
  #   CNV <- CNV[common_cell, ]
  #   response_t <- response_t[common_cell, ]   ;   response_t <- t(response_t)
  #   write.csv( CNV , paste( getwd() ,"/","CNV", sep = "") ) 
  #   write.csv( response_t , paste( getwd() ,"/","RES_CNV", sep = "") ) 
  
}



retrieve_tissue_TCGA <- function(name, tissue_label, features, tissue) { 
  subDir <- name
  dir.create(file.path(result_folder, subDir), showWarnings = FALSE)
  setwd(file.path(result_folder, subDir))
  t <- tissue_label[,1] [   which(tissue_label$tissue %in% tissue )  ] 
  t <- as.character(t)
  GEX <- features[ rownames(features) %in% t  , ]
  GEX <- GEX[unique(rownames(GEX)) , ]
  
  common_cell <- intersect(rownames(GEX), rownames(PANCANCER_CLI))
  GEX <- GEX[ common_cell , ]
  save(GEX, file=paste( getwd() ,"/","GEX.Rdata", sep = ""))
  
  library(xCell) ; x <- t(GEX) ; xCell <- xCellAnalysis(x) ; xCell <- t(xCell)
  write.csv( xCell , paste( getwd() ,"/","xCell", sep = "") ) 
  
#   GEX_by_gene <- c()
#   for(i in 1:length(colnames(GEX))) {  GEX_by_gene <- cbind(GEX_by_gene,  (GEX[ ,i]-min(GEX[ ,i]))/(max(GEX[ ,i])-min(GEX[ ,i]))  )  }
#   rownames(GEX_by_gene) <- rownames(GEX) ; colnames(GEX_by_gene) <- colnames(GEX) 
#   save(GEX_by_gene, file=paste( getwd() ,"/","GEX_by_gene.Rdata", sep = ""))

  load("~/Documents/RWTH_Aachen/PROGENY/model.RData")
  common_gene <- intersect(colnames(GEX) , rownames(model))
  G <- as.matrix(GEX[ ,common_gene ]) ; model <- as.matrix(model[ common_gene, ])
  progeny11 <- G %*% model
  progeny11 <- scale(progeny11) ## scale every column
  write.csv( progeny11 , paste( getwd() ,"/","progeny11", sep = "") ) 
  
  x <- NES_TCGA_GSVA[common_cell, ] ; TF_GSVA <- x
  write.csv( TF_GSVA , paste( getwd() ,"/","TF_GSVA", sep = "") ) 
  
  x <- NES_TCGA_VIPER[common_cell, ] ; TF_VIPER <- x
  write.csv( TF_VIPER , paste( getwd() ,"/","TF_VIPER", sep = "") ) 
  
  progeny <- read.csv( paste( getwd() ,"/","progeny11", sep = ""), row.names = 1 )
  colnames(progeny) <- paste("Progeny_",colnames(progeny),sep = "")
  colnames(TF_GSVA) <- paste("TF_",colnames(TF_GSVA),sep = "")
  progeny_TF_GSVA <- cbind(progeny, TF_GSVA)
  write.csv( progeny_TF_GSVA , paste( getwd() ,"/","progeny11_TF_GSVA", sep = "") ) 
  
  model <- read.csv("/Users/miyang/Documents/RWTH_Aachen/PROGENY/model_14PW.csv", row.names = 1)
  common_gene <- intersect(colnames(GEX) , rownames(model))
  G <- as.matrix(GEX[ ,common_gene ]) ; model <- as.matrix(model[ common_gene, ])
  progeny14 <- G %*% model
  progeny14 <- scale(progeny14) ## scale every column
  write.csv( progeny14 , paste( getwd() ,"/","progeny14", sep = "") ) 
  
  colnames(progeny14) <- paste("Progeny_",colnames(progeny14),sep = "")
  colnames(TF_GSVA) <- paste("TF_",colnames(TF_GSVA),sep = "")
  progeny14_TF_GSVA <- cbind(progeny14, TF_GSVA)
  write.csv( progeny14_TF_GSVA , paste( getwd() ,"/","progeny14_TF_GSVA", sep = "") ) 
  
  x <- PANCANCER_CLI[ common_cell , c("DAYS_TO_LAST_FOLLOWUP", "OS_STATUS")]
  rank <- guan_rank(x) ; write.csv( rank , paste( getwd() ,"/","rank", sep = "") ) 
   
}

load_folder <- function (workdir ) {
  setwd(workdir)
  files <- list.files(path=getwd() , pattern="") # 
  n <- 0
  imp <- list()
  for(file in files)
  {
    n <- n + 1
    imp[[n]] <- read.csv(paste(workdir, "/",files[n], sep = ""), row.names = 1 , check.names = T)
  }
  names(imp) <- files
  return(imp)
}

prepare_for_macau <- function(mat) {
  colnames(mat) <- 0:(length(colnames(mat))-1) ; rownames(mat) <- 0:(length(rownames(mat))-1)
  return(mat)
}

create_Tensor <- function( list_of_mat ) { 
  
  Rows <- c() ; Cols <- c()
  for(i in 1:length(list_of_mat)) {
    Rows <- c(Rows, rownames(list_of_mat[[i]]))
    Cols <- c(Cols, colnames(list_of_mat[[i]]))  
  }
  Rows <- unique(Rows) ; Cols <- unique(Cols)
  Rows <- Rows[order(Rows)] ; Cols <- Cols[order(Cols)]
  
  Tensor <- c()
  for(mat in 1:length(list_of_mat)) {
  df <- list_of_mat[[mat]]
  model_mat <- matrix(nrow=length(Rows),ncol=length(Cols))
  rownames(model_mat) <- Rows ; colnames(model_mat) <- Cols
  coord <- which(is.na(df)==F, arr.ind = T)
  for(i in 1:length(coord[,1])) {
  model_mat[rownames(df)[coord[i,1]],colnames(df)[coord[i,2]]] <- df[coord[i,1],coord[i,2]]
  }
  
  x <- data.frame(rowID=as.vector(row(model_mat)), colID=as.vector(col(model_mat)), mat=as.vector(model_mat)) 
  x$rowID <- x$rowID - 1 ;  x$colID <- x$colID - 1
  x$depth <- rep(mat-1, length(x[,1]))
  x <- x[c("rowID", "colID", "depth", "mat")]
  
  Tensor <- rbind(Tensor, x)
  }
  
  return(list(Tensor,model_mat))
}

storedata <- function(x, y) {
  setwd(y)
  write.csv(x, file = deparse(substitute(x)))
}

stopQuietly <- function(...) {
  blankMsg <- sprintf("\r%s\r", paste(rep(" ", getOption("width")-1L), collapse=" ")); stop(simpleError(blankMsg));
} # stopQuietly()

ClosestMatch2 = function(string, stringVector, maxDist=1){ 
  library(dplyr)
  library(stringdist)
  stringVector[amatch(string, stringVector, maxDist=maxDist)] 
}

ID_to_celllines <- function(x) {
  
  MASTER_LIST <- as.matrix(MASTER_LIST)
  MASTER_LIST[753,2] <- "PC-3-[JPC-3]"
  M <- MASTER_LIST[ ,1:2]
  names <- c()
  for (i in 1:length(rownames(x))) {
    if(length(grep(rownames(x)[i], M[,1])) == 1){
      names <- c(names, as.character(M[,2][grep(rownames(x)[i], M[,1])]) ) } else {names <- c(names, rownames(x)[i])}
  }
  rownames(x) <- names
  return (x)
}

celllines_to_ID <- function(x) {
  
  MASTER_LIST <- as.matrix(MASTER_LIST)
  MASTER_LIST[753,2] <- as.character(MASTER_LIST[753,2])
  MASTER_LIST[753,2] <- "PC-3-[JPC-3]"
  M <- MASTER_LIST[ ,1:2]
  names <- c()
  for (i in 1:length(rownames(x))) {
    if(length(grep(rownames(x)[i], M[,2])) == 1){
      names <- c(names, M[,1][grep(rownames(x)[i], M[,2])] ) } else {names <- c(names, rownames(x)[i])
      }
  }
  names <- trimws(names)
  rownames(x) <- names
  return (x)
}

ID_to_celllines_RNAseq <- function(x) {
  CLI <- CLI[ ,c(16, 3, 11)]
  M <- CLI[ ,1:2]
  names <- c()
  for (i in 1:length(rownames(x))) {
    if(length(grep(rownames(x)[i], M[,1])) == 1){
      names <- c(names, as.character(M[,2][grep(rownames(x)[i], M[,1])]) ) } else {names <- c(names, rownames(x)[i])}
  }
  rownames(x) <- names
  return (x)
}

Hugo_to_EntrezID <- function(gene) {
  library(biomaRt)
  ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
  mapTab <- getBM(attributes = c("hgnc_symbol", "entrezgene"), filters = "hgnc_symbol", values = gene, mart = ensembl, uniqueRows=FALSE)
  dupRows <- union(which(duplicated(mapTab[,1])), which(duplicated(mapTab[,2])))
  entrezIds <- mapTab[-dupRows, 2]
  names(entrezIds) <- mapTab[-dupRows, 1]
  entrezIds <- entrezIds[!is.na(entrezIds)]
  return(entrezIds)
}

EntrezID_to_Hugo <- function(gene) {
  library(biomaRt)
  ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
  mapTab <- getBM(attributes = c("entrezgene", "hgnc_symbol"), filters = "entrezgene", values = gene, mart = ensembl, uniqueRows=FALSE)
  dupRows <- union(which(duplicated(mapTab[,1])), which(duplicated(mapTab[,2])))
  hgnc <- mapTab[-dupRows, 2]
  names(hgnc) <- mapTab[-dupRows, 1]
  hgnc <- hgnc[!is.na(hgnc)]
  return(hgnc)
}

ENS_to_Hugo <- function(gene) {
  library(biomaRt) ## mart = useMart('ensembl') ; listDatasets(mart) ; BiocInstaller::biocLite('grimbough/biomaRt')
  ensembl <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl") ## "hsapiens_gene_ensembl"
  mapTab <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = gene, mart = ensembl, uniqueRows=FALSE)
  dupRows <- union(which(duplicated(mapTab[,1])), which(duplicated(mapTab[,2])))
  hgnc <- mapTab[-dupRows, 2]
  names(hgnc) <- mapTab[-dupRows, 1]
  hgnc <- hgnc[!is.na(hgnc)]
  return(hgnc)
}

 
convert_drugID <- function (x) {
  ii <- c()
  for (i in 1:length(rownames(x))) {
    for (j in 1:length(DRUG_ANALYSIS_SET$DRUG_ID)){
      if (rownames(x)[i] == DRUG_ANALYSIS_SET$DRUG_ID[j]) {
        ii <- c(ii, as.character(DRUG_ANALYSIS_SET$DRUG_NAME)[j])
      }
    }
  }
  rownames(x) <- ii
  return (x)
}

processIteration <- function(drug_temp, n) {
  ## n: number of iterations ## 
  start <- 1-n
  end <- 0
  pcorr <- c() ## pearson correlation
  r2 <- c() ## R2
  sdev <- c() ## standard deviation of pearson correlation
  for (i in 1:(length(drug_temp[ ,1])/n)){
    start <- start + n
    end <- end + n
    temp_corr <- drug_temp[start:end, 1:2]
    temp_corr <- temp_corr[complete.cases(temp_corr[2]), ]
    temp_corr <- temp_corr[ ,1]
    temp_r2 <- temp_corr^2
    pcorr <- c(pcorr, mean(temp_corr, na.rm= T))
    r2 <- c(r2, mean(temp_r2, na.rm= T))
    sdev <- c(sdev,  sd(temp_corr, na.rm= T) )
  }
  pcr2 <- cbind(pcorr, r2, sdev)
  return ( pcr2)
}


create_features <- function(mat, column) {
  
  features <- mat
  col <- which(colnames(features) == column)
  
  features[ ,col] <- as.character(features[ ,col])
  features[is.na(features[ ,col]) , col ] <- "unknown"
  target <- as.character( features[ ,col] ) 
  target <- gsub("[[:blank:]]", "", target)
  target = strsplit(target, ",")  
  target <- unlist(target, recursive = TRUE, use.names = TRUE) ; target <- unique(target)
  
  m <- matrix(nrow= length(target) , ncol = length(features[,1]) )
  rownames(m) <- target  ;  colnames(m) <- rownames(features)
  
  for(i in 1:length(colnames(m))) {
    a <- as.character( features[ ,col][i] ) 
    a <- gsub("[[:blank:]]", "", a)
    a = strsplit(a, ",")  ; a <- a[[1]]
    
    library(stringdist)
    ClosestMatch2 = function(string, stringVector){  stringVector[amatch(string, stringVector, maxDist=Inf)] }
    match <- ClosestMatch2 (a,  rownames(m))
    
    m[ which(rownames(m) %in% match ) , i ] <- 1
  }
  m[is.na(m)] <- 0
  # print(rowSums(m)) ; print(colSums(m))
  m <-t(m)
  if(length(which(colnames(m) == "unknown" )==1)) {
    m <- m[ , -which(colnames(m) == "unknown" ) ]
  }
  return(m)
}

create_features_dummy <- function(features, column) {
  
  col <- which(colnames(features) == column)
  
  features[ ,col] <- as.character(features[ ,col])
  features[is.na(features[ ,col]) , col ] <- "unknown"
  target <- as.character( features[ ,col] ) 
  target <- gsub("[[:blank:]]", "", target)
  target = strsplit(target, ",")  
  target <- unlist(target, recursive = TRUE, use.names = TRUE) ; target <- unique(target)
  
  m <- matrix(nrow= length(target) , ncol = length(features[,1]) )
  rownames(m) <- target  ;  colnames(m) <- rownames(features)
  
  for(i in 1:length(colnames(m))) {
    a <- as.character( features[ ,col][i] ) 
    a <- gsub("[[:blank:]]", "", a)
    a = strsplit(a, ",")  ; a <- a[[1]]
    
    library(stringdist)
    ClosestMatch2 = function(string, stringVector){  stringVector[amatch(string, stringVector, maxDist=Inf)] }
    match <- ClosestMatch2 (a,  rownames(m))
    
    m[ which(rownames(m) %in% match ) , i ] <- 1
  }
  m[is.na(m)] <- 0
  # print(rowSums(m)) ; print(colSums(m))
  m <-t(m)
  if(length(which(colnames(m) == "unknown" )==1)) {
  m <- m[ , -which(colnames(m) == "unknown" ) ]
  }
  reference <- colnames(m)[1] ; m <- m[ ,-1]
  return(list(m,reference))
}


create_features_dummy_MANY <- function(features, range=3:8) {
  
  col <- range
  
  mat_combined <- matrix(data=NA,nrow=length(rownames(features)), ncol = 1)
  for(i in 1:length(col)) {   #   col[1]
    features[ ,col[i]] <- as.character(features[ ,col[i]])
    features[is.na(features[ ,col[i]]) , col[i] ] <- "unknown"
    target <- as.character( features[ ,col[i]] ) 
    target <- gsub("[[:blank:]]", "", target)
    target = strsplit(target, ",")  
    target <- unlist(target, recursive = TRUE, use.names = TRUE) ; target <- unique(target)
    
    m <- matrix(nrow= length(target) , ncol = length(features[,1]) )
    rownames(m) <- target  ;  colnames(m) <- rownames(features)
    
    for(j in 1:length(colnames(m))) {
      a <- as.character( features[ ,col[i]][j] ) 
      a <- gsub("[[:blank:]]", "", a)
      a = strsplit(a, ",")  ; a <- a[[1]]
      
      library(stringdist)
      ClosestMatch2 = function(string, stringVector){  stringVector[amatch(string, stringVector, maxDist=Inf)] }
      match <- ClosestMatch2 (a,  rownames(m))
      
      m[ which(rownames(m) %in% match ) , j ] <- 1
    }
    m[is.na(m)] <- 0
    # print(rowSums(m)) ; print(colSums(m))
    m <-t(m)
    if(length(which(colnames(m) == "unknown" )==1)) {
      m <- m[ , -which(colnames(m) == "unknown" ) ]
    }
    reference <- colnames(m)[1] ; m <- m[ ,-1]
    mat_combined <- cbind(mat_combined, m)
  }
  
  mat_combined <- mat_combined[ , -1]
  return(mat_combined)
}



corr_by_row <- function(mat1, mat2, method="pearson")  {
  mat1 <- as.matrix(mat1)
  mat2 <- as.matrix(mat2) 

  corr_vec <- c()
  for(i in 1:length(mat1[,1]) ) {
    
    c <- rbind(mat1[i, ], mat2[i, ]) ; c <- c[ ,complete.cases(c)]
    if(length(which(apply(c, 1, var) == 0)) > 0) { corr_vec <- c(corr_vec , 0 ) } else 
      {
      temp <- cor.test(mat1[i, ], mat2[i, ],method=method)
      pcorr <- temp$estimate # pearson correlation
      corr_vec <- c(corr_vec , pcorr)
      }
  }
  boxplot(corr_vec ) ; print(paste("mean: ", mean(corr_vec, na.rm = T), sep = "") )   ; print(paste("median: ", median(corr_vec, na.rm = T), sep = "")  ) 
  names(corr_vec) <- rownames(mat1)
  return(corr_vec)
}


rmse_by_row <- function(mat1, mat2)  {
  mat1 <- as.matrix(mat1)
  mat2 <- as.matrix(mat2) 
  
  rmse_vec <- c()
  for(i in 1:length(mat1[,1]) ) {
    rmse_vec <- c(rmse_vec , rmse(mat1[i, ], mat2[i, ], na.rm=TRUE))
  }
  boxplot(rmse_vec ) ; print(paste("mean: ", mean(rmse_vec, na.rm = T), sep = "") )   ; print(paste("median: ", median(rmse_vec, na.rm = T), sep = "")  ) 
  names(rmse_vec) <- rownames(mat1)
  return(rmse_vec)
}

common_full <- function(MAT_LIST)  {
  for (iteration in 1:2) {
  row_list <- list()
  col_list <- list()
  for(i in 1:length(MAT_LIST)) {
    row_list[[i]] <- rownames(MAT_LIST[[i]])
    col_list[[i]] <- colnames(MAT_LIST[[i]])
  }
  common_row <- Reduce(intersect, row_list )
  common_col <- Reduce(intersect, col_list )
  for(i in 1:length(MAT_LIST)) {
    MAT_LIST[[i]] <-MAT_LIST[[i]][common_row, common_col]
    MAT_LIST[[i]] <- MAT_LIST[[i]][rowSums(is.na(MAT_LIST[[i]])) < length(MAT_LIST[[i]][1,]) - 2 , ] # remove rows with too many NAs
   }
  }
  return(MAT_LIST)
}

common_row <- function(MAT_LIST)  {
  for (iteration in 1:2) {
    row_list <- list()
    for(i in 1:length(MAT_LIST)) {
      row_list[[i]] <- rownames(MAT_LIST[[i]])
    }
    common_row <- Reduce(intersect, row_list )
    for(i in 1:length(MAT_LIST)) {
      MAT_LIST[[i]] <-MAT_LIST[[i]][common_row, ]
      MAT_LIST[[i]] <- MAT_LIST[[i]][rowSums(is.na(MAT_LIST[[i]])) < length(MAT_LIST[[i]][1,]) - 2 , ] # remove rows with too many NAs
    }
  }
  return(MAT_LIST)
}
common_col <- function(MAT_LIST)  {
  for (iteration in 1:2) {
    col_list <- list()
    for(i in 1:length(MAT_LIST)) {
      col_list[[i]] <- colnames(MAT_LIST[[i]])
    }
    common_col <- Reduce(intersect, col_list )
    for(i in 1:length(MAT_LIST)) {
      MAT_LIST[[i]] <-MAT_LIST[[i]][ , common_col ]
      MAT_LIST[[i]] <- MAT_LIST[[i]][rowSums(is.na(MAT_LIST[[i]])) < length(MAT_LIST[[i]][1,]) - 2 , ] # remove rows with too many NAs
    }
  }
  return(MAT_LIST)
}

### if 1 NA, then mean is NA
aggregate_mean <- function(x_sep) {
  
  dup <- x_sep[,1][duplicated(x_sep[,1])] ; xtot_dup <- x_sep[x_sep[,1] %in% dup, ] ; xtot_dup[,1] <- as.character(xtot_dup[,1])
  mean <- c()
  for(i in 1:length(table(xtot_dup[,1])) ) {
    v <-  xtot_dup [xtot_dup[,1] %in% names(table(xtot_dup[,1])[i]) , ]
    v <- as.matrix(v[ , 2:length(v[1,]) ])
    mean <- rbind(mean , colMeans(v) ) 
  }
  rownames(mean) <- names(table(xtot_dup[,1]))
  non_dup <- x_sep[,1][!duplicated(x_sep[,1])] ; non_dup <- setdiff(non_dup, dup) ; xtot_non_dup <- x_sep[x_sep[,1]  %in% non_dup, ] 
  xtot_non_dup <-xtot_non_dup[,-1] ; xtot_non_dup <- as.matrix(xtot_non_dup) ; rownames(xtot_non_dup) <- non_dup
  x <- rbind(mean , xtot_non_dup) ; x <- x[order(rownames(x)), ]
  x <- x[rowSums(is.na(x))!=length(x[1,]), ]
  return(x)
}

features_variance <- function( features, RES, RES_number, top) {
  common <- intersect(rownames(features), colnames(RES))
  fea <- features[common, ] ; res <- RES[ , common ]
  key_predictor <- rownames(RES)[RES_number]
  
  col_var <- apply(fea, 2, var) 
  fea <- fea[ , which(col_var > quantile(col_var, 1-top)) ]
  if(length(which(colnames(fea) == key_predictor)) == 0) {
  final_features <- cbind( features[common,key_predictor] , fea ) ; colnames(final_features)[1] <- key_predictor
  } else { col_idx <- grep(key_predictor, colnames(fea)) ; 
           final_features <- fea[, c(col_idx, (1:ncol(fea))[-col_idx])] }
  return(final_features)
}

get_result <- function(vec, result_folder) {
  result <- list()
  for (i in 1:length(vec)) {
    v <- read.csv(paste(path, "/COREAD/", result_folder, "/", vec[i], sep=""))
    result[[i]] <- v[ ,2]
  }
  return(result)
}


make_combo <- function(vector, element=2) {
  combi <- t( combn(x=vector, m=element, FUN = NULL, simplify = TRUE ) ) 
  list_of_combo <- c()
  for(i in 1:length(combi[,1])) {
    pair <- c( combi[i,1] , combi[i,2] ) ; pair <- pair[order(pair)] ;  pair <- paste(pair[1],"_",pair[2],sep="")
    list_of_combo <- c(list_of_combo, pair) 
    }
  list_of_combo <- unique(list_of_combo)
  return(list_of_combo)
}

quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

print_drug_from_gene <- function(w, gene, rank=5) {
  drug <- c()
  for(i in 1:length(w)) { if(length(which(rownames(w[[i]])[1:rank]==gene)) == 1) { drug <- c(drug, i)}  }
  drug <- unique(names(w)[drug])
  return(drug)
}

print_rank_weight <- function (predictor, weight_object, subset) {
  predictor_for_each <- c()
  for(i in 1:length(weight_object)) { predictor_for_each <- c(predictor_for_each, which(names(weight_object[[i]]) == predictor ) ) }
  names(predictor_for_each) <- names(weight_object)
  
  value_for_each <- c()
  for(i in 1:length(weight_object)) { value_for_each <- c(value_for_each, weight_object[[i]][predictor] ) }
  names(value_for_each) <- names(weight_object) 
  
  L <- list(predictor_for_each[subset], value_for_each[subset]) ; names(L) <- c(paste0(predictor," rank"),paste0(predictor," weight"))
  return(L)
}


find_pubmed_association <- function (vector_row, vector_column, constant_term="") {
library(RISmed)
m <- matrix(nrow = length(vector_row), ncol = length(vector_column))
rownames(m) <- vector_row ; colnames(m) <- vector_column
for(i in 1:length(vector_row)) {
  for(j in 1:length(vector_column)) {
    tryCatch({
      query <- paste("(",vector_row[i]," AND ",vector_column[j]," AND ",constant_term,")", sep = "" )
      mindate=2008 ; maxdate=2018
      ngs_search <- EUtilsSummary(query, type="esearch",db = "pubmed",mindate=mindate, maxdate=maxdate, retmax=500)
      m[i,j] <- QueryCount(ngs_search)
    }, error=function(e){})
  }
}
return(m)
}


change_name <- function( to_change , to_change_name , reference_name ) {
  to_change[which(to_change == to_change_name )] <- reference_name
  return(to_change)
}

binarize_by_column <- function(mat, th, binary=F) {
  for(i in 1:length(colnames(mat)) ){ # i=1
    limit <- quantile(mat[ ,i], th)
    mat[ ,i][mat[ ,i]<limit]  <- "low" 
    mat[ ,i][mat[ ,i] %not in% "low"] <- "high" 
  }
  if(binary==T) {
    mat[mat=="high"] <- 1
    mat[mat=="low" ] <- 0
  }
  return(mat)
}


