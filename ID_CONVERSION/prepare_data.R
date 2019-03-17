path <- "~/Documents/RWTH_Aachen"
source("~/Documents/RWTH_Aachen/COREAD/COREAD_Function.R")

############################################# CONVERSION TABLE ###########################################
setwd("~/Documents/RWTH_Aachen/CONVERSION")
prot <- read.csv(paste(path, "/COREAD/FEATURES/prot", sep=""), row.names=1)

uniprot <- read.delim("all_uniprot_human_accessions-MWs.txt", header=TRUE)
uniprot_transformed <- uniprot [, c(1,2,4)]
a <- c() ; b <- c() ; c <- c()
for (i in 1:length(uniprot_transformed[,2])) {
  a <- c(a, as.character(uniprot_transformed[i,1]))
  b <- c(b, as.character(uniprot_transformed[i,2]))
  c <- c(c, as.character(uniprot_transformed[i,3]))
}
uniprot_transformed <- cbind(a,b,c) ; colnames(uniprot_transformed) <- colnames(uniprot)[c(1,2,4)]
for (i in 1:length(uniprot_transformed[,2])) {
  uniprot_transformed[i,2] <- substr(uniprot_transformed[i,2], 1, nchar(uniprot_transformed[i,2])-6) 
}
storedata(uniprot_transformed,"~/Documents/RWTH_Aachen/CONVERSION") ####### an actionable version of Uniprot

####################################### PROTEINS ID to HUGO gene name #####################################
uniprot_transformed <- read.csv("~/Documents/RWTH_Aachen/CONVERSION/uniprot_transformed", row.names=2) ; uniprot_transformed <- uniprot_transformed [,-1]

prot_table <- uniprot_transformed[ colnames(prot) , ]
storedata(prot_table,"~/Documents/RWTH_Aachen/CONVERSION") ####### COREAD protein ID and names

colnames(prot) <- prot_table$Entry.name
storedata(prot,"~/Documents/RWTH_Aachen/CONVERSION") ####### COREAD protein names

######################################### corresponding gene in GEX #####################################
GEX_coread <- read.csv("~/Documents/RWTH_Aachen/COREAD/FEATURES/GEX_coread", row.names=1)
x <- GEX_coread [ , colnames(GEX_coread) %in% colnames(prot) ]

################################## convert all uniprot ID in gene name ##################################
library(clusterProfiler)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
uniprot_transformed <- read.csv("~/Documents/RWTH_Aachen/CONVERSION/uniprot_transformed", check.names = F)

x <- as.character(uniprot_transformed$Entry)
ids <- bitr(x, fromType="UNIPROT", toType="SYMBOL" , annoDb = "org.Hs.eg.db")
write.csv(ids, "/Users/miyang/Documents/RWTH_Aachen/CONVERSION/uniprot_ID_to_HUGO_ALL")

######################################### protein in complexes ##########################################
Prot_coreComplexes <- read.delim("~/Documents/RWTH_Aachen/GENERAL_DATA/Prot_coreComplexes.txt",check.names = F)
prot_name <- as.character(Prot_coreComplexes$`subunits(UniProt IDs)`)
prot_in_cplex <- unique( unlist(strsplit(prot_name, split=";")) ) 

conversion_table <- read.csv("/Users/miyang/Documents/RWTH_Aachen/CONVERSION/uniprot_ID_to_HUGO_ALL",row.names = 1)
prot_in_cplex_NAME <- as.character(conversion_table[conversion_table$UNIPROT %in% prot_in_cplex,"SYMBOL"])
write.csv(prot_in_cplex_NAME, "/Users/miyang/Documents/RWTH_Aachen/CONVERSION/prot_in_cplex_NAME")


