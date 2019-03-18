load("~/Desktop/project/data/pipeline/drug/DRUG_RESPONSE.ro")
source("~/Documents/DRUG_SENSI_GIT/R/balancing/LIB_BALANCING.R")

# ----------------------------------------------------------------------------------------
# Demo for balancing drug response

drugIDX <- 5
res <- OMAUC_RELEASED[drugIDX,]
res <- res[!is.na(res)]

hist(res)
plot(density(res))

hist(balance(res))
plot(density(balance(res)))

