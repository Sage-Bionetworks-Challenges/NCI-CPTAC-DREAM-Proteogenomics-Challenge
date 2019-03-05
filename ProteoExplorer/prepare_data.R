
## swissprot_ID_HGNC got it from the HeLa App

## Training data (retrospective)
path <- "~/Documents/RWTH_Aachen/DREAM_CPTAC/PROT_PREDICTION/DATA/synapse_retrospective/"

retrospective_breast_CNA_sort_common_gene_16884  <- read.csv(paste0(path,"retrospective_breast_CNA_sort_common_gene_16884.txt"), row.names= 1, sep="\t", check.names = F )
retrospective_breast_phospho_sort_common_gene_31981  <- read.csv(paste0(path,"retrospective_breast_phospho_sort_common_gene_31981.txt"), row.names= 1, sep="\t", check.names = F )
retrospective_breast_proteome_sort_common_gene_10006  <- read.csv(paste0(path,"retrospective_breast_proteome_sort_common_gene_10006.txt"), row.names= 1, sep="\t", check.names = F )
retrospective_breast_RNA_sort_common_gene_15107  <- read.csv(paste0(path,"retrospective_breast_RNA_sort_common_gene_15107.txt"), row.names= 1, sep="\t", check.names = F )

retrospective_ova_array_sort_common_gene_15121  <- read.csv(paste0(path,"retrospective_ova_array_sort_common_gene_15121.txt"), row.names= 1, sep="\t", check.names = F )
retrospective_ova_CNA_sort_common_gene_11859  <- read.csv(paste0(path,"retrospective_ova_CNA_sort_common_gene_11859.txt"), row.names= 1, sep="\t", check.names = F )
retrospective_ova_JHU_proteome_sort_common_gene_7061  <- read.csv(paste0(path,"retrospective_ova_JHU_proteome_sort_common_gene_7061.txt"), row.names= 1, sep="\t", check.names = F )
retrospective_ova_phospho_sort_common_gene_10057 <- read.csv(paste0(path,"retrospective_ova_phospho_sort_common_gene_10057.txt"), row.names= 1, sep="\t", check.names = F )
retrospective_ova_PNNL_proteome_sort_common_gene_7061  <- read.csv(paste0(path,"retrospective_ova_PNNL_proteome_sort_common_gene_7061.txt"), row.names= 1, sep="\t", check.names = F )
retrospective_ova_rna_seq_sort_common_gene_15121  <- read.csv(paste0(path,"retrospective_ova_rna_seq_sort_common_gene_15121.txt"), row.names= 1, sep="\t", check.names = F )

Training <- list(retrospective_breast_CNA_sort_common_gene_16884,retrospective_breast_phospho_sort_common_gene_31981,retrospective_breast_proteome_sort_common_gene_10006,
                 retrospective_breast_RNA_sort_common_gene_15107,retrospective_ova_array_sort_common_gene_15121,retrospective_ova_CNA_sort_common_gene_11859,
                 retrospective_ova_JHU_proteome_sort_common_gene_7061,retrospective_ova_phospho_sort_common_gene_10057,retrospective_ova_PNNL_proteome_sort_common_gene_7061,
                 retrospective_ova_rna_seq_sort_common_gene_15121)

names(Training) <- c("retrospective_breast_CNA_sort_common_gene_16884","retrospective_breast_phospho_sort_common_gene_31981","retrospective_breast_proteome_sort_common_gene_10006","retrospective_breast_RNA_sort_common_gene_15107","retrospective_ova_array_sort_common_gene_15121","retrospective_ova_CNA_sort_common_gene_11859","retrospective_ova_JHU_proteome_sort_common_gene_7061","retrospective_ova_phospho_sort_common_gene_10057","retrospective_ova_PNNL_proteome_sort_common_gene_7061","retrospective_ova_rna_seq_sort_common_gene_15121")
save(Training, file="/Users/miyang/Documents/RWTH_Aachen/DREAM_CPTAC/ProteoEstimator/ProteoEstimator/data/Training.Rdata")


## Testing data (prospective)
path <- "/Users/miyang/Documents/RWTH_Aachen/DREAM_CPTAC/DATA_prospective/SC2/"

prospective_breast_CNA_sort_common_gene_16884  <- read.csv(paste0(path,"prospective_breast_CNA_sort_common_gene_16884.txt"), row.names= 1, sep="\t", check.names = F )
prospective_breast_proteome_sort_common_gene_10005  <- read.csv(paste0(path,"prospective_breast_proteome_sort_common_gene_10005.txt"), row.names= 1, sep="\t", check.names = F )
prospective_breast_RNA_sort_common_gene_15107  <- read.csv(paste0(path,"prospective_breast_RNA_sort_common_gene_15107.txt"), row.names= 1, sep="\t", check.names = F )

prospective_ova_CNA_median_sort_common_gene_11859  <- read.csv(paste0(path,"prospective_ova_CNA_median_sort_common_gene_11859.txt"), row.names= 1, sep="\t", check.names = F )
prospective_ova_proteome_sort_common_gene_7061  <- read.csv(paste0(path,"prospective_ova_proteome_sort_common_gene_7061.txt"), row.names= 1, sep="\t", check.names = F )
prospective_ova_rna_seq_sort_common_gene_15121  <- read.csv(paste0(path,"prospective_ova_rna_seq_sort_common_gene_15121.txt"), row.names= 1, sep="\t", check.names = F )

Testing <- list(prospective_breast_CNA_sort_common_gene_16884,prospective_breast_proteome_sort_common_gene_10005,prospective_breast_RNA_sort_common_gene_15107,
                prospective_ova_CNA_median_sort_common_gene_11859,prospective_ova_proteome_sort_common_gene_7061,prospective_ova_rna_seq_sort_common_gene_15121 )
for(i in 1:length(Testing)) { # i=1
  colnames(Testing[[i]]) <- gsub("Participant","Patient",colnames(Testing[[i]]), fixed=TRUE) 
 }
names(Testing) <- c("prospective_breast_CNA_sort_common_gene_16884","prospective_breast_proteome_sort_common_gene_10005","prospective_breast_RNA_sort_common_gene_15107","prospective_ova_CNA_median_sort_common_gene_11859","prospective_ova_proteome_sort_common_gene_7061","prospective_ova_rna_seq_sort_common_gene_15121")
save(Testing, file="/Users/miyang/Documents/RWTH_Aachen/DREAM_CPTAC/ProteoEstimator/ProteoEstimator/data/Testing.Rdata")

scored_breast_protein  <- read.csv(paste0(path,"breast_protein.txt"),  sep="\t", check.names = F ) ; scored_breast_protein <- as.character(scored_breast_protein$IDs)
scored_ovarian_protein  <- read.csv(paste0(path,"ovarian_protein.txt"),  sep="\t", check.names = F ) ; scored_ovarian_protein <- as.character(scored_ovarian_protein$IDs)
scored <- list(scored_breast_protein,scored_ovarian_protein)
names(scored) <- c("scored_breast_protein","scored_ovarian_protein")
save(scored, file="/Users/miyang/Documents/RWTH_Aachen/DREAM_CPTAC/ProteoEstimator/ProteoEstimator/data/scored.Rdata")


## Prediction output of top performing teams
path <- "/Users/miyang/Documents/RWTH_Aachen/DREAM_CPTAC/CHALLENGE_RESULT/finalround/PREDICTIONS/"

load(paste0(path, "sc2_breast/mat_list.Rdata"))
team_Guan_breast <- mat[["Hongyang Li and Yuanfang Guan"]]
team_DMIS_PTG_breast <- mat[["DMIS_PTG"]]
team_hyu_breast <- mat[["hyu"]]
team_DEARGENpg_breast <- mat[["DEARGENpg"]]

load(paste0(path, "sc2_ovarian/mat_list.Rdata"))
team_Guan_ovarian <- mat[["Hongyang Li and Yuanfang Guan"]]
team_DMIS_PTG_ovarian <- mat[["DMIS_PTG"]]
team_hyu_ovarian <- mat[["hyu"]]
team_DEARGENpg_ovarian <- mat[["DEARGENpg"]]

ensemble_top4_breast   <- read.delim(paste0(path,"sc2_breast/predictions_ensemble_cv_top4.tsv"), row.names = 1 )
ensemble_top4_ovarian  <- read.delim(paste0(path,"sc2_ovarian/predictions_ensemble_cv_top4.tsv"), row.names = 1 )

Prediction <- list(team_Guan_breast,team_DMIS_PTG_breast,team_hyu_breast,team_DEARGENpg_breast,team_Guan_ovarian,team_DMIS_PTG_ovarian,team_hyu_ovarian,team_DEARGENpg_ovarian,
                   ensemble_top4_breast,ensemble_top4_ovarian)

for(i in 1:length(Prediction)) { # i=1
  colnames(Prediction[[i]]) <- gsub("Participant","Patient",colnames(Prediction[[i]]), fixed=TRUE) 
}

names(Prediction) <- c("team_Guan_breast","team_DMIS_PTG_breast","team_hyu_breast","team_DEARGENpg_breast","team_Guan_ovarian","team_DMIS_PTG_ovarian","team_hyu_ovarian","team_DEARGENpg_ovarian","ensemble_top4_breast","ensemble_top4_ovarian")
save(Prediction, file="/Users/miyang/Documents/RWTH_Aachen/DREAM_CPTAC/ProteoEstimator/ProteoEstimator/data/Prediction.Rdata")

path <- "/Users/miyang/Documents/RWTH_Aachen/DREAM_CPTAC/CHALLENGE_RESULT/finalround/"
guan_sc2_breast_cor <- read.csv(paste0(path,"reguanteamspredictionfiles/guan_sc2_breast_cor.txt"), sep="")
guan_sc2_ovarian_cor <- read.csv(paste0(path,"reguanteamspredictionfiles/guan_sc2_ovarian_cor.txt"), sep="")

Prediction_performance_per_gene <- list(guan_sc2_breast_cor,guan_sc2_ovarian_cor)
names(Prediction_performance_per_gene) <- c("team_Guan_breast","team_Guan_ovarian")
save(Prediction_performance_per_gene, file="/Users/miyang/Documents/RWTH_Aachen/DREAM_CPTAC/ProteoEstimator/ProteoEstimator/data/Prediction_performance_per_gene.Rdata")


