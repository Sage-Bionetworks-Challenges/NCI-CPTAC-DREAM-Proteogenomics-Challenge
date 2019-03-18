#-----------------------------------------------------------------------------------------------------------------------------------------
# This R script is therefore to generate all kinds of heatmaps. You can choose for which run of the pipeline the heatmaps should be 
# generated, by defining the string "analysis" with the name of the all sorounding folder for that run. By loading either the function
# for crossvalidation or bootstrapp, it's decided what kind of validation for that run has been done. Specify the "calc_type" string 
# in order to decide, which value of performance should be calculated between the observations and predictions. Those values will be 
# drawn in the heatmap. You also need to specify the statistical model and if the trainSet has been balanced or not for obtaining
# these results.
#-----------------------------------------------------------------------------------------------------------------------------------------
setwd("/nfs/nobackup2/saezgrp/homes/chris")

library(Hmisc)
library(pheatmap)
library(hydroGOF)

# ##############################################################################
# Define parameters for which settings a heatmap will be created
# ##############################################################################

# decide here which kind of value the heatmap should include
calc_type <- c("cindex", "PEARSON", "SPEARMAN")

# for each tissue type a heatmap will be created
tissue_type <- c("ALL", "BLCA", "BRCA", "COREAD", "DLBC", "ESCA", "GBM", "HNSC", "KIRC", "LAML", "LUAD", "MESO", "NB", "OV",
                 "PAAD", "SCLC", "SKCM", "STAD")

# a string witht the type of balancing 
balance <- c( "balanced")

# a string containing the statistical models
stat_model <- c( "EN","RF")

# the name of the folder where all drug response folder are included
analysis <- "boot_all_new"

# define the drug response matrice
drug_mat <- c("sanger", "fda", "moa")


#######################################
# create heatmap for crossvalidation
#######################################

# load necessary function
source("script/heatmap_fct/create_heatmap_cross.R")

# plot and save the heatmap
for(f_Idx in 1:length(tissue_type)){
  for(b_Idx in 1:length(balance)){
    for(st_Idx in 1:length(stat_model)){
      for(i in 1:length(calc_type)){
        for(dm in 1:length(drug_mat)){
          
          plot_heatmap(stat_model[st_Idx], calc_type[i], tissue_type[f_Idx], balance[b_Idx], analysis, drug_mat[dm])
          
        }
      }
    }
  }
}


#######################################
# create heatmap for bootstrapping
#######################################

# load necessary function
source("script/heatmap_fct/create_heatmap_boot.R")

# plot and save the heatmap
for(f_Idx in 1:length(tissue_type)){
  for(b_Idx in 1:length(balance)){
    for(st_Idx in 1:length(stat_model)){
      for(i in 1:length(calc_type)){
        for(dm in 1:length(drug_mat)){
          
          plot_heatmap_boot(stat_model[st_Idx], calc_type[i], tissue_type[f_Idx], balance[b_Idx], analysis, drug_mat[dm])
          
        }
      }
    }
  }
}
