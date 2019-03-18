#-----------------------------------------------------------------------------------------------------------------------------------------
# This R script is therefore to generate all kinds of feature importance plots. You can choose for which drug in which kind of analysis,
# the plot will be drawn. 
# This script can be used to plot the following:
#     - feature importance plot for crossvalidation or bootstrapping analysis(depends on the function that is loaded):  
#         draws and saves a bar plot for the feature importance
#     - feature importance csv-file for bootstrapping analysis:  
#         creates a csv-file with the feature importance
#     - feature importance anova plot for bootstrapping analysis:  
#         creates a scatterplot and a density plot for the anova style of feature importance
#     - tendency plot for bootstrapping:
#         plots feature importance in order to compare the feature importance of two features related to Tissue
#-----------------------------------------------------------------------------------------------------------------------------------------
setwd("/nfs/nobackup2/saezgrp/homes/chris")

################################################################
# Load libraries
################################################################
library(Hmisc)
library(beeswarm)
library(hydroGOF)
# #################################################################################
# Choose the parameter for the feature importance plot, csv-file and the anova plot
# #################################################################################
stat_model<-c("RF")

analysis <- "nci_60_outerloop_v3"

drug_id <- c("12198", "701852", "374551", "732517")

featureSet <- c("SWATH")

balance <- c("balanced")

tissue_type <- c("all_tissue")

drug_mat <- c("fda")

# create feature importance plot for bootstrapping of crossvalidation

# bootstrapping
source("script/imp_plot/create_featImp_boot.R")
# crossvalidation
# source("script/imp_plot/create_featImp_cross.R")

for(dm in 1:length(drug_mat)){
  
  for(i in 1:length(featureSet)){
    
    for(dIdx in 1:length(drug_id)){
      for(f_Idx in 1:length(tissue_type)){
        for(b_Idx in 1:length(balance)){
          for(st_Idx in 1:length(stat_model)){
            # This path is the path to a directory in which a directory with the name of "featureSet" needs to exist
            search_path <- (paste(analysis, "/",  drug_mat[dm], "/", tissue_type[f_Idx], "/", balance[b_Idx], "/", featureSet[i], sep=""))
            # create FeatureImpPlot
            create_feature_imp(drug_id[dIdx], search_path,
                               paste("drug", drug_id[dIdx], "_", stat_model[st_Idx], "_", featureSet[i], ".eps",sep=""), 
                               stat_model[st_Idx], analysis, tissue_type[f_Idx], balance[b_Idx], drug_mat[dm])
          }
        }
      }
    }
  }
}


# create feature importance csv-file for boostrapping
# FUNCTION for feature importance csv-file
source("script/imp_plot/create_csv_file_featImp_boot.R")

for(dm in 1:length(drug_mat)){
  
  for(i in 1:length(featureSet)){
    
    for(dIdx in 1:length(drug_id)){
      for(f_Idx in 1:length(tissue_type)){
        for(b_Idx in 1:length(balance)){
          for(st_Idx in 1:length(stat_model)){
            # This path is the path to a directory in which a directory with the name of "featureSet" needs to exist
            search_path <- (paste(analysis, "/",  drug_mat[dm], "/", tissue_type[f_Idx], "/", balance[b_Idx], "/", featureSet[i], sep=""))
            # create FeatureImpPlot
            create_csv_file(drug_id[dIdx], search_path,
                               paste("drug", drug_id[dIdx], "_", stat_model[st_Idx], "_", featureSet[i], ".csv",sep=""), 
                               stat_model[st_Idx], analysis, tissue_type[f_Idx], balance[b_Idx], drug_mat[dm])
          }
        }
      }
    }
  }
}

# create an anova importance plot
# FUNCTION for anova plot
source("script(imp_plot/create_anova_plot_boot.R")

for(dm in 1:length(drug_mat)){
  
  for(i in 1:length(featureSet)){
    
    for(dIdx in 1:length(drug_id)){
      for(f_Idx in 1:length(tissue_type)){
        for(b_Idx in 1:length(balance)){
          for(st_Idx in 1:length(stat_model)){
            # This path is the path to a directory in which a directory with the name of "featureSet" needs to exist
            search_path <- (paste(analysis, "/",  drug_mat[dm], "/", tissue_type[f_Idx], "/", balance[b_Idx], "/", featureSet[i], sep=""))
            # create FeatureImpPlot
            create_anova_plot(drug_id[dIdx], search_path, 
                            stat_model[st_Idx], analysis, tissue_type[f_Idx], balance[b_Idx], drug_mat[dm], featureSet[i])
          }
        }
      }
    }
  }
}

# #################################################################################
# Choose the parameter for the tendency plot for bootstrapping
# #################################################################################

stat_model<-c("EN", "RF")

analysis <- "nci_60_outerloop_v3"

featureSet <- c("SWATH:highConf","GEX:NCI60","TISSUE:NCI60")

balance <- c("balanced")

tissue_type <- c("all_tissue")

drug_mat <- c("fda")

calc_tpye<-c("cindex","PEARSON","SPEARMAN")


# create tendency plot
# FUNCTION for tendency plot
source("script/imp_plot/create_tendency_plot_boot.R")

for(dm in 1:length(drug_mat)){
  for(cIdx in 1:length(calc_tpye)){
    for(f_Idx in 1:length(tissue_type)){
      for(b_Idx in 1:length(balance)){
        for(st_Idx in 1:length(stat_model)){
          tendency(stat_model[st_Idx], calc_tpye[cIdx], tissue_type[f_Idx], balance[b_Idx], analysis, featureSet, drug_mat[dm])
        }
      }
    }
  }
}