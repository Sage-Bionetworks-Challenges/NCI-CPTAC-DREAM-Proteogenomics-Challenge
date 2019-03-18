#------------------------------------------------------------------------------------------------------------
#DESC:    This R-script is a demo script for the crossvalidation mechanism.
#         It runs all the three kinds of crossvalidation: SingleDrug, MultiDrug "random" and 
#         MultiDrug "leave_out_idep_entities". The demo returns three lists, each one specific for
#         one type of crossvalidation.
#IN:      The script runs the three type of crossvalidation for given exmaple values:
#             - "entity_single" and "nFold_single" for the SingleDrug CV
#             - "idep_entities", "dep_entities" and "nFold_multi" for MultiDrug CV
#OUT:     The script creates three lists:
#             - list "single" for SingleDrug CV
#             - list "multi_random" for MultiDrug "random" CV
#             - list "multi_lo1" for MultiDrug "leave_out_idep_entities" CV 
#------------------------------------------------------------------------------------------------------------

#load neccessary functions
source("~/../Dropbox/withMicha/R/crossvalidate_data/02_single_drug_CV.R")
source("~/../Dropbox/withMicha/R/crossvalidate_data/03_multi_drug_CV.R")

#TEST-section for SingleDrug

#example values for SingleDrug
entity_single<-1000:1100
entity_single<-as.numeric(entity_single)
nFold_single<-10

#call function "singlecross" for a "nFold_single" crossvalidation on "entity_single"
single<-singlecross(entity_single,nFold_single)


#TEST-section for MultiDrug 

#example values for MultiDrug
idep_entities<- 1:1100
idep_entities<-as.numeric(idep_entities)
dep_entities<- 1:400
dep_entities<-as.numeric(dep_entities)
nFold_multi<-10

#call function "multicross" for a "nFold_multi" random crossvalidation 
multi_random<-multicross(idep_entities,dep_entities,nFold_multi,"random")

#call function "multicross" for a "nFold_multi" leave_out_idep_entities crossvalidation
multi_lo1<-multicross(idep_entities,dep_entities,nFold_multi,"leave_out_idep_entities")

#remove unimportant files from workspace
rm(dep_entities,idep_entities,entity_single)