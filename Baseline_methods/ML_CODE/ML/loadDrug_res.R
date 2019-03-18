loadDrug_res <-  function(drug_mat){
  
  # load challice drug_response data
  if(drug_mat=="challice"){
    load("data/pipeline/DREAM/challice.ro")
    drug_res <- challice
  }
  
  # load geneData_loewe drug_response data
  if(drug_mat=="geneData_loewe"){
    load("data/pipeline/DREAM/geneData_loewe.ro")
    drug_res <- geneData_loewe
  }
  
  # load geneData_bliss drug_response data
  if(drug_mat=="geneData_bliss"){
    load("data/pipeline/DREAM/geneData_bliss.ro")
    drug_res <- geneData_bliss
  }
  
  
  # load fda drug_response data
  if(drug_mat=="fda"){
    load("data/pipeline/drug/fda_drug_response.ro")
    drug_res <- fda
  }
  
  # load fda drug_response data
  if(drug_mat=="fda_v2"){
    load("data/pipeline/NCI60_final_v2/drugRes.ro")
    drug_res <- drugRes
  }
    
  # load moa response data
  if(drug_mat=="moa"){
    load("data/pipeline/drug/moa_drug_response.ro")
    drug_res <- moa
  }

  # load sanger response data v15
  if(drug_mat=="sanger_NCI60_IC50"){
    load("data/pipeline/NCI60_final/GDSC_IC50.ro")
    drug_res <- GDSC_IC50
  }
  
  # load sanger response data v15
  if(drug_mat=="sanger_NCI60_OMAUC"){
    load("data/pipeline/NCI60_final/GDSC_OMAUC.ro")
    drug_res <- GDSC_OMAUC
  }
  
  # load sanger response data v15
  if(drug_mat=="sanger"){
    load("data/pipeline/drug/DRUG_RESPONSE.ro")
    drug_res <- OMAUC_RELEASED
  }
  
  # load sanger response data v15
  if(drug_mat=="sanger_v15_IC50"){
    load("data/pipeline/drug/DRUG_RESPONSE.ro")
    drug_res <- IC50_RELEASED
  }
  
  # load sanger response data v16
  if(drug_mat=="sanger_v16"){
    load("data/pipeline/drug/DRUG_RESPONSE_v16.ro")
    drug_res <- OMAUC_CLASSIFIED
  }
  
  # load sanger response data v16
  if(drug_mat=="sanger_v16_IC50"){
    load("data/pipeline/drug/DRUG_RESPONSE_v16.ro")
    drug_res <- IC50_CLASSIFIED
  }
  
  # load sanger response data v16 (only public drugs)
  if(drug_mat=="sanger_v16_pub"){
    load("data/pipeline/drug/DRUG_RESPONSE_v16.ro")
    drug_res <- OMAUC_RELEASED
  }
  
  # load sanger response data v16 (only public drugs)
  if(drug_mat=="sanger_v17_OMAUC") {
    load("data/pipeline/v17/OMAUC_v17.ro")
    drug_res <- OMAUC_v17
  }
  
  return(drug_res)
}