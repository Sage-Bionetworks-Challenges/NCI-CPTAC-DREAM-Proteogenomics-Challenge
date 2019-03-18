#--------------------------------------------------------------------------------------------------------------------------
#DESC:    This function performs a MultiDrug "nFold" crossvalidation. It?s doing a "random" crossvalidation
#         or a "leave_out_idep_entities" crossvalidation. It?s decided by calling this function, which type will
#         be performed. On default a "random" one is done. first->independent secound->dependent 
#IN:      idep_entities   => vector with independent entities/id?s, which will only appear in "test", "xTrain"
#                            or "train" set. None of the independent entities will be shared among these sets 
#         dep_entities    => the dependent entities/id`s will be duplicated in all three sets("test","xTrain","train")
#                            however the combination of independent and dependent entities is considered as a new
#                            unique entity
#         nFold           => numeric value in charge of which tpye of "nFold" crossvalidation will 
#                            be applied to "idep_entities"
#         names           => input in charge of the column-names of every set in each "fold_XXX"
#                            you need to commit two values or none, otherwise there will be an ERROR or the default case
#                            will be done(see ERROR-section below)
#                            on default the columns are named "idep_entities" and "dep_entities"
#         method          => string in charge of which type of MultiDrug CV will be performed
#                            on default a "random" MultiDrug CV is done, in order to do a "leave_out_idep_entities"
#                            MultiDrug CV the input for "method" needs to look like "leave_out_idep_entities"
#OUT:     Returns the list "cv", which contains as many lists as the size of "nFold". 
#         Each element in "cv" is called "fold_XXX".
#         Each of these "fold_XXX" contains three sets: "test", "xTrain" and "train".
#         Each of these sets being of the type data frame, and consisting out of two columns. The first and second column
#         comprises the independent and dependent entity identifiers, respectively. 
#ERROR:   This function throws an ERROR if:
#           - "dep_entities" isn?t unique or not of the type "numeric"
#           - "idep_entities" isn?t unique or not of the type "numeric"
#           - the input for "method" isn?t "random" or "leave_out_idep_entities" taking no account of use of 
#             capital and small letters. Only if there?s anything commited otherwise "method" == "random"(default)
#           - there are more or less than two arguments committed for "names". Only if there?s anything commited otherwise
#             "names" == c("idep_entities","dep_entities") (default)
#           - there is an ERROR occuring in the function "crossvaildation.R"
#--------------------------------------------------------------------------------------------------------------------------
multicross<-function(idep_entities, dep_entities, nFold, method="random", names=c("idep_entities", "dep_entities")){
  
  #TEST that dep_entities is unique
  if(length(unique(dep_entities))!=length(dep_entities)) stop("dep_entities is NOT unique",call. =FALSE)
  #TEST that dep_entities is numeric
  if(class(dep_entities)!="numeric") stop("dep_entities is NOT numeric",call. =FALSE)
  
  #TEST that idep_entities is unique
  if(length(unique(idep_entities))!=length(idep_entities)) stop("idep_entities is NOT unique",call. =FALSE)
  #TEST that idep_entities is numeric
  if(class(idep_entities)!="numeric") stop("idep_entities is NOT numeric",call. =FALSE)
  
  #TEST if the input for "method" is correct 
  method<-toupper(method)
  if(!((method=="RANDOM")||(method=="LEAVE_OUT_IDEP_ENTITIES"))) stop("input for 'method' isn?t correct",call. =FALSE)
  
  #TEST if both values for "names" are given otherwise "default case"
  if(length(names)!=2) stop("commit two or none values for 'name'",call. =FALSE)
  
  # "random" MultiDrug cross validation
  if (method == "RANDOM"){
    
    #create crossproduct of idep_entities and dep_entities
    cv<-expand.grid(idep_entities,dep_entities)
    
    #paste the two colomns of "cv" to one value seperated by "_"
    cv <- paste(cv[,1], cv[,2], sep="_")
    
    #perform "nFold" crossvalidation on "cv" and testing it
    cv<-crossvalidation(cv,nFold,TRUE)
    
    #split values of "cv" into two columns and create new dataframe "cv"
    for(i in 1:nFold){
      for (j in 1:length(cv[[i]])) {
        tmp <- strsplit(cv[[i]][[j]], "_")
        cv[[i]][[j]] <- data.frame(as.numeric(sapply(tmp, function(x) x[1])), 
                                   as.numeric(sapply(tmp, function(x) x[2])))
        
        #call the columns of "cv" after the input for "names"
        colnames(cv[[i]][[j]])<-names
      }
    }
  }
  
  
  # "leave_out_idep_entities" MultiDrug cross validation
  else if (method=="LEAVE_OUT_IDEP_ENTITIES"){
    
    #perform "nFOld" crossvalidation on "idep_entities" and testing it
    cv<-crossvalidation(idep_entities,nFold,TRUE)
    
    #doing crossproduct between the values from each set of each fold and dep_entities 
    for(i in 1:nFold){
      
      test_lo1<-expand.grid(cv[[i]]$test,dep_entities) 
      #call the columns after the input for "names" 
      colnames(test_lo1)<-names 
      
      xTrain_lo1<-expand.grid(cv[[i]]$xTrain,dep_entities)
      #call the columns after the input for "names" 
      colnames(xTrain_lo1)<-names
      
      train_lo1<-expand.grid(cv[[i]]$train,dep_entities)
      #call the columns after the input for "names" 
      colnames(train_lo1)<-names
      
      #creating the final crossvalidation list
      cv[[paste("fold_", i, sep="")]]<-list(test=test_lo1,xTrain=xTrain_lo1,train=train_lo1) 
    }
  }
  return(cv)
}