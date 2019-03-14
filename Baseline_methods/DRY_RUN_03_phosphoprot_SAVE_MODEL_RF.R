test.I<-function(sem){ # --- run in parallel
  
Y.phospho<-as.matrix(read.table(paste(".../Data/retrospective_phospho.txt",sep=""),stringsAsFactors =FALSE,sep="\t",header=TRUE))
Y.pro<-as.matrix(read.table(paste(".../Data/retrospective_proteome.txt",sep=""),stringsAsFactors =FALSE,sep="\t",header=TRUE))
phospho.name<-Y.phospho[,1]; global.name<-Y.pro[,1]

# -- match samples
phospho.ID<-colnames(Y.phospho); global.ID<-colnames(Y.pro)
m.g<-match(phospho.ID,global.ID)
Y.pro<-Y.pro[,m.g[is.na(m.g)==FALSE]]
Y.phospho<-Y.phospho[,is.na(m.g)==FALSE]
Y.phospho<-Y.phospho[,-1]; Y.pro<-Y.pro[,-1]

index<-(rowSums(is.na(Y.phospho)))
phospho.name<-phospho.name[index<dim(Y.phospho)[2]-3]
Y.phospho<-Y.phospho[index<dim(Y.phospho)[2]-3,]

class(Y.phospho)<-"numeric"
class(Y.pro)<-"numeric"

# -- proteomics filtering
protein.subset<-read.table(paste(".../Data/no_na_prospective_proteome_id.txt",sep=""),stringsAsFactors =FALSE,sep="\t",header=TRUE)[,1]
m.g<-match(protein.subset,global.name); Y.pro<-Y.pro[m.g,]

index<-(rowSums(is.na(Y.pro))==0)
Y.pro<-Y.pro[index,]; protein.subset<-protein.subset[index]

library("randomForest")

iter.max<-100
init<-c((sem-1)*iter.max+1,sem*iter.max)
if (init[2]>dim(Y.phospho)[1]) init[2]<-dim(Y.phospho)[1]

out<-vector("list", init[2]-init[1])
ss=0
for (k in init[1]:init[2]){
  ss=ss+1
  index<-(is.na(Y.phospho[k,])==FALSE)
  out[[ss]]<-randomForest(x=t(Y.pro[,index]),y=Y.phospho[k,index],nTree=1000)
  print(k)
} 

phospho.reg<-phospho.name[seq(init[1],init[2])]
save(out,phospho.reg, file=paste(".../Models/Model_",sem,".rda",sep=""))

}
