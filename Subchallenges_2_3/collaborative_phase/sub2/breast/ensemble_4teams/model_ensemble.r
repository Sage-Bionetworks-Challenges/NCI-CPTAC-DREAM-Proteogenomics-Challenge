## name: model_ensemble.r
## date: 02/15/2018

lb=c(0.4882,0.4611,0.4283,0.4293)
names(lb)=c("hyyf","hyu","dmis","deargen")

ensemble=as.matrix(read.delim(paste0("predictions_",names(lb)[1],".tsv"),sep="\t",check.names=F,row.names=1))
ensemble=ensemble[rownames(ensemble)!="VGF",]
ensemble=ensemble*lb[1]
for(i in 2:length(lb)){
	pred=as.matrix(read.delim(paste0("predictions_",names(lb)[i],".tsv"),sep="\t",check.names=F,row.names=1))
	pred=pred[rownames(ensemble),colnames(ensemble)]
	ensemble=ensemble + pred*lb[i]
}
ensemble=ensemble/sum(lb)


output=cbind(proteinID=rownames(ensemble),ensemble)
write.table(output,file="prediction_ensemble.tsv",quote=F,sep="\t",row.names=F)

