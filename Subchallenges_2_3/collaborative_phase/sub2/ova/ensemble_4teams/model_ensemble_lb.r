## name: model_ensemble_lb.r
## date: 02/15/2018

lb=c(0.4882,0.4611,0.4283,0.4293)
names(lb)=c("hyyf","hyu","dmis","deargen")

ensemble=as.matrix(read.delim(paste0("predictions_",names(lb)[1],".tsv"),sep="\t",check.names=F,row.names=1))
#ensemble=ensemble[rownames(ensemble)!="VGF",]
## for re-center ##
ref=ensemble
d1=dim(ref)[1];d2=dim(ref)[2]
mat_avg1=matrix(rep(apply(ref,1,mean,na.rm=T),times=d2),nrow=d1)
mat_sd1=matrix(rep(apply(ref,1,sd,na.rm=T),times=d2),nrow=d1)
###################
ensemble=ensemble*lb[1]
sum_lb=lb[1]
for(i in 2:length(lb)){
    pred=as.matrix(read.delim(paste0("predictions_",names(lb)[i],".tsv"),sep="\t",check.names=F,row.names=1))
    pred=pred[rownames(ensemble),colnames(ensemble)]
    ## re-center ##
    mat_avg2=matrix(rep(apply(pred,1,mean,na.rm=T),times=d2),nrow=d1)
    mat_sd2=matrix(rep(apply(pred,1,sd,na.rm=T),times=d2),nrow=d1)
    mat_sd2[mat_sd2==0]=1 # if sd==0, set sd to 1; e.g. dmis predicts ABHD5 all 0s
    pred=(pred-mat_avg2)/mat_sd2*mat_sd1+mat_avg1
    ###############
    ensemble=ensemble + pred*lb[i]

    sum_lb=sum_lb+lb[i]
    output=ensemble/sum_lb
    output=cbind(proteinID=rownames(ensemble),ensemble)
    write.table(output,file=paste0("predictions_ensemble_lb_top",i,".tsv"),quote=F,sep="\t",row.names=F)
}


