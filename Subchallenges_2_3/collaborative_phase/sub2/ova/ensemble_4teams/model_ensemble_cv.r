## name: model_ensemble_cv.r 
## date: 03/08/2018

p1=read.delim("cv_ova_hyyf.tsv",row.names=1,header=T)[,1:5]
p2=read.csv("cv_ova_hyu.csv",row.names=1)[,1:5]
p3=read.delim("cv_ova_sunkyu.tsv",row.names=1,header=T)[,1:5]
p4=read.csv("cv_ova_deargen.csv",row.names=1)[,1:5]

p2[is.na(p1)]=NA;p3[is.na(p1)]=NA;p4[is.na(p1)]=NA # exclude cases that have <5 observations

## 1. global weights (all genes same weight)
w1=mean(apply(p1,2,mean,na.rm=T)) #0.4646781
w2=mean(apply(p2,2,mean,na.rm=T)) #0.4928081
w3=mean(apply(p3,2,mean,na.rm=T)) #0.4801788
w4=mean(apply(p4,2,mean,na.rm=T)) #0.4453607

lb=c(w1,w2,w3,w4) # scores from training data cross-validation, not leaderboard
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
    write.table(output,file=paste0("predictions_ensemble_cv_top",i,".tsv"),quote=F,sep="\t",row.names=F)
}
ensemble=ensemble/sum(lb)

output=cbind(proteinID=rownames(ensemble),ensemble)
write.table(output,file="predictions_ensemble_global.tsv",quote=F,sep="\t",row.names=F)
###########################################

## 2. individual weights (one gene one weight)
w1=apply(p1,1,mean,na.rm=T);w1[is.na(w1)]=1
w2=apply(p2,1,mean,na.rm=T);w2[is.na(w2)]=1
w3=apply(p3,1,mean,na.rm=T);w3[is.na(w3)]=1
w4=apply(p4,1,mean,na.rm=T);w4[is.na(w4)]=1

lb=list(w1,w2,w3,w4) # scores from training data cross-validation, not leaderboard
names(lb)=c("hyyf","hyu","dmis","deargen")

ensemble=as.matrix(read.delim(paste0("predictions_",names(lb)[1],".tsv"),sep="\t",check.names=F,row.names=1))
#ensemble=ensemble[rownames(ensemble)!="VGF",]
## for re-center ##
ref=ensemble
d1=dim(ref)[1];d2=dim(ref)[2]
mat_avg1=matrix(rep(apply(ref,1,mean,na.rm=T),times=d2),nrow=d1)
mat_sd1=matrix(rep(apply(ref,1,sd,na.rm=T),times=d2),nrow=d1)
###################
w=matrix(rep(lb[[1]],time=dim(ensemble)[2]),nrow=length(lb[[1]]),ncol=dim(ensemble)[2]) # weight matrix
ensemble=ensemble*w
w_sum=w
for(i in 2:length(lb)){
    pred=as.matrix(read.delim(paste0("predictions_",names(lb)[i],".tsv"),sep="\t",check.names=F,row.names=1))
    pred=pred[rownames(ensemble),colnames(ensemble)]
    ## re-center ##
    mat_avg2=matrix(rep(apply(pred,1,mean,na.rm=T),times=d2),nrow=d1)
    mat_sd2=matrix(rep(apply(pred,1,sd,na.rm=T),times=d2),nrow=d1)
    mat_sd2[mat_sd2==0]=1 # if sd==0, set sd to 1; e.g. dmis predicts ABHD5 all 0s
    pred=(pred-mat_avg2)/mat_sd2*mat_sd1+mat_avg1
    ###############
    w=matrix(rep(lb[[i]],time=dim(ensemble)[2]),nrow=length(lb[[i]]),ncol=dim(ensemble)[2])
    ensemble=ensemble + pred*w
    w_sum=w_sum+w
}
ensemble=ensemble/w_sum

output=cbind(proteinID=rownames(ensemble),ensemble)
write.table(output,file="predictions_ensemble_individual.tsv",quote=F,sep="\t",row.names=F)

