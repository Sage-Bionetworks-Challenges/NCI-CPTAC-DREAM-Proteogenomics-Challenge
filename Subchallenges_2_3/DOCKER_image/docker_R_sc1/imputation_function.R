require('pamr');

my.normlize = function(data)
{
  median.temp = apply(data,1,median,na.rm=T);
  sd.temp = apply(data,1,sd,na.rm=T);
  # apply(t(t(data)-median.temp),2,median,na.rm=T);
  # apply(t((t(data)-median.temp)/sd.temp),2,sd,na.rm=T);
  data.new = ((data-median.temp)/sd.temp);
  return(list(data = data.new, median = median.temp, sd=sd.temp));
}

my.normlize.rev = function(data, median.old, sd.old)
{
  data.new = (data*sd.old+median.old);
}

my.imputation = function(data,k=10)
{
  norm.temp = my.normlize(data)
  
  data.new = norm.temp[[1]]
  median.new = norm.temp[[2]]
  sd.new = norm.temp[[3]]
  
  #set.seed(1387)
  impu.temp = pamr.knnimpute(list(x =  as.matrix(data.new),y = rep(1,dim(data)[2])),k,rowmax = 0.95,colmax = 0.95)
  
  X1b.new = my.normlize.rev(impu.temp[[1]], median.new, sd.new)
  X1b.new[X1b.new<0] = 0
  return(X1b.new)
}


# library('missForest');
# 
# my.normlize = function(data)
# {
#   median.temp = apply(data,1,median,na.rm=T);
#   sd.temp = apply(data,1,sd,na.rm=T);
#   # apply(t(t(data)-median.temp),2,median,na.rm=T);
#   # apply(t((t(data)-median.temp)/sd.temp),2,sd,na.rm=T);
#   data.new = ((data-median.temp)/sd.temp);
#   return(list(data = data.new, median = median.temp, sd=sd.temp));
# }
# 
# my.normlize.rev = function(data, median.old, sd.old)
# {
#   data.new = (data*sd.old+median.old);
# }
# 
# my.imputation = function(data,k=10)
# {
#   print(paste("Max Num of Missing Values per protein",max(colSums(is.na(data)))))
#   impu.temp<-missForest(data, maxiter = 10, ntree = 100, variablewise = FALSE)
#   
#   return(impu.temp$ximp);  
# }
