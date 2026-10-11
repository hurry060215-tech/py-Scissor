root<-commandArgs(trailingOnly=TRUE)[1];source(file.path(root,'compile.R'))
suppressPackageStartupMessages(library(Matrix));suppressPackageStartupMessages(library(jsonlite))
plain<-function(z){if(inherits(z,'Matrix'))return(as.matrix(z));if(is.data.frame(z))return(z);if(is.list(z))return(lapply(z,plain));z}
results<-list(); n<-25;p<-6
for(seed in 1:10)for(family in c('gaussian','binomial','cox'))for(alpha in c(.4,1)){
set.seed(seed);x<-matrix(runif(n*p,-.8,.8),n,p);graph<-matrix(0,p,p);graph[1,2]<-graph[2,1]<-1;graph[2,3]<-graph[3,2]<-1;graph[3,4]<-graph[4,3]<-1
foldid<-rep(c(1,2,3),length.out=n);foldid[1:3]<-1
 y<-switch(family,gaussian=2+x[,1]+.4*rnorm(n),binomial=as.numeric(x[,1]+rnorm(n)>.0),cox=cbind(time=sample(1:8,n,replace=TRUE),status=as.numeric(seq_len(n)%%3!=0)))
# avoid zero held-out events at fold 3
if(family=='cox')y[c(6,12,18,24),2]<-1
inpt<-list(x=x,y=y,Omega=graph,family=family,alpha=alpha,lambda=c(.3,.1,.03,.01),thresh=1e-7,maxit=100000L,foldid=foldid,inzero=FALSE,keep.beta=TRUE)
o<-tryCatch(plain(do.call(APML1,inpt)),error=function(e)list(error=conditionMessage(e)))
results[[length(results)+1]]<-list(seed=seed,inputs=inpt,output=o)
}
write_json(results,file.path(root,'review','unequal_fold_cv.json'),auto_unbox=TRUE,digits=16,na='string',null='null',matrix='rowmajor')
