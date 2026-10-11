root <- commandArgs(trailingOnly=TRUE)[1]
source(file.path(root,'compile.R'))
suppressPackageStartupMessages(library(Matrix)); suppressPackageStartupMessages(library(jsonlite))
plain <- function(z) {if(inherits(z,'Matrix'))return(as.matrix(z));if(is.data.frame(z))return(z);if(is.list(z))return(lapply(z,plain));z}
results<-list()
for(seed in 1:8) for(p in c(1,6)) for(family in c('gaussian','binomial','cox')) for(alpha in c(0,.4,1)) {
set.seed(seed);n<-24;x<-matrix(runif(n*p,-.8,.8),n,p);graph<-matrix(0,p,p)
if(p>3){graph[1,2]<-graph[2,1]<-1;graph[2,3]<-graph[3,2]<-1;graph[3,4]<-graph[4,3]<-1}
y<-switch(family,gaussian=2+x[,1]+.4*rnorm(n),binomial=as.numeric(x[,1]+rnorm(n)>.0),cox=cbind(time=sample(1:8,n,replace=TRUE),status=sample(c(0,1),n,replace=TRUE)))
inpt<-list(x=x,y=y,Omega=graph,family=family,alpha=alpha,lambda=c(.3,.1,.03,.01),thresh=1e-7,maxit=100000L,nfolds=1,inzero=FALSE,keep.beta=TRUE)
o<-tryCatch(plain(do.call(APML1,inpt)),error=function(e)list(error=conditionMessage(e)))
results[[length(results)+1]]<-list(seed=seed,inputs=inpt,output=o)
}
write_json(results,file.path(root,'review','random_paths.json'),auto_unbox=TRUE,digits=16,na='string',null='null',matrix='rowmajor')
