root <- dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[[1]])))
source(file.path(root,'compile.R'))
library(jsonlite)
library(Matrix)
x<-matrix(seq(-.6,.6,length.out=24),ncol=1); y<-1.7*x[,1]+.2*cos(seq_len(24))
train<-seq_len(24)%%4!=0
g<-OmegaC(matrix(0,1,1),as.integer(1));g$loc<-g$loc+1
run<-function(thresh,thresh2) cvNetLmC(x[train,,drop=FALSE],y[train],.4,.04,1L,1,g$Omega,g$loc,g$nadj,sum(train),1L,thresh,100000L,x[!train,,drop=FALSE],y[!train],sum(!train),thresh2)
xc<-x[train,1]-mean(x[train,1]);yc<-y[train]-mean(y[train]);z<-mean(xc*yc)
exact<-sign(z)*max(abs(z)-.04*.4,0)/(mean(xc^2)+.04*.6)
out<-list(inputs=list(x=x,y=y,train=train,alpha=.4,lambda=.04),closed_form_coefficient=exact,official_default=run(1e-7,1e-5),official_APML_tighter_thresh=run(1e-12,1e-5),audit_all_tolerances_tight=run(1e-14,1e-14))
jsonlite::write_json(out,file.path(root,'fixtures','gaussian_cv_curvature_audit.json'),pretty=TRUE,auto_unbox=TRUE,digits=16,matrix='rowmajor')
