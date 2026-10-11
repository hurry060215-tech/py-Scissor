root <- dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[[1]])))
source(file.path(root,'compile.R'))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(jsonlite))
dir.create(file.path(root,'fixtures'),showWarnings=FALSE)
plain <- function(z) {
  if (inherits(z,'Matrix')) return(as.matrix(z))
  if (is.data.frame(z)) return(z)
  if (is.list(z)) return(lapply(z,plain))
  z
}
emit <- function(name,z) {
  jsonlite::write_json(plain(z),file.path(root,'fixtures',paste0(name,'.json')),pretty=TRUE,auto_unbox=TRUE,digits=16,na='string',null='null',matrix='rowmajor')
  cat('SAVED',name,'\n')
}
safe <- function(expr) tryCatch(list(ok=TRUE,result=plain(expr)),error=function(e) list(ok=FALSE,error=conditionMessage(e)))
set.seed(731)
n <- 24L; p <- 6L
x <- matrix(runif(n*p,-.7,.8),n,p)
colnames(x) <- paste0('cell_',seq_len(p)); rownames(x)<-paste0('sample_',seq_len(n))
graph <- matrix(0,p,p)
graph[1,2]<-graph[2,1]<-1
graph[1,3]<-graph[3,1]<-1
graph[2,3]<-graph[3,2]<-1
graph[2,4]<-graph[4,2]<-1
# Last two original nodes are isolated, avoiding the two observed index offsets.
dimnames(graph)<-list(colnames(x),colnames(x))
foldid <- rep(1:4,length.out=n)
ys <- list(gaussian=2+1.6*x[,1]-.9*x[,3]+.15*sin(seq_len(n)),
  binomial=as.numeric(1.6*x[,1]-.9*x[,3]+.4*sin(seq_len(n))>0),
  cox=cbind(time=rep(seq_len(n/2),each=2),status=as.numeric(seq_len(n)%%5!=0)))
params <- list(alpha=.4,lambda=c(.12,.06,.03,.015),nfolds=4L,foldid=foldid,inzero=FALSE,keep.beta=TRUE,maxit=100000L)
for (family in names(ys)) for (precision in c('default','tight')) {
  threshold <- if(precision=='default') 1e-7 else 1e-12
  inpt <- c(list(x=x,y=ys[[family]],family=family,Omega=graph),params,list(thresh=threshold))
  out <- safe(do.call(APML1,inpt))
  emit(paste(family,precision,'explicit_cv',sep='_'),list(provenance=list(upstream_sha='311560a1f160f665917e7da0fcf289914ee99773',R=R.version.string,modified_upstream=FALSE,scope='APML1 Net with explicit lambda/foldid; unused L0 trimming disabled'),inputs=inpt,output=out))
}
# The high-level default APML1 hard-threshold stage is separately observed, never silently replaced.
for(family in names(ys)) {
  inpt<-c(list(x=x,y=ys[[family]],family=family,Omega=graph),params,list(thresh=1e-7))
  inpt$inzero<-TRUE
  emit(paste0(family,'_default_with_trimming'),list(inputs=inpt,output=safe(do.call(APML1,inpt))))
}
# Actual R graph construction, and bounds diagnostics before ever calling an unsafe C++ access.
padded<-rbind(0,cbind(0,graph)); W<-OmegaC(padded,as.integer(rep(1,p+1))); W$loc<-W$loc+1
aligned<-OmegaC(graph,as.integer(rep(1,p))); aligned$loc<-aligned$loc+1
unsafe_graph<-graph; unsafe_graph[1,p]<-unsafe_graph[p,1]<-1
UW<-OmegaC(rbind(0,cbind(0,unsafe_graph)),as.integer(rep(1,p+1))); UW$loc<-UW$loc+1
emit('graph_index_audit',list(p=p,graph=graph,padded_normalized=as.matrix(W$Omega),padded_nadj=W$nadj,padded_loc=W$loc,aligned_normalized=as.matrix(aligned$Omega),unsafe_graph=unsafe_graph,unsafe_padded_nadj=UW$nadj,unsafe_padded_loc=UW$loc,unsafe_runtime_invoked=FALSE,note='C++ accesses beta[loc-1] and Gaussian inactive check beta[loc]; beta length is p. Out-of-bounds is diagnosed, never executed.'))
# Stable direct kernels distinguish graph alignment from the R wrapper's extra node.
for(family in c('gaussian','cox')) {
  y<-ys[[family]]; L<-params$lambda
  call_kernel<-function(G) {
    if(family=='gaussian') return(NetLmC(x,y,.4,L,length(L),0,rep(1,p),G$Omega,G$loc,G$nadj,p,n,1e-12,100000L,1e-5))
    prep<-PrepCox(x,y)
    NetCoxC(prep$x,prep$tevent,.4,L,length(L),0,rep(1,p),G$Omega,G$loc,G$nadj,prep$N,prep$nevent,prep$nevent1,prep$loc1,prep$n,p,n,1e-12,100000L,0L)
  }
  emit(paste0(family,'_graph_alignment_audit'),list(inputs=list(x=x,y=y,graph=graph,alpha=.4,lambda=L),official_padded=safe(call_kernel(W)),audit_aligned=safe(call_kernel(aligned)),note='Only direct-kernel graph input changed; upstream source was not patched. Aligned result is an audit, not a high-level official oracle.'))
}
# One predictor: closed-form Gaussian optimum is known, independent of any Python solver.
xx<-matrix(seq(-.6,.6,length.out=24),ncol=1)
yy<-1.7*xx[,1]+.2*cos(seq_len(24))
a<-.4; lam<-.04
xc<-xx[,1]-mean(xx[,1]); yc<-yy-mean(yy)
z<-mean(xc*yc); denom<-mean(xc^2)+lam*(1-a)
exact<-sign(z)*max(abs(z)-lam*a,0)/denom
gg<-OmegaC(matrix(0,1,1),as.integer(1));gg$loc<-gg$loc+1
full_default<-NetLmC(xx,yy,a,lam,1L,0L,1,gg$Omega,gg$loc,gg$nadj,1L,24L,1e-7,100000L,1e-5)
full_tight<-NetLmC(xx,yy,a,lam,1L,0L,1,gg$Omega,gg$loc,gg$nadj,1L,24L,1e-14,100000L,1e-14)
emit('gaussian_curvature_audit',list(inputs=list(x=xx,y=yy,alpha=a,lambda=lam),centered_second_moment=mean(xc^2),closed_form_coefficient=exact,official_default=full_default,official_tight=full_tight,note='Isolated-node I-Omega penalty; tests real stopping behavior with a known closed-form optimum.'))
# Cox internally rescales a by-value lambda, while wrapper uses original ratios in its CV.
prep<-PrepCox(x,ys$cox); xs<-sweep(prep$x,2,colMeans(prep$x))
lmax<-maxLambdaCoxC(xs,prep$tevent,prep$N,prep$nevent,prep$nevent1,prep$loc1,prep$n,.4,rep(1,p),n,p)
ratios<-c(1,.2,.04,.008)
invoke_cox<-function(L,automatic) NetCoxC(prep$x,prep$tevent,.4,L,length(L),as.integer(automatic),rep(1,p),W$Omega,W$loc,W$nadj,prep$N,prep$nevent,prep$nevent1,prep$loc1,prep$n,p,n,1e-10,100000L,0L)
emit('cox_lambda_audit',list(inputs=list(x=x,y=ys$cox,graph=graph,alpha=.4,ratios=ratios),lambda_max=lmax,actual_lambda=ratios*lmax,automatic=invoke_cox(ratios,TRUE),explicit_scaled=invoke_cox(ratios*lmax,FALSE),explicit_unscaled=invoke_cox(ratios,FALSE),wrapper_auto=safe(APML1(x,ys$cox,family='cox',Omega=graph,alpha=.4,nlambda=4,rlambda=.008,foldid=foldid,inzero=FALSE,keep.beta=TRUE,ifast=FALSE)),note='Compare automatic coefficient path to explicit scaled vs original-ratio paths. R wrapper lambda labels/CV grid use original ratio vector.'))
writeLines(capture.output(sessionInfo()),file.path(root,'sessionInfo.txt'))
