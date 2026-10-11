script <- normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]))
here <- dirname(script); root <- dirname(dirname(here))
Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',R_MAKEVARS_USER=file.path(here,'Makevars'))
.libPaths(c(file.path(root,'reference/library'),.libPaths()))
suppressPackageStartupMessages(library(Matrix)); suppressPackageStartupMessages(library(jsonlite))
source(file.path(here,'oracle.R'))
dir.create(file.path(here,'cache'),recursive=TRUE,showWarnings=FALSE)
Rcpp::sourceCpp(file.path(here,'corrected_oracle.cpp'),cacheDir=file.path(here,'cache'),rebuild=FALSE)
out <- list(R=R.version.string,scope='small mathematical self-tests of independent oracle; these are not full LUAD validation')
x <- matrix(c(-2,-1,0,1,2)*0.7,ncol=1); y <- 2+3*x[,1]; W <- sparseMatrix(i=integer(),j=integer(),dims=c(1,1))
for(alpha in c(0,0.4,1)) {
  lambda <- c(1,0.3,0.01); fit <- fit_corrected(x,y,W,lambda,alpha,'gaussian',1e-10)
  expected <- pmax(mean(x[,1]*y)-lambda*alpha,0)/(mean(x^2)+lambda*(1-alpha))
  err <- max(abs(fit$coefficients-expected),abs(fit$intercepts-2))
  stopifnot(all(fit$converged),err<1e-8)
  out[[paste0('closed_form_gaussian_alpha_',alpha)]] <- list(max_error=err,kkt=fit$kkt_residual)
}
# Independent base-R smooth minimization and exact R risk-set certificates.
set.seed(491); n <- 48; p <- 6; X <- matrix(rnorm(n*p),n,p)
A <- sparseMatrix(i=c(1:5,2:6),j=c(2:6,1:5),x=1,dims=c(p,p)); W <- normalize_adjacency(A)
responses <- list(gaussian=1+X[,1]-X[,3]+rnorm(n)*0.2,
                  binomial=as.numeric(runif(n)<plogis(0.3+X[,1]-0.7*X[,2])),
                  cox=cbind(time=sample(1:15,n,replace=TRUE),status=rep(c(0,1,1),length.out=n)))
for(family in names(responses)) for(alpha in c(0,0.4,1)) {
  y <- responses[[family]]; lambda <- 0.15
  fit <- fit_corrected(X,y,W,lambda,alpha,family,1e-9)
  stopifnot(all(fit$converged),fit$certificate_kkt_disagreement<1e-10,fit$certificate_objective_disagreement<1e-10)
  # Split b into nonnegative b+ and b-, optimizing with the base-R bound solver.
  decode <- function(z) list(b=z[seq_len(p)]-z[p+seq_len(p)],a=if(family=='cox') 0 else z[2*p+1])
  objective <- function(z) { d <- decode(z); certificate(X,y,W,d$b,d$a,lambda,alpha,family)$objective-lambda*alpha*sum(abs(d$b))+lambda*alpha*sum(z[seq_len(2*p)]) }
  start <- c(rep(0,2*p),if(family!='cox') 0)
  reference <- optim(start,objective,method='L-BFGS-B',lower=c(rep(0,2*p),if(family!='cox') -Inf),
                     control=list(factr=1e3,pgtol=1e-9,maxit=10000,ndeps=rep(1e-6,length(start))))
  delta <- fit$certificates[[1]]$objective-reference$value
  cat('SELFTEST',family,alpha,'delta',delta,'ref convergence',reference$convergence,'KKT',fit$kkt_residual,'\n'); flush.console()
  stopifnot(abs(delta)<1e-7)
  out[[paste(family,alpha,sep='_')]] <- list(objective_difference_vs_base_R=delta,base_R_convergence=reference$convergence,
       kkt=fit$kkt_residual,certificate_kkt_disagreement=fit$certificate_kkt_disagreement,
       certificate_objective_disagreement=fit$certificate_objective_disagreement)
}
duplicated_example <- matrix(1:8,4,2,dimnames=list(c('A','B','A','C'),c('s1','s2')))
common <- intersect(rownames(duplicated_example),c('A','B','D'))
stopifnot(identical(common,c('A','B')),identical(unname(duplicated_example[common,,drop=FALSE]),unname(duplicated_example[1:2,,drop=FALSE])))
out$R_character_indexing <- 'intersect is unique in first-input order; duplicate matrix row names select first matching source row'
destination <- file.path(root,'validation-data/luad/r-oracle'); dir.create(destination,recursive=TRUE,showWarnings=FALSE)
json(out,file.path(destination,'selftest.json')); cat('ALL INDEPENDENT ORACLE SELFTESTS PASSED\n')
