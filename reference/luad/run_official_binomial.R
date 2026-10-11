args <- commandArgs(trailingOnly=TRUE)
script <- normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]))
here <- dirname(script); root <- dirname(dirname(here))
input <- if(length(args)>=1) normalizePath(args[1]) else file.path(root,'validation-data/luad/prepared')
output <- if(length(args)>=2) args[2] else file.path(root,'validation-data/luad/r-oracle/official-binomial')
Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',R_MAKEVARS_USER=file.path(here,'Makevars'))
.libPaths(c(file.path(root,'reference/library'),.libPaths()))
suppressPackageStartupMessages(library(Matrix)); suppressPackageStartupMessages(library(jsonlite))
source(file.path(here,'oracle.R'))
source(file.path(here,'compile_official.R')); compile_official(here,root); log_kernel_calls()
manifest <- read_json(file.path(input,'manifest.json'),simplifyVector=FALSE)
config <- manifest$families$binomial; cohort <- manifest$cohorts[[config$cohort]]
X <- read_f64(file.path(input,config$cohort,'X.bin'),unlist(cohort$X_shape)); y <- as.numeric(unlist(config$y)); foldid <- as.integer(unlist(config$foldid))
A <- as(readMM(file.path(input,'graph.mtx')),'dgCMatrix'); W <- normalize_adjacency(A)
stopifnot(ncol(X)==nrow(A),all(apply(X,2,sd)>0))
dir.create(output,recursive=TRUE,showWarnings=FALSE); writeLines(capture.output(sessionInfo()),file.path(output,'sessionInfo.txt'))
extract <- function(fit,lambdas,alpha) {
  b <- as.matrix(fit$Beta); B <- b[-1,,drop=FALSE]; intercept <- b[1,]
  list(coefficients=B,intercepts=intercept,upstream_flag=fit$flag,
       certificates=lapply(seq_len(ncol(B)),function(j) certificate(X,y,W,B[,j],intercept[j],lambdas[j],alpha,'binomial')))
}
for(path in config$paths) {
  alpha <- path$alpha; lambdas <- as.numeric(unlist(path$lambdas)); destination <- file.path(output,paste0('alpha-',format(alpha,scientific=FALSE)))
  dir.create(destination,recursive=TRUE,showWarnings=FALSE)
  cat('START unchanged APML1 binomial',alpha,'\n'); flush.console(); begin <- proc.time()
  original <- APML1(X,y,family='binomial',penalty='Net',Omega=A,alpha=alpha,lambda=lambdas,
      foldid=foldid,nfolds=length(unique(foldid)),inzero=FALSE,keep.beta=TRUE,thresh=1e-7,maxit=100000L)
  full <- extract(original,lambdas,alpha)
  # The unchanged R wrapper discards its fold objects. Re-run its exact direct
  # CV kernel calls to expose those coefficients and flags, and audit CV equality.
  padded <- rbind(0,cbind(0,A)); net <- OmegaSC(padded,rep(1L,ncol(X)+1L)); net$loc <- net$loc+1
  folds <- list(); losses <- matrix(NA_real_,length(unique(foldid)),length(lambdas))
  for(k in sort(unique(foldid))) {
    train <- foldid!=k; xtrain <- cbind(1,X[train,,drop=FALSE]); xtest <- cbind(1,X[!train,,drop=FALSE])
    stopifnot(all(apply(X[train,,drop=FALSE],2,sd)>0))
    z <- cvNetLogC(xtrain,y[train],alpha,lambdas,length(lambdas),c(0,rep(1,ncol(X))),c(0,rep(1,ncol(X))),
        net$Omega,net$loc,net$nadj,ncol(X)+1L,sum(train),1e-7,100000L,xtest,y[!train],sum(!train),1e-5)
    q <- z$nlambda; folds[[k]] <- list(coefficients=z$Beta[-1,seq_len(q),drop=FALSE],intercepts=z$Beta[1,seq_len(q)],upstream_flag=z$flag[seq_len(q)],effective_path_length=q)
    folds[[k]]$certificates <- lapply(seq_len(q),function(j) certificate(X[train,,drop=FALSE],y[train],W,folds[[k]]$coefficients[,j],folds[[k]]$intercepts[j],lambdas[j],alpha,'binomial'))
    losses[k,seq_len(q)] <- -2*z$LLF[seq_len(q)]
    cat('unchanged binomial',alpha,'fold',k,'path length',q,'\n'); flush.console()
  }
  weights <- as.numeric(table(foldid)); cvm <- colSums(losses*weights)/sum(weights)
  cvse <- sqrt(colSums(sweep(losses,2,cvm,'-')^2*weights)/sum(weights)/(length(weights)-1))
  stopifnot(max(abs(cvm-original$fit$cvm))<1e-12,max(abs(cvse-original$fit$cvse))<1e-12)
  refit_original <- APML1(X,y,family='binomial',penalty='Net',Omega=A,alpha=alpha,lambda=original$lambda.min,
      nfolds=1,inzero=FALSE,keep.beta=TRUE,thresh=1e-7,maxit=100000L)
  refit <- extract(refit_original,original$lambda.min,alpha)
  out <- list(full=full,folds=folds,refit=refit,fit=original$fit,fold_losses=losses,fold_weights=weights,
      cvm=cvm,cvse=cvse,lambda_min=original$lambda.min,
      metadata=list(family='binomial',alpha=alpha,lambdas=lambdas,foldid=foldid,thresh=1e-7,maxit=100000,
        n_samples=nrow(X),n_cells=ncol(X),R=R.version.string,source_commit='311560a1f160f665917e7da0fcf289914ee99773',
        semantics='unchanged APML1 binomial',upstream_source_patched=FALSE,full_dataset=TRUE,no_sampling=TRUE,
        compiler_flags='-O2 -g0; no fast-math',logging='external generated-Rcpp-entrypoint wrappers; arguments/results unchanged',
        cv_fold_export='exact published direct-kernel CV calls rerun; aggregate equality checked against unchanged wrapper',
        elapsed_seconds=unname((proc.time()-begin)['elapsed'])))
  export_cv(out,destination)
  cat('DONE unchanged APML1 binomial',alpha,'seconds',out$metadata$elapsed_seconds,'\n'); flush.console()
}
