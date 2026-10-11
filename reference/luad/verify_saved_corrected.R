# Re-open every saved coefficient, recompute every certificate and CV statistic.
# This is an output verification, never an implicit refit or sampled check.
args <- commandArgs(trailingOnly=TRUE); family <- if(length(args)) args[1] else 'gaussian'
script <- normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]))
here <- dirname(script); root <- dirname(dirname(here))
Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
.libPaths(c(file.path(root,'reference/library'),.libPaths()))
suppressPackageStartupMessages(library(Matrix)); suppressPackageStartupMessages(library(jsonlite))
source(file.path(here,'oracle.R'))
input <- file.path(root,'validation-data/luad/prepared'); output <- file.path(root,'validation-data/luad/r-oracle')
manifest <- read_json(file.path(input,'manifest.json'),simplifyVector=FALSE); config <- manifest$families[[family]]
X <- read_f64(file.path(input,config$cohort,'X.bin'),unlist(manifest$cohorts[[config$cohort]]$X_shape))
y <- if(family=='cox') do.call(rbind,config$y) else unlist(config$y); storage.mode(y) <- 'double'
foldid <- unlist(config$foldid); W <- normalize_adjacency(readMM(file.path(input,'graph.mtx')))
verification <- list(status='passed',family=family,R=R.version.string,scope='every saved full/fold/refit coefficient and every lambda; independent R certificates and CV recomputed',alphas=list())
for(path in config$paths) {
  alpha <- path$alpha; lambdas <- unlist(path$lambdas)
  directory <- file.path(output,'corrected',family,paste0('alpha-',alpha)); result <- read_json(file.path(directory,'result.json'),simplifyVector=FALSE)
  verify_fit <- function(saved,xx,yy,ll) {
    shape <- unlist(saved$coefficient_shape); b <- read_f64(file.path(directory,saved$coefficient_file),shape)
    a <- unlist(saved$intercepts); stopifnot(identical(as.integer(shape),c(ncol(xx),length(ll))),all(is.finite(b)),all(is.finite(a)),all(unlist(saved$converged)))
    cert <- lapply(seq_along(ll),function(j) certificate(xx,yy,W,b[,j],a[j],ll[j],alpha,family))
    kkt <- vapply(cert,`[[`,numeric(1),'kkt_residual'); objective <- vapply(cert,`[[`,numeric(1),'objective')
    stopifnot(all(kkt<=config$thresh),max(abs(objective-unlist(saved$objective)))<1e-12,
              max(abs(kkt-unlist(saved$kkt_residual)))<1e-12)
    list(coefficients=b,intercepts=a,max_kkt=max(kkt),coefficient_entries=length(b),max_objective_difference=max(abs(objective-unlist(saved$objective))))
  }
  full <- verify_fit(result$full,X,y,lambdas); q <- length(lambdas); levels <- sort(unique(foldid)); losses <- matrix(NA_real_,length(levels),q)
  checks <- list(full=full[c('max_kkt','coefficient_entries','max_objective_difference')])
  for(k in seq_along(levels)) {
    train <- foldid!=levels[k]; yy <- if(family=='cox') y[train,,drop=FALSE] else y[train]
    fit <- verify_fit(result$folds[[k]],X[train,,drop=FALSE],yy,lambdas)
    checks[[paste0('fold_',k)]] <- fit[c('max_kkt','coefficient_entries','max_objective_difference')]
    losses[k,] <- held_out_loss(X,y,train,fit,family)
  }
  weights <- unlist(result$fold_weights); cvm <- colSums(losses*weights)/sum(weights)
  cvse <- sqrt(colSums(sweep(losses,2,cvm,'-')^2*weights)/sum(weights)/(length(levels)-1))
  stopifnot(max(abs(cvm-unlist(result$cvm)))<1e-12,max(abs(cvse-unlist(result$cvse)))<1e-12,
            all(unlist(result$valid)),lambdas[which.min(cvm)]==result$lambda_min)
  refit <- verify_fit(result$refit,X,y,result$lambda_min)
  checks$refit <- refit[c('max_kkt','coefficient_entries','max_objective_difference')]
  stopifnot(identical(which(refit$coefficients[,1]>0),as.integer(unlist(result$selection$plus))),
            identical(which(refit$coefficients[,1]<0),as.integer(unlist(result$selection$minus))))
  checks$CV_max_abs_error <- max(abs(cvm-unlist(result$cvm)),abs(cvse-unlist(result$cvse)))
  checks$lambda_min <- result$lambda_min; checks$selection_count <- sum(refit$coefficients!=0)
  verification$alphas[[as.character(alpha)]] <- checks
  cat('SAVED OUTPUTS REVERIFIED',family,alpha,'\n'); flush.console()
}
json(verification,file.path(output,paste0('saved-',family,'-reverification.json')))
