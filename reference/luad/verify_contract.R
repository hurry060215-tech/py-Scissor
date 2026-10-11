# Independent null-score/lambda-path and fixed-fold checks on complete designs.
script <- normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]))
here <- dirname(script); root <- dirname(dirname(here))
Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
.libPaths(c(file.path(root,'reference/library'),.libPaths()))
suppressPackageStartupMessages(library(Matrix)); suppressPackageStartupMessages(library(jsonlite))
source(file.path(here,'oracle.R'))
input <- file.path(root,'validation-data/luad/prepared')
manifest <- read_json(file.path(input,'manifest.json'),simplifyVector=FALSE)
results <- list(R=R.version.string,sampling='none',families=list())
for(family in names(manifest$families)) {
  config <- manifest$families[[family]]; cohort <- manifest$cohorts[[config$cohort]]
  X <- read_f64(file.path(input,config$cohort,'X.bin'),unlist(cohort$X_shape))
  y <- if(family=='cox') do.call(rbind,config$y) else unlist(config$y); storage.mode(y) <- 'double'
  folds <- unlist(config$foldid); stopifnot(all(is.finite(X)),all(is.finite(y)),identical(as.integer(sort(unique(folds))),1:5))
  if(family=='cox') {
    s <- cox_reference(X,y,rep(0,ncol(X)),TRUE)
    score <- as.vector(crossprod(X,s$score_eta))/nrow(X)
  } else score <- as.vector(crossprod(sweep(X,2,colMeans(X),'-'),y-mean(y)))/nrow(X)
  paths <- list()
  for(path in config$paths) {
    maximum <- max(abs(score))/path$alpha
    expected <- maximum*0.01^((0:19)/19); actual <- unlist(path$lambdas)
    error <- max(abs(expected-actual)); stopifnot(length(actual)==20,all(diff(actual)<0),error<1e-12)
    paths[[as.character(path$alpha)]] <- list(independent_lambda_max=maximum,frozen_lambda_max=actual[1],full_path_max_abs_error=error)
  }
  fold_summary <- lapply(1:5,function(k) {
    train <- folds!=k; stopifnot(sum(train)>=2)
    out <- list(fold=k,training_samples=sum(train),validation_samples=sum(!train))
    if(family=='cox') {
      out$training_events <- sum(y[train,2]); out$validation_events <- sum(y[!train,2]); stopifnot(out$training_events>0,out$validation_events>0)
    } else if(family=='binomial') {
      out$training_classes <- as.list(table(y[train])); out$validation_classes <- as.list(table(y[!train])); stopifnot(length(out$training_classes)==2,length(out$validation_classes)==2)
    }
    out
  })
  result <- list(samples=nrow(X),cells=ncol(X),paths=paths,folds=fold_summary)
  if(family=='gaussian') {
    expected_y <- 2+2*X[,1]-1.5*X[,138]+X[,778]+0.02*sin(seq_len(nrow(X)))
    result$known_signal_max_abs_error <- max(abs(y-expected_y)); stopifnot(result$known_signal_max_abs_error<1e-14)
  }
  results$families[[family]] <- result
  cat('FULL FIXED CONTRACT VERIFIED',family,'\n'); flush.console()
}
json(results,file.path(root,'validation-data/luad/r-oracle/fixed-contract-verification.json'))
