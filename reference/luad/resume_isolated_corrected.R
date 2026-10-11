# GPL-3.0-or-later. See ../../LICENSE.
# Isolated recovery variant of frozen resume_corrected.R; original writer lifecycle stays unknown.
# Usage: Rscript resume_isolated_corrected.R FAMILY INPUT OUTPUT [validate|resume]
#        [EXPECTED_MANIFEST_SHA256] [confirmed-isolated-new-output]
# This separately named variant requires a hash-verified fresh-root snapshot.
# It never asserts the original writer stopped. New scheduling pauses if the old
# log advances; all writes remain confined to the independent snapshot root.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<3) stop('Need FAMILY INPUT OUTPUT; default mode is validate')
family <- args[1]; input <- normalizePath(args[2]); output <- args[3]
mode <- if(length(args)>=4) args[4] else 'validate'; stopifnot(mode %in% c('validate','resume'),family %in% c('gaussian','binomial','cox'))
expected_sha <- if(length(args)>=5) args[5] else '001f3cd90dd02953ed7379e01827699bc4da7a621ef055500289d4e908991eee'
if(mode=='resume' && (length(args)<6 || args[6]!='confirmed-isolated-new-output')) stop('Isolated resume requires confirmed-isolated-new-output; original source writer may remain unknown')
script <- normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]))
here <- dirname(script); root <- dirname(dirname(here))
Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',RCPP_PARALLEL_NUM_THREADS='1',R_MAKEVARS_USER=file.path(here,'Makevars'))
.libPaths(c(file.path(root,'reference/library'),.libPaths()))
suppressPackageStartupMessages(library(Matrix)); suppressPackageStartupMessages(library(jsonlite))
source(file.path(here,'oracle.R')); source(file.path(here,'resume_helpers.R'))
snapshot_path <- file.path(output,'isolation-snapshot.json')
snapshot <- jsonlite::read_json(snapshot_path,simplifyVector=FALSE)
stopifnot(normalizePath(output)==snapshot$output_root,normalizePath(output)!=snapshot$source_root,
          snapshot$input_manifest_sha256==expected_sha,snapshot$original_writer_lifecycle=='unknown_after_executor_key_change')
check_old_writer <- function() {
  changed_time <- abs(as.numeric(file.info(snapshot$source_log)$mtime)-snapshot$source_log_mtime_ns/1e9)>1e-6
  if(changed_time || sha256_file(snapshot$source_log)!=snapshot$source_log_sha256) stop('Old source log advanced: pausing new isolated scheduling; do not kill unknown original writer')
}
check_old_writer()
for(name in names(snapshot$files)) {
  record <- snapshot$files[[name]]; path <- file.path(output,name)
  stopifnot(file.exists(path),file.info(path)$size==record$bytes,sha256_file(path)==record$sha256)
}
contract <- verify_contract_files(input,family,expected_sha); manifest <- contract$manifest; config <- contract$config
contract$oracle_source_sha256 <- as.list(setNames(vapply(file.path(here,c('oracle.R','corrected_oracle.cpp')),sha256_file,character(1)),c('oracle.R','corrected_oracle.cpp')))
cohort <- manifest$cohorts[[config$cohort]]; X <- read_f64(file.path(input,config$cohort,'X.bin'),unlist(cohort$X_shape))
y <- if(family=='cox') do.call(rbind,config$y) else unlist(config$y); storage.mode(y) <- 'double'
foldid <- as.integer(unlist(config$foldid)); levels <- sort(unique(foldid))
W <- normalize_adjacency(as(as(readMM(file.path(input,'graph.mtx')),'dMatrix'),'generalMatrix'))
stopifnot(nrow(X)==length(foldid),ncol(X)==nrow(W),length(cohort$sample_ids)==nrow(X),length(cohort$cell_ids)==ncol(X))
if(mode=='resume') {
  Sys.setenv(SCISSOR_ORACLE_PROGRESS='1')
  Rcpp::sourceCpp(file.path(here,'corrected_oracle.cpp'),cacheDir=file.path(here,'cache'),rebuild=FALSE)
}
report <- list(mode=mode,family=family,manifest_sha256=expected_sha,alphas=list())
for(index in seq_along(config$paths)) {
  path <- config$paths[[index]]; alpha <- path$alpha; lambdas <- as.numeric(unlist(path$lambdas)); tol <- config$thresh
  directory <- file.path(output,family,paste0('alpha-',format(alpha,scientific=FALSE,trim=TRUE)))
  audit <- file.path(output,'resume-audit',family,paste0('alpha-',format(alpha,scientific=FALSE,trim=TRUE)))
  statuses <- list(); reused <- character(); computed <- character()
  obtain <- function(role,rows,ll) {
    yy <- if(family=='cox') y[rows,,drop=FALSE] else y[rows]
    parameters <- parameter_record(contract,family,index,alpha,ll,role,rows,tol)
    saved <- find_saved_fit(directory,role)
    fit <- verify_saved_fit(saved,directory,role,X[rows,,drop=FALSE],yy,W,ll,alpha,family,tol,parameters,audit)
    if(!is.null(fit)) { reused <<- c(reused,role); statuses[[role]] <<- 'verified_existing'; return(fit) }
    statuses[[role]] <<- 'missing_or_incomplete'
    if(mode=='validate') return(NULL)
    check_old_writer()
    cat('RESUME COMPUTING MISSING',family,alpha,role,'\n'); flush.console()
    fit <- fit_corrected(X[rows,,drop=FALSE],yy,W,ll,alpha,family,tol)
    fit <- store_fit_atomic(fit,directory,role)
    if(!all(fit$converged)) stop('Newly computed checkpoint did not converge; saved for diagnosis, no failed result accepted: ',role)
    fit <- verify_saved_fit(find_saved_fit(directory,role),directory,role,X[rows,,drop=FALSE],yy,W,ll,alpha,family,tol,parameters,audit)
    computed <<- c(computed,role); statuses[[role]] <<- 'computed_and_verified'
    fit
  }
  full <- obtain('full',seq_len(nrow(X)),lambdas); folds <- vector('list',length(levels))
  for(k in seq_along(levels)) folds[k] <- list(obtain(paste0('fold-',k),which(foldid!=levels[k]),lambdas))
  all_present <- !is.null(full) && all(vapply(folds,function(z) !is.null(z),logical(1)))
  if(all_present) {
    losses <- matrix(NA_real_,length(levels),length(lambdas)); weights <- numeric(length(levels))
    for(k in seq_along(levels)) {
      train <- foldid!=levels[k]; losses[k,] <- held_out_loss(X,y,train,folds[[k]],family)
      weights[k] <- if(family=='cox') sum(y[!train,2]) else sum(!train)
    }
    cvm <- colSums(losses*weights)/sum(weights); cvse <- sqrt(colSums(sweep(losses,2,cvm,'-')^2*weights)/sum(weights)/(length(levels)-1))
    best <- which.min(cvm); one_se <- which(cvm<=cvm[best]+cvse[best])[1]
    refit <- obtain('refit',seq_len(nrow(X)),lambdas[best])
    if(mode=='resume' && !is.null(refit) && length(computed)>0) {
      result <- list(full=without_coefficients(full),folds=lapply(folds,without_coefficients),refit=without_coefficients(refit),
        fold_levels=levels,fold_weights=weights,fold_losses=losses,cvm=cvm,cvse=cvse,valid=rep(TRUE,length(lambdas)),
        best_index_1based=best,one_se_index_1based=one_se,lambda_min=lambdas[best],lambda_1se=lambdas[one_se],
        selection=list(fraction=mean(refit$coefficients!=0),cutoff_reached=mean(refit$coefficients!=0)<config$cutoff,
          plus=which(refit$coefficients[,1]>0),minus=which(refit$coefficients[,1]<0),near_zero=which(refit$coefficients[,1]!=0 & abs(refit$coefficients[,1])<1e-5)),
        metadata=list(family=family,alpha=alpha,lambdas=lambdas,foldid=foldid,thresh=tol,cutoff=config$cutoff,
          semantics='corrected',n_samples=nrow(X),n_cells=ncol(X),sample_ids=cohort$sample_ids,cell_ids=cohort$cell_ids,
          no_sampling=TRUE,full_dataset=TRUE,input_manifest_sha256=expected_sha,R=R.version.string,
          resumed=TRUE,isolated_resume=TRUE,original_writer_lifecycle='unknown_after_executor_key_change',
          isolation_snapshot_sha256=sha256_file(snapshot_path),reused_components=reused,new_components=computed,
          oracle='same independent proximal-Newton coordinate-descent kernel; each reused node independently recertified'))
      atomic_json(result,file.path(directory,'result.json'))
    }
  }
  report$alphas[[as.character(alpha)]] <- list(components=statuses,all_path_components_available=all_present,reused_components=reused,new_components=computed)
  atomic_json(report,file.path(output,'resume-audit',paste0(family,'-',mode,'-summary.json')))
  cat('CHECKPOINT',toupper(mode),family,alpha,'verified',length(reused),'new',length(computed),'\n'); flush.console()
}
