args <- commandArgs(trailingOnly=TRUE)
script <- normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]))
here <- dirname(script); root <- dirname(dirname(here))
input <- if(length(args)>=1) normalizePath(args[1]) else file.path(root,'validation-data/luad/prepared')
output <- if(length(args)>=2) args[2] else file.path(root,'validation-data/luad/r-oracle/corrected')
requested <- if(length(args)>=3) args[3] else 'all'
Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',RCPP_PARALLEL_NUM_THREADS='1',R_MAKEVARS_USER=file.path(here,'Makevars'))
Sys.setenv(SCISSOR_ORACLE_PROGRESS='1')
.libPaths(c(file.path(root,'reference/library'),.libPaths()))
suppressPackageStartupMessages(library(Matrix)); suppressPackageStartupMessages(library(jsonlite))
source(file.path(here,'oracle.R'))
dir.create(file.path(here,'cache'),recursive=TRUE,showWarnings=FALSE)
Rcpp::sourceCpp(file.path(here,'corrected_oracle.cpp'),cacheDir=file.path(here,'cache'),rebuild=FALSE)
manifest <- read_json(file.path(input,'manifest.json'),simplifyVector=FALSE)
A <- as(readMM(file.path(input,'graph.mtx')),'dgCMatrix'); W <- normalize_adjacency(A)
dir.create(output,recursive=TRUE,showWarnings=FALSE)
writeLines(capture.output(sessionInfo()),file.path(output,'sessionInfo.txt'))
results <- list()
for(family in names(manifest$families)) {
  if(requested!='all' && family!=requested) next
  config <- manifest$families[[family]]; cohort <- manifest$cohorts[[config$cohort]]
  X <- read_f64(file.path(input,config$cohort,'X.bin'),unlist(cohort$X_shape))
  y <- if(family=='cox') do.call(rbind,config$y) else unlist(config$y)
  storage.mode(y) <- 'double'; folds <- as.integer(unlist(config$foldid))
  stopifnot(nrow(X)==length(folds),ncol(X)==nrow(A),nrow(X)==length(cohort$sample_ids),ncol(X)==length(cohort$cell_ids))
  for(path in config$paths) {
    alpha <- path$alpha; lambdas <- as.numeric(unlist(path$lambdas)); tag <- paste0('alpha-',format(alpha,scientific=FALSE))
    destination <- file.path(output,family,tag)
    cat('START corrected',family,alpha,'samples',nrow(X),'cells',ncol(X),'lambdas',length(lambdas),'\n'); flush.console()
    begin <- proc.time()
    out <- cv_corrected(X,y,W,lambdas,alpha,family,folds,config$thresh,destination)
    out$metadata <- list(family=family,alpha=alpha,lambdas=lambdas,semantics='corrected',
       oracle='independent proximal Newton with exact-objective Armijo search and cyclic coordinate descent',
       n_samples=nrow(X),n_cells=ncol(X),sample_ids=cohort$sample_ids,cell_ids=cohort$cell_ids,
       foldid=folds,cutoff=config$cutoff,thresh=config$thresh,R=R.version.string,
       elapsed_seconds=unname((proc.time()-begin)['elapsed']),
       no_sampling=TRUE,full_dataset=TRUE,upstream_source_patched=FALSE,
       gaussian_response=if(family=='gaussian') 'synthetic known signal on full real LUAD similarity; no biological interpretation' else NULL)
    if(!is.null(out$refit)) out$selection <- list(fraction=mean(out$refit$coefficients!=0),
       cutoff_reached=mean(out$refit$coefficients!=0)<config$cutoff,
       plus=which(out$refit$coefficients[,1]>0),minus=which(out$refit$coefficients[,1]<0),
       near_zero=which(out$refit$coefficients[,1]!=0 & abs(out$refit$coefficients[,1])<1e-5))
    export_cv(out,destination)
    results[[length(results)+1]] <- list(family=family,alpha=alpha,path=file.path(family,tag,'result.json'),
        max_kkt=max(c(out$full$kkt_residual,unlist(lapply(out$folds,`[[`,'kkt_residual')))),
        valid_candidates=sum(out$valid),lambda_min=out$lambda_min,elapsed_seconds=out$metadata$elapsed_seconds)
    json(results,file.path(output,paste0('summary-',requested,'.json')))
    cat('DONE corrected',family,alpha,'valid',sum(out$valid),'seconds',out$metadata$elapsed_seconds,'\n'); flush.console()
    rm(out); gc()
  }
  rm(X); gc()
}
