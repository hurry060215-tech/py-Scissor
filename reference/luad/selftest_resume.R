# GPL-3.0-or-later. See ../../LICENSE.
# Small recovery-mechanism test; never counted as full-LUAD evidence.
script <- normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]))
here <- dirname(script); root <- dirname(dirname(here))
.libPaths(c(file.path(root,'reference/library'),.libPaths()))
suppressPackageStartupMessages(library(Matrix)); suppressPackageStartupMessages(library(jsonlite))
Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',R_MAKEVARS_USER=file.path(here,'Makevars'))
source(file.path(here,'oracle.R')); source(file.path(here,'resume_helpers.R'))
Rcpp::sourceCpp(file.path(here,'corrected_oracle.cpp'),cacheDir=file.path(here,'cache'),rebuild=FALSE)
base <- tempfile('resume-selftest-',tmpdir=file.path(root,'validation-data/luad/r-oracle'))
input <- file.path(base,'input'); output <- file.path(base,'output')
dir.create(file.path(input,'small'),recursive=TRUE); dir.create(output)
set.seed(1701); X <- matrix(rnorm(25*6),25,6); y <- 1+X[,1]-.4*X[,3]+.01*sin(seq_len(25))
A <- sparseMatrix(i=c(1:5,2:6),j=c(2:6,1:5),x=2,dims=c(6,6)); W <- normalize_adjacency(A)
foldid <- rep(1:5,5); alpha <- .2; lambdas <- c(.4,.2,.1,.04); tol <- 1e-8
write_f64(X,file.path(input,'small/X.bin')); writeMM(A,file.path(input,'graph.mtx'))
files <- lapply(c('small/X.bin','graph.mtx'),function(name) list(bytes=file.info(file.path(input,name))$size,sha256=sha256_file(file.path(input,name))))
names(files) <- c('small/X.bin','graph.mtx')
manifest <- list(cohorts=list(small=list(X_shape=dim(X),sample_ids=paste0('sample-',1:25),cell_ids=paste0('cell-',1:6))),
                 families=list(gaussian=list(cohort='small',y=y,foldid=foldid,paths=list(list(alpha=alpha,lambdas=lambdas)),thresh=tol,cutoff=.2)),files=files)
json(manifest,file.path(input,'manifest.json')); input_sha <- sha256_file(file.path(input,'manifest.json'))
contract <- verify_contract_files(input,'gaussian',input_sha)
contract$oracle_source_sha256 <- as.list(setNames(vapply(file.path(here,c('oracle.R','corrected_oracle.cpp')),sha256_file,character(1)),c('oracle.R','corrected_oracle.cpp')))
directory <- file.path(output,'gaussian/alpha-0.2'); dir.create(directory,recursive=TRUE)
full <- store_fit_atomic(fit_corrected(X,y,W,lambdas,alpha,'gaussian',tol),directory,'full')
train <- foldid!=1
fold1 <- store_fit_atomic(fit_corrected(X[train,,drop=FALSE],y[train],W,lambdas,alpha,'gaussian',tol),directory,'fold-1')
# Simulated interruption: completed full/fold 1, unfinished fold 2 binary with
# no completed checkpoint JSON. Recovery must overwrite only unfinished work.
writeBin(as.double(123),file.path(directory,'fold-2-coefficients.f64'))
protected <- c('full-coefficients.f64','full-checkpoint.json','fold-1-coefficients.f64','fold-1-checkpoint.json')
before <- vapply(file.path(directory,protected),sha256_file,character(1))
invoke <- function(mode,confirmed=FALSE) {
  arguments <- c(file.path(here,'resume_corrected.R'),'gaussian',input,output,mode,input_sha,if(confirmed) 'confirmed-stopped')
  log <- file.path(base,paste0(mode,if(confirmed) '-confirmed' else '','.log'))
  status <- system2(file.path(R.home('bin'),'Rscript'),shQuote(arguments),stdout=log,stderr=log)
  list(status=status,log=log)
}
validation <- invoke('validate'); stopifnot(validation$status==0)
summary <- read_json(file.path(output,'resume-audit/gaussian-validate-summary.json'),simplifyVector=FALSE)$alphas[['0.2']]
stopifnot(!summary$all_path_components_available,identical(unlist(summary$reused_components),c('full','fold-1')),
          summary$components[['fold-2']]=='missing_or_incomplete',length(summary$components)==6)
unconfirmed <- invoke('resume'); stopifnot(unconfirmed$status!=0)
resumed <- invoke('resume',TRUE); stopifnot(resumed$status==0)
after <- vapply(file.path(directory,protected),sha256_file,character(1)); stopifnot(identical(before,after))
summary <- read_json(file.path(output,'resume-audit/gaussian-resume-summary.json'),simplifyVector=FALSE)$alphas[['0.2']]
stopifnot(identical(unlist(summary$reused_components),c('full','fold-1')),
          identical(unlist(summary$new_components),c('fold-2','fold-3','fold-4','fold-5','refit')))
uninterrupted <- cv_corrected(X,y,W,lambdas,alpha,'gaussian',foldid,tol)
saved <- read_json(file.path(directory,'result.json'),simplifyVector=FALSE)
fits <- c(list(full=uninterrupted$full),setNames(uninterrupted$folds,paste0('fold-',1:5)),list(refit=uninterrupted$refit))
errors <- vapply(names(fits),function(role) {
  fit <- find_saved_fit(directory,role)$fit
  actual <- read_f64(file.path(directory,fit$coefficient_file),unlist(fit$coefficient_shape))
  max(abs(actual-fits[[role]]$coefficients))
},numeric(1))
stopifnot(all(errors==0),max(abs(unlist(saved$cvm)-uninterrupted$cvm))<1e-14,
          saved$lambda_min==uninterrupted$lambda_min,saved$lambda_1se==uninterrupted$lambda_1se)
expect_failure <- function(expression) {
  message <- tryCatch({force(expression); NA_character_},error=function(e) conditionMessage(e))
  stopifnot(!is.na(message)); message
}
parameters <- parameter_record(contract,'gaussian',1,alpha,lambdas,'full',seq_len(nrow(X)),tol)
audit <- file.path(output,'resume-audit/gaussian/alpha-0.2')
negative <- list()
negative$manifest_hash <- expect_failure(verify_contract_files(input,'gaussian',paste(rep('0',64),collapse='')))
bad <- find_saved_fit(directory,'full'); bad$fit$coefficient_shape <- c(6,3)
negative$shape <- expect_failure(verify_saved_fit(bad,directory,'full',X,y,W,lambdas,alpha,'gaussian',tol,parameters,audit))
bad <- find_saved_fit(directory,'full'); bad$fit$converged <- NULL
negative$missing_convergence <- expect_failure(verify_saved_fit(bad,directory,'full',X,y,W,lambdas,alpha,'gaussian',tol,parameters,audit))
changed <- parameter_record(contract,'gaussian',1,alpha,lambdas,'full',seq_len(nrow(X)),tol/100)
negative$threshold_signature <- expect_failure(verify_saved_fit(find_saved_fit(directory,'full'),directory,'full',X,y,W,lambdas,alpha,'gaussian',tol/100,changed,audit))
bad <- full$coefficients; bad[1,1] <- bad[1,1]+1
write_f64(bad,file.path(directory,'full-coefficients.f64'))
negative$corrupt_coefficients <- expect_failure(verify_saved_fit(find_saved_fit(directory,'full'),directory,'full',X,y,W,lambdas,alpha,'gaussian',tol,parameters,audit))
# Restore this test's intentional corruption so its retained outputs remain readable.
write_f64(full$coefficients,file.path(directory,'full-coefficients.f64'))
report <- list(status='passed',scope='small mechanism test; not LUAD numerical evidence',
               directory=base,simulated_interruption='full and fold 1 completed, fold 2 partial binary, remaining nodes absent',
               reused_components=summary$reused_components,new_components=summary$new_components,
               original_completed_files_unchanged=identical(before,after),max_coefficient_error_vs_uninterrupted=max(errors),
               missing_fold_not_accepted=TRUE,unconfirmed_resume_refused=TRUE,negative_rejections=negative,R=R.version.string)
json(report,file.path(root,'validation-data/luad/r-oracle/resume-selftest.json'))
cat('RESUME SELFTEST PASSED: unchanged completed files; all recomputed coefficients bit-identical; five integrity refusals\n')
