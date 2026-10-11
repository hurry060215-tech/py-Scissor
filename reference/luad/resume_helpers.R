# GPL-3.0-or-later. See ../../LICENSE.
# Checkpoint recovery helpers. Existing model coefficients are never trusted
# from a filename alone: input hashes, exact parameters and every KKT/objective
# certificate are independently checked before any reuse.
sha256_file <- function(path) {
  stopifnot(file.exists(path))
  value <- system2('sha256sum',shQuote(normalizePath(path)),stdout=TRUE)
  stopifnot(length(value)==1,is.null(attr(value,'status')))
  strsplit(value,' ',fixed=TRUE)[[1]][1]
}
sha256_raw <- function(bytes) {
  path <- tempfile(); on.exit(unlink(path)); con <- file(path,'wb'); writeBin(bytes,con); close(con)
  sha256_file(path)
}
atomic_json <- function(value,path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  temporary <- tempfile(pattern='.checkpoint-',tmpdir=dirname(path)); on.exit(unlink(temporary))
  json(value,temporary); stopifnot(file.rename(temporary,path))
}
parameter_record <- function(contract,family,path_index,alpha,lambdas,role,rows,tol) {
  numeric_hex <- paste(sprintf('%02x',as.integer(writeBin(as.double(c(alpha,lambdas,tol)),raw(),size=8,endian='little'))),collapse='')
  text <- paste(contract$manifest_sha256,contract$X_sha256,contract$graph_sha256,
                paste(unlist(contract$oracle_source_sha256),collapse=','),family,path_index,
                role,paste(rows,collapse=','),numeric_hex,sep='\n')
  list(parameter_sha256=sha256_raw(charToRaw(text)),manifest_sha256=contract$manifest_sha256,
       X_sha256=contract$X_sha256,graph_sha256=contract$graph_sha256,
       oracle_source_sha256=contract$oracle_source_sha256,family=family,
       alpha_path_index=path_index,alpha=alpha,lambdas=lambdas,fit_role=role,
       training_rows_1based=rows,thresh=tol,
       signature_format='UTF-8 identifiers/rows plus exact little-endian float64 alpha, lambdas and threshold')
}
verify_contract_files <- function(input,family,expected_manifest_sha256) {
  stopifnot(grepl('^[0-9a-f]{64}$',expected_manifest_sha256))
  actual <- sha256_file(file.path(input,'manifest.json'))
  if(actual!=expected_manifest_sha256) stop('Input manifest SHA-256 mismatch; checkpoint reuse refused')
  manifest <- jsonlite::read_json(file.path(input,'manifest.json'),simplifyVector=FALSE)
  config <- manifest$families[[family]]; stopifnot(!is.null(config))
  required <- c(file.path(config$cohort,'X.bin'),'graph.mtx')
  hashes <- vapply(required,function(name) {
    expected <- manifest$files[[name]]; stopifnot(!is.null(expected))
    path <- file.path(input,name); stopifnot(file.exists(path),file.info(path)$size==expected$bytes)
    hash <- sha256_file(path); if(hash!=expected$sha256) stop('Design/graph SHA-256 mismatch: ',name)
    hash
  },character(1))
  list(manifest=manifest,config=config,manifest_sha256=actual,X_sha256=unname(hashes[1]),graph_sha256=unname(hashes[2]))
}
find_saved_fit <- function(directory,role) {
  checkpoint <- file.path(directory,paste0(role,'-checkpoint.json'))
  if(file.exists(checkpoint)) return(list(path=checkpoint,fit=jsonlite::read_json(checkpoint,simplifyVector=FALSE),container=NULL))
  result_path <- file.path(directory,'result.json')
  if(!file.exists(result_path)) return(NULL)
  result <- jsonlite::read_json(result_path,simplifyVector=FALSE)
  fit <- if(role=='full') result$full else if(role=='refit') result$refit else result$folds[[as.integer(sub('fold-','',role,fixed=TRUE))]]
  if(is.null(fit)) return(NULL)
  list(path=result_path,fit=fit,container=result)
}
verify_saved_fit <- function(saved,directory,role,X,y,W,lambdas,alpha,family,tol,parameters,audit_directory) {
  if(is.null(saved)) return(NULL)
  json_hash <- sha256_file(saved$path); fit <- saved$fit
  expected_name <- paste0(role,'-coefficients.f64')
  stopifnot(identical(fit$coefficient_file,expected_name))
  binary <- file.path(directory,expected_name); shape <- as.integer(unlist(fit$coefficient_shape))
  if(!identical(shape,c(ncol(X),length(lambdas)))) stop('Checkpoint coefficient shape mismatch: ',role)
  if(!file.exists(binary) || file.info(binary)$size!=prod(shape)*8) stop('Incomplete checkpoint binary: ',role)
  binary_hash <- sha256_file(binary); B <- read_f64(binary,shape); intercepts <- as.numeric(unlist(fit$intercepts))
  stopifnot(length(intercepts)==length(lambdas),length(unlist(fit$converged))==length(lambdas),
            length(unlist(fit$objective))==length(lambdas),all(is.finite(B)),all(is.finite(intercepts)),
            all(unlist(fit$converged)))
  if(!is.null(saved$container)) {
    metadata <- saved$container$metadata
    if(metadata$family!=family || metadata$alpha!=alpha || metadata$thresh!=tol) stop('Completed-result parameters do not match requested family/alpha/threshold')
  }
  audit_path <- file.path(audit_directory,paste0(role,'.json'))
  previous <- if(file.exists(audit_path)) jsonlite::read_json(audit_path,simplifyVector=FALSE) else NULL
  if(!is.null(previous) && previous$parameters$parameter_sha256!=parameters$parameter_sha256) stop('Checkpoint parameter signature mismatch: ',role)
  if(!is.null(previous) && previous$coefficient_sha256!=binary_hash) stop('Previously verified checkpoint coefficient SHA-256 changed: ',role)
  certificates <- lapply(seq_along(lambdas),function(j) certificate(X,y,W,B[,j],intercepts[j],lambdas[j],alpha,family))
  residuals <- vapply(certificates,`[[`,numeric(1),'kkt_residual')
  objectives <- vapply(certificates,`[[`,numeric(1),'objective')
  if(any(!is.finite(residuals)) || any(residuals>tol)) stop('Checkpoint fails the requested absolute KKT threshold: ',role,' max=',max(residuals),' target=',tol)
  objective_error <- max(abs(objectives-as.numeric(unlist(fit$objective))))
  if(!is.finite(objective_error) || objective_error>1e-12*max(1,max(abs(objectives)))) stop('Checkpoint objective does not match the requested parameters: ',role)
  if(sha256_file(saved$path)!=json_hash || sha256_file(binary)!=binary_hash) stop('Checkpoint changed during verification; retry after its writer stops')
  audit <- list(status='verified',parameters=parameters,
    source_json=normalizePath(saved$path),source_json_sha256=json_hash,
    coefficient_file=normalizePath(binary),coefficient_sha256=binary_hash,coefficient_shape=shape,
    coefficient_entries=length(B),max_kkt=max(residuals),max_objective_error=objective_error,
    independently_recomputed_every_lambda=TRUE,legacy_checkpoint_enrolled_by_recomputation=is.null(previous),
    verified_at_utc=format(Sys.time(),tz='UTC',usetz=TRUE),R=R.version.string)
  atomic_json(audit,audit_path)
  fit$coefficients <- B; fit$intercepts <- intercepts
  for(name in c('objective','kkt_residual','outer_iterations','coordinate_sweeps')) fit[[name]] <- unlist(fit[[name]])
  fit$converged <- as.logical(unlist(fit$converged)); fit$termination_reason <- as.character(unlist(fit$termination_reason))
  fit$resume_verification <- audit
  fit
}
store_fit_atomic <- function(fit,directory,role) {
  dir.create(directory,recursive=TRUE,showWarnings=FALSE)
  name <- paste0(role,'-coefficients.f64'); temporary <- tempfile(pattern='.coefficients-',tmpdir=directory)
  on.exit(unlink(temporary)); write_f64(fit$coefficients,temporary)
  stopifnot(file.rename(temporary,file.path(directory,name)))
  fit$coefficient_file <- name; fit$coefficient_shape <- dim(fit$coefficients)
  saved <- fit; saved$coefficients <- NULL
  atomic_json(saved,file.path(directory,paste0(role,'-checkpoint.json')))
  fit
}
without_coefficients <- function(fit) { fit$coefficients <- NULL; fit }
