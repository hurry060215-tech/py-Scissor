# Compile unchanged original sources; only compiler optimization differs from
# the earlier small-fixture -O0 build. No fast-math and no upstream source edits.
compile_official <- function(here,root) {
  Sys.setenv(R_MAKEVARS_USER=file.path(here,'Makevars'))
  upstream <- file.path(root,'reference/upstream')
  cache <- file.path(here,'cache-upstream-O2'); dir.create(cache,recursive=TRUE,showWarnings=FALSE)
  Rcpp::sourceCpp(file.path(upstream,'src/Scissor.cpp'),cacheDir=cache,rebuild=FALSE,env=.GlobalEnv)
  for(f in c('APML1.R','net_lm.R','net_logit.R','net_cox.R')) source(file.path(upstream,'R',f),local=.GlobalEnv)
}
log_kernel_calls <- function() {
  # Generated Rcpp entrypoints are wrapped externally for logging only. Every
  # argument and the complete returned object are passed through unchanged.
  logged <- function(fn,label) {
    force(fn); force(label)
    function(...) {
      cat('ENTER',label,format(Sys.time(),tz='UTC',usetz=TRUE),'\n'); flush.console(); start <- proc.time()
      result <- fn(...)
      cat('EXIT',label,'seconds',unname((proc.time()-start)['elapsed']),'path_length',result$nlambda,
          'flags',paste(result$flag,collapse=','),'\n'); flush.console()
      result
    }
  }
  for(name in c('NetLogC','cvNetLogC')) assign(name,logged(get(name,.GlobalEnv),name),envir=.GlobalEnv)
}
