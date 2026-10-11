script <- normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]))
here <- dirname(script); root <- dirname(dirname(here))
Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
.libPaths(c(file.path(root,'reference/library'),.libPaths()))
suppressPackageStartupMessages(library(Matrix)); suppressPackageStartupMessages(library(jsonlite))
source(file.path(here,'oracle.R')); source(file.path(here,'compile_official.R')); compile_official(here,root)
results <- list()
files <- c('binomial_default_explicit_cv','binomial_tight_explicit_cv','binomial_default_with_trimming','binomial_default_automatic_cv')
for(name in files) {
  fixture <- read_json(file.path(root,'reference/fixtures',paste0(name,'.json')),simplifyVector=TRUE)
  input <- fixture$inputs; expected <- fixture$output
  if(!is.null(expected$ok)) { stopifnot(expected$ok); expected <- expected$result }
  got <- do.call(APML1,input)
  numeric_differences <- list()
  for(field in c('Beta','Beta0','fit','fit0','lambda.min','lambda.opt','flag')) {
    if(is.null(expected[[field]])) next
    a <- got[[field]]; b <- expected[[field]]
    if(is.data.frame(a) || is.data.frame(b)) {
      a <- a[,vapply(a,is.numeric,logical(1)),drop=FALSE]
      b <- b[,vapply(b,is.numeric,logical(1)),drop=FALSE]
    }
    a <- as.matrix(a); b <- as.matrix(b); stopifnot(identical(dim(a),dim(b)))
    difference <- max(abs(a-b)); stopifnot(is.finite(difference),difference<1e-10)
    numeric_differences[[field]] <- difference
  }
  results[[name]] <- numeric_differences
  cat('UNCHANGED -O2 GOLDEN VERIFIED',name,'max difference',max(unlist(numeric_differences)),'\n'); flush.console()
}
json(list(R=R.version.string,compiler_flags='-O2 -g0; no fast-math',source_patched=FALSE,fixtures=results),
     file.path(root,'validation-data/luad/r-oracle/official-O2-golden-verification.json'))
