# Full-data R preprocessing verification. No genes, cells, samples or entries are sampled.
args <- commandArgs(trailingOnly=TRUE)
script <- normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]))
here <- dirname(script); root <- dirname(dirname(here))
input <- if(length(args)>=1) normalizePath(args[1]) else file.path(root,'validation-data/luad/prepared')
output <- if(length(args)>=2) args[2] else file.path(root,'validation-data/luad/r-oracle/preprocessing')
original <- file.path(dirname(input),'original')
Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
.libPaths(c(file.path(root,'reference/library'),.libPaths()))
suppressPackageStartupMessages(library(Matrix)); suppressPackageStartupMessages(library(jsonlite)); library(preprocessCore)
source(file.path(here,'oracle.R'))
manifest <- read_json(file.path(input,'manifest.json'),simplifyVector=FALSE)
dir.create(output,recursive=TRUE,showWarnings=FALSE)
writeLines(capture.output(sessionInfo()),file.path(output,'sessionInfo.txt'))
read_object <- function(name) { e <- new.env(); load(file.path(original,name),e); stopifnot(length(ls(e))==1); e[[ls(e)[1]]] }
summary <- list(provenance=list(R=R.version.string,preprocessCore=as.character(packageVersion('preprocessCore')),
    preprocessCore_commit='5f4522fcbce0820b4983f453e3d56341da62da85',threading=FALSE,sampling='none',
    graph='uses the frozen explicit all-cell graph; no claim of Seurat/SNN equivalence'),cohorts=list())

# This is the published R character-index semantics for duplicate gene names.
# It is an explicit validation adapter, not a relaxation of the strict Python API.
sc <- read_object('scRNA-seq.RData'); stopifnot(identical(dim(sc),c(33694L,4102L)))
stopifnot(!anyDuplicated(rownames(sc)),!anyDuplicated(colnames(sc)))
library_size <- colSums(sc); stopifnot(all(library_size>0))
for(start in seq(1,ncol(sc),by=32)) { cols <- start:min(start+31,ncol(sc)); sc[,cols] <- log1p(sweep(sc[,cols,drop=FALSE],2,library_size[cols],'/')*10000) }
for(cohort_name in names(manifest$cohorts)) {
  config <- manifest$cohorts[[cohort_name]]; bulk <- read_object(paste0('TCGA_LUAD_',cohort_name,'.RData'))
  common <- intersect(rownames(bulk),rownames(sc)); nb <- ncol(bulk); nc <- ncol(sc)
  stopifnot(identical(common,unlist(config$common_gene_ids)),identical(colnames(bulk),unlist(config$sample_ids)),
            identical(colnames(sc),unlist(config$cell_ids)),sum(duplicated(rownames(bulk)))==config$duplicate_rows)
  first_indices <- match(common,rownames(bulk)); stopifnot(identical(unname(bulk[common,1]),unname(bulk[first_indices,1])))
  con <- file(file.path(input,cohort_name,'joint.bin'),'rb'); max_error <- 0; comparisons <- 0
  for(start in seq(1,nb+nc,by=32)) {
    cols <- start:min(start+31,nb+nc)
    expected <- matrix(readBin(con,'double',n=length(common)*length(cols),size=8,endian='little'),length(common),length(cols))
    actual <- matrix(0,length(common),length(cols)); bc <- cols<=nb; cc <- !bc
    if(any(bc)) actual[,bc] <- bulk[common,cols[bc],drop=FALSE]
    if(any(cc)) actual[,cc] <- sc[common,cols[cc]-nb,drop=FALSE]
    max_error <- max(max_error,max(abs(actual-expected))); comparisons <- comparisons+length(actual)
  }
  stopifnot(length(readBin(con,'raw',n=1))==0); close(con)
  summary$cohorts[[cohort_name]] <- list(raw_bulk_shape=dim(bulk),raw_sc_shape=dim(sc),
    bulk_duplicate_rows=sum(duplicated(rownames(bulk))),shared_genes=length(common),
    duplicate_policy='unique intersect in bulk order; R character indexing selects first source row',
    raw_to_joint_entries_compared=comparisons,raw_to_joint_max_abs_error=max_error)
  stopifnot(max_error<1e-12,comparisons==prod(unlist(config$joint_shape)))
  rm(bulk,actual,expected); gc(); cat('RAW INPUTS VERIFIED',cohort_name,'entries',comparisons,'max error',max_error,'\n'); flush.console()
}
rm(sc); gc()
survival <- read_object('TCGA_LUAD_survival.RData'); tp53 <- read_object('TCGA_LUAD_TP53_mutation.RData')
stopifnot(identical(as.character(survival$TCGA_patient_barcode),unlist(manifest$cohorts$exp1$sample_ids)),
          identical(as.numeric(survival$OS_time),vapply(manifest$families$cox$y,`[[`,numeric(1),1L)),
          identical(as.numeric(survival$Status),vapply(manifest$families$cox$y,`[[`,numeric(1),2L)),
          identical(as.numeric(tp53),as.numeric(unlist(manifest$families$binomial$y))))
summary$phenotype_source_verification <- TRUE

for(cohort_name in names(manifest$cohorts)) {
  config <- manifest$cohorts[[cohort_name]]; shape <- unlist(config$joint_shape); nb <- unlist(config$X_shape)[1]; nc <- unlist(config$X_shape)[2]
  begin <- proc.time(); joint <- read_f64(file.path(input,cohort_name,'joint.bin'),shape)
  cat('NORMALIZING ALL VALUES',cohort_name,shape,'\n'); flush.console()
  # copy=FALSE performs the identical preprocessCore computation in place.
  joint <- preprocessCore::normalize.quantiles(joint,copy=FALSE)
  con <- file(file.path(input,cohort_name,'normalized.bin'),'rb'); max_error <- 0; total_error <- 0; above_1e8 <- 0; comparisons <- 0
  for(start in seq(1,ncol(joint),by=32)) {
    cols <- start:min(start+31,ncol(joint)); expected <- matrix(readBin(con,'double',n=shape[1]*length(cols),size=8,endian='little'),shape[1],length(cols))
    delta <- abs(joint[,cols,drop=FALSE]-expected); max_error <- max(max_error,max(delta)); total_error <- total_error+sum(delta)
    above_1e8 <- above_1e8+sum(delta>1e-8); comparisons <- comparisons+length(delta)
  }
  stopifnot(length(readBin(con,'raw',n=1))==0); close(con); rm(expected,delta); gc()
  similarity <- matrix(0,nb,nc); bulk <- joint[,seq_len(nb),drop=FALSE]
  for(start in seq(1,nc,by=32)) {
    cols <- start:min(start+31,nc); similarity[,cols] <- cor(bulk,joint[,nb+cols,drop=FALSE])
  }
  expected_X <- read_f64(file.path(input,cohort_name,'X.bin'),c(nb,nc)); x_error <- max(abs(similarity-expected_X))
  write_f64(similarity,file.path(output,paste0(cohort_name,'-R-X.bin')))
  summary$cohorts[[cohort_name]]$quantile <- list(entries_compared=comparisons,max_abs_error=max_error,
       mean_abs_error=total_error/comparisons,entries_above_1e8=above_1e8)
  summary$cohorts[[cohort_name]]$pearson <- list(entries_compared=length(similarity),max_abs_error=x_error,shape=dim(similarity),
       all_finite=all(is.finite(similarity)),R_X_file=paste0(cohort_name,'-R-X.bin'))
  summary$cohorts[[cohort_name]]$elapsed_seconds <- unname((proc.time()-begin)['elapsed'])
  json(summary,file.path(output,'result.json'))
  stopifnot(comparisons==prod(shape),max_error<1e-8,x_error<1e-12,all(is.finite(similarity)))
  rm(joint,bulk,similarity,expected_X); gc()
  cat('FULL PREPROCESSING VERIFIED',cohort_name,'qnorm max error',max_error,'Pearson max error',x_error,'\n'); flush.console()
}
