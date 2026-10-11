# GPL-3.0-or-later. Test the actual isolated driver's scheduling guard only.
script <- normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]))
here <- dirname(script); root <- dirname(dirname(here))
source(file.path(here,'resume_helpers.R'))
expressions <- parse(file.path(here,'resume_isolated_corrected.R'))
guard <- Filter(function(x) is.call(x) && identical(x[[1]],as.name('<-')) && identical(x[[2]],as.name('check_old_writer')),as.list(expressions))
stopifnot(length(guard)==1); eval(guard[[1]])
path <- tempfile('guard-',tmpdir=file.path(root,'validation-data/luad/r-oracle'))
writeLines('frozen old log',path)
snapshot <- list(source_log=path,source_log_sha256=sha256_file(path),source_log_mtime_ns=as.numeric(file.info(path)$mtime)*1e9)
check_old_writer()
expect_refusal <- function() {
  message <- tryCatch({check_old_writer();NA_character_},error=function(e) conditionMessage(e))
  stopifnot(!is.na(message),grepl('pausing new isolated scheduling',message,fixed=TRUE)); message
}
writeLines('advanced old log',path); changed_content <- expect_refusal()
snapshot$source_log_sha256 <- sha256_file(path); snapshot$source_log_mtime_ns <- as.numeric(file.info(path)$mtime)*1e9
check_old_writer(); Sys.setFileTime(path,as.POSIXct(snapshot$source_log_mtime_ns/1e9,origin='1970-01-01')+2)
changed_mtime <- expect_refusal(); unlink(path)
.libPaths(c(file.path(root,'reference/library'),.libPaths()))
jsonlite::write_json(list(status='passed',unchanged_source_accepted=TRUE,content_advance_refused=changed_content,
                         mtime_only_advance_refused=changed_mtime,scope='small scheduling guard test; no model fitting'),
                    file.path(root,'validation-data/luad/r-oracle/isolation-guard-selftest.json'),pretty=TRUE,auto_unbox=TRUE)
cat('ISOLATION GUARD PASSED: unchanged allowed; changed content and mtime both pause scheduling\n')
