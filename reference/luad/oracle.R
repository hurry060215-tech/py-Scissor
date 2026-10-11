# Independent mathematical R oracle, GPL-3.0-or-later. See ../../LICENSE.
# All returned certificates are recomputed here, independently of C++ state.

read_f64 <- function(path, shape) {
  stopifnot(file.info(path)$size == prod(shape)*8)
  con <- file(path, 'rb'); on.exit(close(con))
  matrix(readBin(con, 'double', n=prod(shape), size=8, endian='little'), shape[1], shape[2])
}
write_f64 <- function(x, path) {
  con <- file(path, 'wb'); on.exit(close(con))
  writeBin(as.double(x), con, size=8, endian='little')
}
json <- function(x,path) jsonlite::write_json(x,path,pretty=TRUE,auto_unbox=TRUE,digits=16,na='null',null='null',matrix='rowmajor')

normalize_adjacency <- function(A) {
  A <- as(A,'dgCMatrix'); stopifnot(nrow(A)==ncol(A),all(is.finite(A@x)),all(A@x>=0))
  diag(A) <- 0; A <- Matrix::drop0(A); A@x[] <- 1
  stopifnot(isTRUE(all.equal(A,Matrix::t(A))))
  degree <- Matrix::rowSums(A); inv <- ifelse(degree>0,1/sqrt(pmax(degree,1)),0)
  as(Matrix::Diagonal(x=inv) %*% A %*% Matrix::Diagonal(x=inv),'dgCMatrix')
}

cox_reference <- function(X,y,beta,gradient=FALSE) {
  # Deliberately explicit risk sets, independent of C++ suffix-risk arithmetic.
  eta <- as.vector(X %*% beta); event_times <- sort(unique(y[y[,2]==1,1]))
  ll <- 0; score_eta <- y[,2]
  for(t in event_times) {
    deaths <- which(y[,1]==t & y[,2]==1); risk <- which(y[,1]>=t)
    e <- eta[risk]; m <- max(e); weight <- exp(e-m); total <- sum(weight)
    ll <- ll + sum(eta[deaths]-m-log(total))
    if(gradient) score_eta[risk] <- score_eta[risk]-length(deaths)*weight/total
  }
  if(gradient) list(loglik=ll,score_eta=score_eta) else ll
}

certificate <- function(X,y,W,b,a,lambda,alpha,family) {
  n <- nrow(X); eta <- as.vector(X %*% b+a); Z <- sweep(X,2,colMeans(X),'-')
  if(family=='gaussian') { error <- eta-y; loss <- mean(error^2)/2; g <- as.vector(crossprod(Z,error))/n; ia <- abs(mean(error)) }
  else if(family=='binomial') {
    error <- plogis(eta)-y
    loss <- mean(pmax(eta,0)+log1p(exp(-abs(eta)))-y*eta)
    g <- as.vector(crossprod(Z,error))/n; ia <- abs(mean(error))
  } else {
    keep <- y[,1]>=min(y[y[,2]==1,1]); Z <- sweep(X[keep,,drop=FALSE],2,colMeans(X[keep,,drop=FALSE]),'-')
    s <- cox_reference(X[keep,,drop=FALSE],y[keep,,drop=FALSE],b,TRUE)
    loss <- -s$loglik/n; g <- -as.vector(crossprod(Z,s$score_eta))/n; ia <- 0
  }
  lap <- as.vector(b-W%*%b); g <- g+lambda*(1-alpha)*lap
  residual <- pmax(abs(g)-lambda*alpha,0); active <- b!=0
  residual[active] <- abs(g[active]+lambda*alpha*sign(b[active]))
  list(objective=loss+lambda*alpha*sum(abs(b))+lambda*(1-alpha)*sum(b*lap)/2,
       kkt_residual=max(ia,residual),loss=loss,active=sum(active),
       plus=sum(b>0),minus=sum(b<0),near_zero=sum(active & abs(b)<1e-5))
}

fit_corrected <- function(X,y,W,lambdas,alpha,family,tol=1e-8) {
  W <- as(W,'dgCMatrix')
  original_n <- nrow(X); time <- numeric()
  if(family=='cox') {
    keep <- y[,1]>=min(y[y[,2]==1,1]); X0 <- X[keep,,drop=FALSE]; y0 <- y[keep,,drop=FALSE]
    order <- order(y0[,1],-y0[,2]); X0 <- X0[order,,drop=FALSE]; y0 <- y0[order,,drop=FALSE]
    time <- y0[,1]; response <- y0[,2]
  } else { X0 <- X; response <- y }
  means <- colMeans(X0); Z <- sweep(X0,2,means,'-')
  out <- corrected_oracle_path(Z,as.numeric(response),as.numeric(time),W@p,W@i,W@x,
      lambdas,alpha,match(family,c('gaussian','binomial','cox'))-1L,original_n,tol)
  if(family!='cox') out$intercepts <- out$intercepts-as.vector(crossprod(means,out$coefficients))
  else out$intercepts[] <- 0
  out$certificates <- lapply(seq_along(lambdas),function(j) certificate(X,y,W,out$coefficients[,j],out$intercepts[j],lambdas[j],alpha,family))
  cert_kkt <- vapply(out$certificates,`[[`,numeric(1),'kkt_residual')
  cert_obj <- vapply(out$certificates,`[[`,numeric(1),'objective')
  out$certificate_kkt_disagreement <- max(abs(cert_kkt-out$kkt_residual))
  out$certificate_objective_disagreement <- max(abs(cert_obj-out$objective))
  out$converged <- out$converged & cert_kkt<=tol
  out
}

held_out_loss <- function(X,y,training,fit,family) {
  held <- !training; B <- fit$coefficients; q <- ncol(B)
  if(family=='cox') return(vapply(seq_len(q),function(j)
    -(cox_reference(X,y,B[,j])-cox_reference(X[training,,drop=FALSE],y[training,,drop=FALSE],B[,j]))/sum(y[held,2]),numeric(1)))
  prediction <- sweep(X[held,,drop=FALSE]%*%B,2,fit$intercepts,'+')
  if(family=='gaussian') return(colMeans((prediction-y[held])^2))
  p <- pmin(pmax(plogis(prediction),1e-5),1-1e-5)
  -2*colMeans(y[held]*log(p)+(1-y[held])*log1p(-p))
}

cv_corrected <- function(X,y,W,lambdas,alpha,family,foldid,tol=1e-8,checkpoint_directory=NULL) {
  full_start <- proc.time()
  full <- fit_corrected(X,y,W,lambdas,alpha,family,tol)
  full$elapsed_seconds <- unname((proc.time()-full_start)['elapsed'])
  cat('corrected',family,'alpha',alpha,'FULL max KKT',max(full$kkt_residual),'outer',sum(full$outer_iterations),'sweeps',sum(full$coordinate_sweeps),'seconds',full$elapsed_seconds,'\n'); flush.console()
  if(!is.null(checkpoint_directory)) {
    dir.create(checkpoint_directory,recursive=TRUE,showWarnings=FALSE)
    json(export_fit(full,checkpoint_directory,'full'),file.path(checkpoint_directory,'full-checkpoint.json'))
  }
  levels <- sort(unique(foldid)); folds <- vector('list',length(levels)); losses <- matrix(NA_real_,length(levels),length(lambdas)); weights <- numeric(length(levels))
  for(i in seq_along(levels)) {
    training <- foldid!=levels[i]; yt <- if(family=='cox') y[training,,drop=FALSE] else y[training]
    folds[[i]] <- fit_corrected(X[training,,drop=FALSE],yt,W,lambdas,alpha,family,tol)
    if(!is.null(checkpoint_directory)) json(export_fit(folds[[i]],checkpoint_directory,paste0('fold-',i)),file.path(checkpoint_directory,paste0('fold-',i,'-checkpoint.json')))
    losses[i,] <- held_out_loss(X,y,training,folds[[i]],family)
    weights[i] <- if(family=='cox') sum(y[!training,2]) else sum(!training)
    cat('corrected',family,'alpha',alpha,'fold',levels[i],'max KKT',max(folds[[i]]$kkt_residual),'\n'); flush.console()
  }
  valid <- full$converged & apply(vapply(folds,`[[`,logical(length(lambdas)),'converged'),1,all)
  cvm <- colSums(losses*weights)/sum(weights)
  cvse <- sqrt(colSums(sweep(losses,2,cvm,'-')^2*weights)/sum(weights)/(length(levels)-1))
  best <- if(any(valid)) which.min(ifelse(valid,cvm,Inf)) else NA_integer_
  one_se <- if(any(valid)) which(valid & cvm<=cvm[best]+cvse[best])[1] else NA_integer_
  refit <- if(any(valid)) fit_corrected(X,y,W,lambdas[best],alpha,family,tol) else NULL
  list(full=full,folds=folds,refit=refit,fold_levels=levels,fold_weights=weights,fold_losses=losses,
       cvm=cvm,cvse=cvse,valid=valid,best_index_1based=best,one_se_index_1based=one_se,
       lambda_min=if(any(valid)) lambdas[best] else NA_real_,lambda_1se=if(any(valid)) lambdas[one_se] else NA_real_)
}

export_fit <- function(fit,directory,name) {
  if(is.null(fit)) return(NULL)
  file <- paste0(name,'-coefficients.f64'); write_f64(fit$coefficients,file.path(directory,file))
  fit$coefficient_file <- file; fit$coefficient_shape <- dim(fit$coefficients); fit$coefficients <- NULL
  fit
}
export_cv <- function(result,directory) {
  dir.create(directory,recursive=TRUE,showWarnings=FALSE)
  result$full <- export_fit(result$full,directory,'full')
  result$folds <- lapply(seq_along(result$folds),function(i) export_fit(result$folds[[i]],directory,paste0('fold-',i)))
  result$refit <- export_fit(result$refit,directory,'refit')
  json(result,file.path(directory,'result.json'))
}
