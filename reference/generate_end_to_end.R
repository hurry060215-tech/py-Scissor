root <- dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[[1]])))
source(file.path(root,'compile.R'))
library(Matrix);library(jsonlite)
.libPaths(c(file.path(root,'library'),.libPaths()))
library(preprocessCore)
plain<-function(z){if(inherits(z,'Matrix'))return(as.matrix(z));if(is.data.frame(z))return(z);if(is.list(z))return(lapply(z,plain));z}
emit<-function(name,z){write_json(plain(z),file.path(root,'fixtures',paste0(name,'.json')),pretty=TRUE,auto_unbox=TRUE,digits=16,na='string',null='null',matrix='rowmajor');cat('SAVED',name,'\n')}
safe<-function(expr)tryCatch(list(ok=TRUE,result=plain(expr)),error=function(e)list(ok=FALSE,error=conditionMessage(e)))
# Auto-lambda counterparts of the already verified fixed-input CV cases.
for(family in c('gaussian','binomial')){
  d<-read_json(file.path(root,'fixtures',paste0(family,'_default_explicit_cv.json')),simplifyVector=TRUE)
  inp<-d$inputs;inp$lambda<-NULL;inp$nlambda<-20L
  emit(paste0(family,'_default_automatic_cv'),list(provenance=list(upstream_sha='311560a1f160f665917e7da0fcf289914ee99773',scope='APML1 automatic lambda with explicit folds; no L0 trimming'),inputs=inp,output=safe(do.call(APML1,inp))))
}
set.seed(1907)
ng<-32L;nb<-24L;nc<-6L
genes<-paste0('gene_',seq_len(ng));samples<-paste0('sample_',seq_len(nb));cells<-paste0('cell_',seq_len(nc))
all_expr<-matrix(rpois(ng*(nb+nc),lambda=rep(2+(seq_len(ng)%%7),nb+nc)),ng,nb+nc)
bulk<-rbind(all_expr[,seq_len(nb)],rep(11,nb));rownames(bulk)<-c(genes,'bulk_only');colnames(bulk)<-samples
sc<-rbind(all_expr[,nb+seq_len(nc)],rep(7,nc));rownames(sc)<-c(genes,'sc_only');colnames(sc)<-cells
sc<-sc[sample(seq_len(nrow(sc))),,drop=FALSE]
common<-intersect(rownames(bulk),rownames(sc));joint<-cbind(bulk[common,,drop=FALSE],sc[common,,drop=FALSE])
normalized<-normalize.quantiles(joint);X<-cor(normalized[,seq_len(nb)],normalized[,nb+seq_len(nc)])
graph<-matrix(0,nc,nc);graph[1,2]<-graph[2,1]<-graph[1,3]<-graph[3,1]<-graph[2,3]<-graph[3,2]<-graph[2,4]<-graph[4,2]<-1
ys<-list(gaussian=2+2*X[,1]-X[,3]+.1*cos(seq_len(nb)),binomial=as.numeric(2*X[,1]-X[,3]+.5*sin(seq_len(nb))>0),cox=cbind(time=rep(seq_len(nb/2),each=2),status=as.numeric(seq_len(nb)%%5!=0)))
foldid<-rep(1:4,length.out=nb);alphas<-c(.05,.2,.6);cutoff<-.4;nlambda<-20L
pheno_order<-c(24:13,1:12)
for(family in names(ys)){
  Y<-ys[[family]];history<-list();selected<-NULL;failure<-NULL
  for(alpha in alphas){
    cv<-safe(APML1(X,Y,family=family,penalty='Net',alpha=alpha,Omega=graph,nlambda=nlambda,foldid=foldid))
    if(!cv$ok){failure<-cv$error;history[[length(history)+1]]<-list(alpha=alpha,cv=cv);break}
    lam<-cv$result$lambda.min
    refit<-safe(APML1(X,Y,family=family,penalty='Net',alpha=alpha,Omega=graph,lambda=lam))
    if(!refit$ok){failure<-refit$error;history[[length(history)+1]]<-list(alpha=alpha,cv=cv,refit=refit);break}
    coefs<-as.numeric(refit$result$Beta);if(family=='binomial')coefs<-coefs[-1]
    frac<-mean(coefs!=0);labels<-ifelse(coefs>0,'Scissor+',ifelse(coefs<0,'Scissor-','Background'))
    selected<-list(alpha=alpha,lambda=lam,coefficients=coefs,labels=labels,selected_fraction=frac,cutoff_reached=frac<cutoff)
    history[[length(history)+1]]<-list(alpha=alpha,cv=cv,refit=refit,selected=selected)
    if(frac<cutoff)break
  }
  pheno<-if(is.matrix(Y))Y[pheno_order,,drop=FALSE] else Y[pheno_order]
  emit(paste0(family,'_expression_harness'),list(provenance=list(scope='Deterministic expression-to-selection harness following Scissor.R; explicit graph/folds bypass Seurat; not a whole Scissor() execution',Scissor_sha='311560a1f160f665917e7da0fcf289914ee99773',preprocessCore=as.character(packageVersion('preprocessCore')),R=R.version.string),inputs=list(bulk=bulk,bulk_gene_ids=rownames(bulk),sample_ids=samples,sc=sc,sc_gene_ids=rownames(sc),cell_ids=cells,phenotype=pheno,phenotype_sample_ids=samples[pheno_order],graph=graph,graph_cell_ids=cells,family=family,alpha=alphas,cutoff=cutoff,nlambda=nlambda,foldid=foldid),intermediates=list(common_gene_ids=common,joint=joint,normalized=normalized,similarity=X,phenotype_aligned=Y),alpha_history=history,result=selected,error=failure))
}
