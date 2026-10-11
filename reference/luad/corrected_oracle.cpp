// Independent exact-objective numerical oracle for Scissor Net validation.
// Copyright (c) 2026. GPL-3.0-or-later; see ../../LICENSE.
// This is not a patch or translation of the published Scissor optimizer.
// The R caller supplies centered X and a sparse symmetric normalized graph.
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <numeric>
#include <vector>
using namespace Rcpp;

namespace {
double soft(double z, double t) { return z > t ? z-t : (z < -t ? z+t : 0.0); }
double logadd(double a, double b) {
  if (!std::isfinite(a)) return b;
  if (!std::isfinite(b)) return a;
  double m=std::max(a,b); return m+std::log1p(std::exp(-std::abs(a-b)));
}
struct Stats { double loss; std::vector<double> g, h; };
Stats eta_statistics(const std::vector<double>& eta, const NumericVector& y,
            const NumericVector& time, int family, double divisor) {
  int n=eta.size(); Stats s{0.0,std::vector<double>(n),std::vector<double>(n)};
  if (family == 0) {
    for (int i=0;i<n;++i) { double r=eta[i]-y[i]; s.loss+=0.5*r*r/divisor;
      s.g[i]=r/divisor; s.h[i]=1.0/divisor; }
  } else if (family == 1) {
    for(int i=0;i<n;++i) { double e=eta[i];
      double mu=e>=0 ? 1.0/(1.0+std::exp(-e)) : std::exp(e)/(1.0+std::exp(e));
      s.loss+=(std::max(e,0.0)+std::log1p(std::exp(-std::abs(e)))-y[i]*e)/divisor;
      s.g[i]=(mu-y[i])/divisor; s.h[i]=std::max(mu*(1.0-mu),1e-12)/divisor; }
  } else {
    // Sorted ascending by time. Risk sets include every tied censor/event.
    std::vector<double> log_risk(n); double r=-INFINITY;
    for(int i=n-1;i>=0;--i) { r=logadd(r,eta[i]); log_risk[i]=r; }
    double log_hazard=-INFINITY;
    for(int i=0;i<n;) { int end=i+1; while(end<n && time[end]==time[i]) ++end;
      double events=0.0; for(int k=i;k<end;++k) events+=y[k];
      if(events>0) { log_hazard=logadd(log_hazard,std::log(events)-log_risk[i]);
        for(int k=i;k<end;++k) if(y[k]!=0) s.loss+=(log_risk[i]-eta[k])/divisor; }
      for(int k=i;k<end;++k) { double expected=std::isfinite(log_hazard) ? std::exp(eta[k]+log_hazard) : 0.0;
        s.g[k]=(expected-y[k])/divisor;
        // Diagonal risk contribution majorizes the current eta Hessian.
        s.h[k]=std::max(expected,1e-12)/divisor; }
      i=end;
    }
  }
  return s;
}
double adj(int j, const std::vector<double>& b, const IntegerVector& ptr,
           const IntegerVector& row, const NumericVector& val) {
  double z=0; for(int k=ptr[j];k<ptr[j+1];++k) z+=val[k]*b[row[k]]; return z;
}
double objective(const Stats& s, const std::vector<double>& b, double l1, double l2,
                 const IntegerVector& ptr,const IntegerVector& row,const NumericVector& val) {
  double out=s.loss; for(unsigned j=0;j<b.size();++j) out+=l1*std::abs(b[j])+0.5*l2*b[j]*(b[j]-adj(j,b,ptr,row,val));
  return out;
}
double kkt_one(double g,double b,double l1) { return b==0 ? std::max(std::abs(g)-l1,0.0) : std::abs(g+(b>0?l1:-l1)); }
}

// [[Rcpp::export]]
List corrected_oracle_path(NumericMatrix X, NumericVector y, NumericVector time,
    IntegerVector ptr, IntegerVector row, NumericVector val, NumericVector lambdas,
    double alpha, int family, double divisor, double tol=1e-8,
    int max_outer=2000, int max_inner=20000) {
  int n=X.nrow(),p=X.ncol(),nl=lambdas.size();
  if(y.size()!=n || ptr.size()!=p+1 || divisor<n || (family==2 && time.size()!=n)) stop("invalid oracle dimensions");
  NumericMatrix B(p,nl); NumericVector intercept(nl),obj(nl),errors(nl);
  IntegerVector outer(nl),sweeps(nl); LogicalVector converged(nl); CharacterVector reason(nl);
  std::vector<double> b(p,0),eta(n,0),grad(p),curvature(p);
  double a=0; if(family<2) { double mean=std::accumulate(y.begin(),y.end(),0.0)/n;
    a=family==0 ? mean : std::log(mean/(1-mean)); std::fill(eta.begin(),eta.end(),a); }
  for(int l=0;l<nl;++l) {
    double l1=lambdas[l]*alpha,l2=lambdas[l]*(1-alpha),err=INFINITY;
    std::string why="max_outer"; int total_sweeps=0,iteration=0;
    Stats current=eta_statistics(eta,y,time,family,divisor); double oldobj=objective(current,b,l1,l2,ptr,row,val);
    for(iteration=0;iteration<max_outer;++iteration) {
      checkUserInterrupt(); err=0; double intercept_grad=0;
      for(double z:current.g) intercept_grad+=z;
      if(family<2) err=std::abs(intercept_grad);
      for(int j=0;j<p;++j) { double g=0,h=0; const double* x=&X(0,j);
        for(int i=0;i<n;++i) { g+=x[i]*current.g[i]; h+=x[i]*x[i]*current.h[i]; }
        grad[j]=g+l2*(b[j]-adj(j,b,ptr,row,val)); curvature[j]=h+l2;
        err=std::max(err,kkt_one(grad[j],b[j],l1)); }
      if(err<=tol) { why="kkt"; break; }
      double qtol=std::max(tol*0.05,std::min(1e-4,err*0.03));
      std::vector<double> candidate=b,u=current.g; double ca=a;
      std::vector<int> active; std::vector<bool> admitted(p,false);
      for(int j=0;j<p;++j) if(b[j]!=0 || std::abs(grad[j])>l1) { active.push_back(j); admitted[j]=true; }
      double sumh=std::accumulate(current.h.begin(),current.h.end(),0.0);
      for(int sweep=0;sweep<max_inner;++sweep) {
        ++total_sweeps;
        if(family<2) { double gu=std::accumulate(u.begin(),u.end(),0.0),d=-gu/sumh; ca+=d;
          for(int i=0;i<n;++i) u[i]+=current.h[i]*d; }
        for(int j:active) { const double* x=&X(0,j); double g=0;
          for(int i=0;i<n;++i) g+=x[i]*u[i];
          g+=l2*(candidate[j]-adj(j,candidate,ptr,row,val));
          double next=curvature[j]>0 ? soft(candidate[j]-g/curvature[j],l1/curvature[j]) : 0;
          double d=next-candidate[j]; candidate[j]=next;
          if(d!=0) for(int i=0;i<n;++i) u[i]+=current.h[i]*x[i]*d;
        }
        if(sweep%5==4 || sweep==max_inner-1) {
          double qerr=family<2 ? std::abs(std::accumulate(u.begin(),u.end(),0.0)) : 0;
          for(int j=0;j<p;++j) { const double* x=&X(0,j); double g=0;
            for(int i=0;i<n;++i) g+=x[i]*u[i];
            g+=l2*(candidate[j]-adj(j,candidate,ptr,row,val));
            double e=kkt_one(g,candidate[j],l1); qerr=std::max(qerr,e);
            if(!admitted[j] && e>qtol) { active.push_back(j); admitted[j]=true; }
          }
          if(qerr<=qtol) break;
        }
      }
      std::vector<double> delta(p),deta(n,ca-a); double directional=(ca-a)*intercept_grad;
      for(int j=0;j<p;++j) { delta[j]=candidate[j]-b[j]; directional+=grad[j]*delta[j]+l1*(std::abs(candidate[j])-std::abs(b[j]));
        if(delta[j]!=0) { const double* x=&X(0,j); for(int i=0;i<n;++i) deta[i]+=x[i]*delta[j]; } }
      bool accepted=false; double step=1,newobj=INFINITY; Stats trial=current; std::vector<double> trialb(p),trialeta(n);
      for(int back=0;back<60;++back) {
        for(int j=0;j<p;++j) trialb[j]=step==1 ? candidate[j] : b[j]+step*delta[j];
        for(int i=0;i<n;++i) trialeta[i]=eta[i]+step*deta[i];
        trial=eta_statistics(trialeta,y,time,family,divisor); newobj=objective(trial,trialb,l1,l2,ptr,row,val);
        if(std::isfinite(newobj) && newobj<=oldobj+1e-4*step*std::min(0.0,directional)+2e-15*std::max(1.0,std::abs(oldobj))) { accepted=true; break; }
        step*=0.5;
      }
      if(!accepted) { why="line_search_failed"; break; }
      b=trialb; eta=trialeta; a+=step*(ca-a); current=trial; oldobj=newobj;
    }
    // Recompute certificate from returned coefficients, not stale inner state.
    err=family<2 ? std::abs(std::accumulate(current.g.begin(),current.g.end(),0.0)) : 0;
    for(int j=0;j<p;++j) { double g=0; const double* x=&X(0,j); for(int i=0;i<n;++i) g+=x[i]*current.g[i];
      err=std::max(err,kkt_one(g+l2*(b[j]-adj(j,b,ptr,row,val)),b[j],l1)); B(j,l)=b[j]; }
    converged[l]=err<=tol; if(converged[l]) why="kkt";
    intercept[l]=a; obj[l]=oldobj; errors[l]=err; outer[l]=iteration; sweeps[l]=total_sweeps; reason[l]=why;
    if(std::getenv("SCISSOR_ORACLE_PROGRESS")) Rcout << "path family=" << family << " n=" << n << " alpha=" << alpha << " lambda_index=" << l+1 << " kkt=" << err << " outer=" << iteration << " sweeps=" << total_sweeps << " reason=" << why << "\n";
  }
  return List::create(_["coefficients"]=B,_["intercepts"]=intercept,_["objective"]=obj,
      _["kkt_residual"]=errors,_["converged"]=converged,_["outer_iterations"]=outer,
      _["coordinate_sweeps"]=sweeps,_["termination_reason"]=reason);
}
