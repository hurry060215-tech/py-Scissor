import numpy as np
from scipy import sparse
from scipy.optimize import minimize
from scipy.special import expit,logsumexp
from scissor_py import fit_path,normalize_graph
bad=[]; count=0
for seed in range(10):
 r=np.random.default_rng(seed);n=30;p=4;X=r.normal(size=(n,p));X[:,1]=.8*X[:,0]+.2*X[:,1]
 G=sparse.diags([np.ones(p-1),np.ones(p-1)],[-1,1],format='csr');W=normalize_graph(G)
 outcomes={'gaussian':X[:,0]-.7*X[:,2]+r.normal(size=n)*.5,'binomial':r.binomial(1,expit(X[:,0]-.7*X[:,2])).astype(float),'cox':np.column_stack([r.integers(1,15,n),r.binomial(1,.8,n)])}
 for family,y in outcomes.items():
  for alpha in [0,.4,1]:
   a=fit_path(X,y,G,family=family,alpha=alpha,semantics='corrected',lambdas=[.1,.001],thresh=1e-9,maxit=10000)
   for k,lam in enumerate(a.lambdas):
    l1=lam*alpha;l2=lam*(1-alpha)
    def objective(z):
     beta=z[1:p+1]-z[p+1:];eta=X@beta+z[0];lap=beta-W@beta
     if family=='gaussian':loss=np.mean((eta-y)**2)/2;g=(eta-y)/n;grad=X.T@g;g0=g.sum()
     elif family=='binomial':loss=np.mean(np.logaddexp(0,(1-2*y)*eta));g=(expit(eta)-y)/n;grad=X.T@g;g0=g.sum()
     else:
      loss=0.;grad=np.zeros(p);g0=0.
      for t in np.unique(y[y[:,1]==1,0]):
       event=(y[:,0]==t)&(y[:,1]==1);risk=y[:,0]>=t;d=event.sum();ls=logsumexp(eta[risk]);loss+=d*ls-eta[event].sum();grad+=d*(np.exp(eta[risk]-ls)@X[risk])-X[event].sum(axis=0)
      loss/=n;grad/=n
     val=loss+l1*np.sum(z[1:])+l2/2*(beta@lap);grad+=l2*lap
     return val,np.r_[g0,grad+l1,-grad+l1]
    o=minimize(objective,np.zeros(1+2*p),jac=True,method='L-BFGS-B',bounds=[(None,None)]+[(0,None)]*(2*p),options={'ftol':1e-14,'gtol':1e-10,'maxiter':10000})
    b=a.coefficients[:,k];z=np.r_[a.intercepts[k],np.maximum(b,0),np.maximum(-b,0)];val=objective(z)[0];diff=val-o.fun
    if not a.converged[k] or diff>1e-7:bad.append((seed,family,alpha,lam,bool(a.converged[k]),diff,a.metadata['kkt_residual'][k],o.success))
    count+=1
print('checked',count,'bad',len(bad));print(*bad,sep='\n')
