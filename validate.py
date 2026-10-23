"""Independent numerical checks and output invariants; no real data."""
import json
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import expit
from simulate import ROOT, design, fit_logit

r=ROOT/'results'; d=pd.read_csv(r/'synthetic_cohort.csv');d=d[d.retained]
X=design(d);y=d.preterm.to_numpy()
b,se,p=fit_logit(X,y)
def loss(v): return np.sum(np.logaddexp(0,X@v)-y*(X@v))
def grad(v):return X.T@(expit(X@v)-y)
opt=minimize(loss,np.zeros(X.shape[1]),jac=grad,method='BFGS',options={'gtol':1e-7})
assert np.max(np.abs(b-opt.x))<1e-5
assert np.max(np.abs(grad(opt.x)))<1e-4
power=pd.read_csv(r/'design_operating_characteristics.csv')
assert len(power)==24
assert (power.successful+power.failures==300).all()
assert (power.failures==0).all()
assert power.interaction_rejection.between(0,1).all()
assert power.interaction_coverage.between(0,1).all()
assert np.isfinite(power.select_dtypes('number')).all().all()
mol=pd.read_csv(r/'molecular_estimates.csv')
assert len(mol)==500 and mol.q.between(0,1).all()
assert (mol.true_beta>0).sum()==25
assert len(pd.read_csv(r/'synthetic_molecular_visits.csv'))==900
assert d.synthetic_id.is_unique
s=json.loads((r/'summary.json').read_text())
assert s['retained']==len(d) and s['preterm_events']==int(y.sum())
result={'status':'PASS','checks':['independent BFGS logistic agreement','24 scenarios x 300 repetitions','zero fit failures','finite estimates','500 features and BH bounds','900 molecular visits','summary consistency'],'maximum_coefficient_difference':float(np.max(np.abs(b-opt.x)))}
(r/'validation.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
