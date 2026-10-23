"""Numerical, temporal, causal-algebra and sampler checks for the extension."""
import json
import numpy as np
import pandas as pd
from scipy.special import ndtr
from scipy.stats import norm
from extend_simulation import OUT, P, normal_precision

def main():
    tests=[];d=dict(np.load(OUT/'analysis_inputs.npz'));s=dict(np.load(OUT/'joint_population_draws.npz'))
    truth=json.loads((OUT/'generating_truth.json').read_text());rng=np.random.default_rng(777)
    long=pd.read_csv(OUT/'synthetic_longitudinal.csv');clinical=pd.read_csv(OUT/'synthetic_clinical.csv')
    assert len(long)==truth['n']*4 and len(clinical)==truth['n']
    assert (long.outcome_day>long.sample_day).all()
    merged=long[long.visit<=3].merge(clinical,on='id')
    assert (merged.outcome_day<merged.ascertainment_day).all()
    assert np.allclose(d['Z'][:,-1],d['mediator'][:,:3].mean(axis=1))
    tests.append('Temporal order and no postpartum leakage into complication predictors')
    # Independent analytic mediator posterior moment checks.
    X=d['X'].reshape(-1,P);m=d['mediator'].ravel();L=X.T@X+np.eye(P)
    V=np.linalg.inv(L);mu=V@X.T@m
    shape=3+len(m)/2;rate=2+.5*(m@m-mu@L@mu)
    drawmean=s['a'].mean(axis=(0,1));expectedvar=rate/(shape-1)*np.diag(V)
    mediator_z=(drawmean-mu)/np.sqrt(expectedvar/(s['a'].shape[0]*s['a'].shape[1]))
    assert np.max(np.abs(mediator_z))<5
    assert abs((s['sm']**2).mean()-rate/(shape-1))<.01
    tests.append('Mediator draws agree with independent normal-inverse-gamma posterior moments')
    # Precision sampler must return covariance P^-1, not P^-2.
    prec=np.array([[3.,.7],[.7,2.]]);rhs=np.array([.2,-.3]);vv=np.linalg.inv(prec)
    a=np.array([normal_precision(rng,prec,rhs) for _ in range(20000)])
    assert np.max(np.abs(a.mean(axis=0)-vv@rhs))<.02
    assert np.max(np.abs(np.cov(a.T)-vv))<.025
    tests.append('Gaussian precision sampling independently matches target moments')
    # Standardisation identity: total = direct + indirect for each posterior summary mean.
    ef=pd.read_csv(OUT/'joint_causal_effects.csv')
    for _,group in ef.groupby(['outcome','genotype']):
        z=group.set_index('effect')['mean'];assert abs(z['total']-z['direct']-z['indirect'])<1e-10
    sens=pd.read_csv(OUT/'joint_mediation_sensitivity.csv')
    for _,row in sens[np.isclose(sens.rho,0)].iterrows():
        base=ef[(ef.outcome==row.outcome)&(ef.genotype=='population')&(ef.effect=='indirect')].iloc[0]
        assert abs(base['mean']-row['mean'])<1e-10
    tests.append('Effect decomposition and rho=0 sensitivity reproduce base estimates')
    # The probit correlated-error sensitivity denominator is invariant algebraically.
    sm=.7/np.sqrt(3);bc=.35
    for rho in [-.6,-.3,0,.3,.6]:
        sy=1/np.sqrt(1-rho*rho);bs=bc-rho*sy/sm
        totalvar=bs*bs*sm*sm+sy*sy+2*bs*rho*sm*sy
        assert np.isclose(totalvar,1+bc*bc*sm*sm)
    tests.append('Probit correlated-error sensitivity variance identity')
    # Analytic marginal probit risk versus independent Monte Carlo integration.
    eta=.2;D=.49;load=.65;medcoef=.35;sm=.7/np.sqrt(3)
    analytic=ndtr(eta/np.sqrt(1+load**2*D+medcoef**2*sm**2))
    monte=ndtr(eta+load*rng.normal(0,np.sqrt(D),200000)+medcoef*rng.normal(0,sm,200000)).mean()
    assert abs(analytic-monte)<.003
    tests.append('Marginal probit integration verified by independent Monte Carlo')
    # Diagnostics cover all parameters and latent subject effects (no selective reporting).
    summaries={}
    for label in ['joint','wide_prior','noisy_module']:
        diag=pd.read_csv(OUT/f'{label}_diagnostics.csv')
        assert np.isfinite(diag.select_dtypes('number')).all().all()
        assert len(diag[diag.parameter.str.startswith('participant_')])==truth['n']*3
        summaries[label]={'max_rhat':float(diag.rhat.max()),'min_bulk_ess':float(diag.ess_bulk.min()),'min_tail_ess':float(diag.ess_tail.min())}
        assert diag.rhat.max()<1.01, f'{label}: Rhat threshold failed'
        assert diag.ess_bulk.min()>400 and diag.ess_tail.min()>400, f'{label}: ESS threshold failed'
    tests.append('All three four-chain fits pass Rhat<1.01 and bulk/tail ESS>400')
    scan=pd.read_csv(OUT/'transcriptomic_scan.csv');assert len(scan)==500 and scan.q.between(0,1).all()
    tests.append('Transcriptomic dimensions and FDR-adjustment bounds')
    (OUT/'validation.json').write_text(json.dumps({'status':'PASS','tests':tests,'diagnostics':summaries,'max_mediator_mean_MC_z':float(np.max(np.abs(mediator_z)))},indent=2)+'\n')
    print(json.dumps({'status':'PASS','tests':tests,'diagnostics':summaries},indent=2))

if __name__=='__main__':main()
