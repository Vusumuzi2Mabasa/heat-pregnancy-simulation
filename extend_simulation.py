"""Purpose-built synthetic joint model. No real patients, genes or climate data.

Blocked Gibbs sampling: correlated Gaussian random intercepts link three
longitudinal outcomes; an Albert-Chib probit endpoint shares the vascular
random intercept. A measured synthetic expression module precedes outcomes.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import argparse, json, time, platform
import numpy as np
import pandas as pd
from scipy.special import ndtr, ndtri, expit
from scipy.stats import invwishart, t as student_t
from scipy.linalg import cho_factor, cho_solve, solve_triangular
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results'/'extended'; OUT.mkdir(parents=True,exist_ok=True)
SEED=2026102302
DOMAINS=['vascular','depression','anxiety']
FEATURES=['intercept','heat','genotype','heat_x_genotype','age','deprivation',
          'ancestry','baseline_distress','site_2','site_3','season_sin','season_cos',
          'visit_2','visit_3','postpartum']
P=len(FEATURES)

def normal_precision(rng, precision, rhs):
    L=np.linalg.cholesky(precision)
    mean=cho_solve((L,True),rhs)
    return mean+solve_triangular(L.T,rng.normal(size=len(rhs)),lower=False)

def create_data(n=600,seed=SEED):
    rng=np.random.default_rng(seed)
    weather=pd.read_csv(ROOT/'results'/'synthetic_daily_weather.csv.gz',parse_dates=['date'])
    lookup={key:g.sort_values('date').tmax_C.to_numpy() for key,g in weather.groupby(['site','neighbourhood'])}
    site=np.arange(n)%3; rng.shuffle(site); hood=rng.integers(0,6,n)
    concept=rng.integers(0,730,n)+int((pd.Timestamp('2027-01-01')-pd.Timestamp('2010-01-01')).days)
    # Fixed term delivery at day 280; postpartum sample at day 322. No survival model.
    sample_days=np.array([98,140,182,322]); dates=concept[:,None]+sample_days
    temp=np.zeros((n,4))
    for i in range(n):
        for v in range(4):temp[i,v]=lookup[(site[i],hood[i])][dates[i,v]-7:dates[i,v]].mean()
    hs=(temp-temp.mean())/temp.std()
    pc=rng.normal(.3*(site-1),1,n)
    G=rng.binomial(2,expit(np.log(.3/.7)+.2*pc)); gc=G-.6
    age=rng.normal(size=n);dep=rng.normal(.25*site,1,n)
    baseline=rng.normal(.25*dep,1,n)
    X=np.zeros((n,4,P));X[:,:,0]=1;X[:,:,1]=hs;X[:,:,2]=gc[:,None];X[:,:,3]=hs*gc[:,None]
    for j,a in [(4,age),(5,dep),(6,pc),(7,baseline),(8,site==1),(9,site==2)]:X[:,:,j]=a[:,None]
    X[:,:,10]=np.sin(2*np.pi*dates/365.25);X[:,:,11]=np.cos(2*np.pi*dates/365.25)
    X[:,1,12]=1;X[:,2,13]=1;X[:,3,14]=1
    a=np.zeros(P);a[[1,2,3,5,6,7,12,13,14]]=[.4,.10,.12,.12,.08,.10,.05,.10,.08]
    sm=.7; med=np.einsum('nvp,p->nv',X,a)+rng.normal(0,sm,(n,4))
    B=np.zeros((P+1,3))
    B[1]=[.18,.15,.12];B[2]=[.08,.06,.05];B[3]=[.08,.06,.05]
    B[4]=[.12,.04,.03];B[5]=[.10,.16,.15];B[6]=[.06,.05,.05];B[7]=[.10,.40,.35]
    B[8]=[.05,.04,.03];B[9]=[.10,.08,.06];B[10]=[.06,.07,.06]
    B[12]=[.08,.05,.04];B[13]=[.15,.10,.08];B[14]=[-.08,.18,.14]
    B[-1]=[.30,.40,.35]
    D=np.array([[.49,.21,.14],[.21,.64,.32],[.14,.32,.49]])
    b=rng.multivariate_normal(np.zeros(3),D,n);sd=np.array([.65,.80,.75])
    W=np.concatenate([X,med[:,:,None]],axis=2)
    y=np.einsum('nvp,pk->nvk',W,B)+b[:,None,:]+rng.normal(size=(n,4,3))*sd
    # Clinical endpoint follows all three antenatal pairs, never uses postpartum values.
    Z=np.concatenate([X[:,:3].mean(axis=1),med[:,:3].mean(axis=1)[:,None]],axis=1)
    # visit indicators are constant after averaging and redundant with intercept; omit them.
    clinical_cols=list(range(12))+[P]
    Z=Z[:,clinical_cols]
    c=np.zeros(Z.shape[1]);c[0]=-1.35;c[1]=.20;c[2]=.06;c[3]=.10
    c[4]=.15;c[5]=.12;c[6]=.06;c[7]=.08;c[-1]=.35
    loading=.65
    latent=Z@c+loading*b[:,0]+rng.normal(size=n);event=(latent>0).astype(int)
    ids=np.array([f'EXT-{i+1:04}' for i in range(n)])
    long=pd.DataFrame({'id':np.repeat(ids,4),'site':np.repeat(site,4),'neighbourhood':np.repeat(hood,4),
         'visit':np.tile(np.arange(1,5),n),'sample_gestational_week':np.tile([14,20,26,46],n),
         'sample_day':dates.ravel(),'outcome_day':(dates+7).ravel(),'temperature_C':temp.ravel(),
         'module':med.ravel(),'vascular_z':y[:,:,0].ravel(),'depression_z':y[:,:,1].ravel(),
         'anxiety_z':y[:,:,2].ravel()})
    for j,name in enumerate(FEATURES):long[name]=X[:,:,j].ravel()
    long.to_csv(OUT/'synthetic_longitudinal.csv',index=False)
    pd.DataFrame({'id':ids,'complication':event,'ascertainment_day':concept+280,'G':G,
                  'vascular_intercept_truth':b[:,0]}).to_csv(OUT/'synthetic_clinical.csv',index=False)
    truth=dict(seed=seed,n=n,visits=4,features=FEATURES,domains=DOMAINS,a=a.tolist(),B=B.tolist(),
       D=D.tolist(),residual_sd=sd.tolist(),mediator_sd=sm,clinical_c=c.tolist(),loading=loading,
       clinical_columns=clinical_cols,temperature_mean=float(temp.mean()),temperature_sd=float(temp.std()),
       complications=int(event.sum()),statement='All values and effects are synthetic assumptions; fixed term delivery.')
    (OUT/'generating_truth.json').write_text(json.dumps(truth,indent=2)+'\n')
    np.savez_compressed(OUT/'analysis_inputs.npz',X=X,W=W,Z=Z,y=y,mediator=med,event=event)
    return dict(X=X,W=W,Z=Z,y=y,mediator=med,event=event),truth

def run_chain(data,chain,draws,warmup,prior_scale=1.,independent=False):
    """Exact conditional Gibbs updates; all reported iterations are post warmup."""
    rng=np.random.default_rng(SEED+1000+chain+int(prior_scale*100))
    X,W,Z,y,med,event=[data[k] for k in ['X','W','Z','y','mediator','event']]
    n,V,K=y.shape;p=W.shape[-1];q=Z.shape[1]
    xf=X.reshape(-1,P);wf=W.reshape(-1,p);yf=y.reshape(-1,K)
    xtx=xf.T@xf;wtw=wf.T@wf;z0=Z.T@Z
    priorprec=np.eye(p)/prior_scale**2;priorprec[0,0]=1/(2*prior_scale)**2
    cprec=np.eye(q+1)/prior_scale**2;cprec[0,0]=1/(2*prior_scale)**2
    # Conjugate NIG mediator posterior. No mediator random intercept in the DGP.
    Lm=xtx+np.eye(P)/prior_scale**2
    Vm=np.linalg.inv(Lm);am=Vm@xf.T@med.ravel();Lmchol=np.linalg.cholesky(Vm)
    am_shape=3+n*V/2
    am_rate=2+.5*(med.ravel()@med.ravel()-am@Lm@am)
    B=rng.normal(0,.2,(p,K));D=np.eye(K)*(.25+chain*.20)
    sig=np.ones(K);b=rng.normal(0,.15,(n,K));c=np.zeros(q);c[0]=-1
    loading=.5;z=np.where(event,1.,-1.)
    storage={ 'B':np.zeros((draws,p,K)), 'D':np.zeros((draws,K,K)), 'sd':np.zeros((draws,K)),
      'a':np.zeros((draws,P)), 'sm':np.zeros(draws), 'c':np.zeros((draws,q)),
      'loading':np.zeros(draws),'b':np.zeros((draws,n,K),dtype=np.float32)}
    for it in range(warmup+draws):
        # Complete-data Gaussian regression for longitudinal outcomes.
        adjusted=(y-b[:,None,:]).reshape(-1,K)
        for k in range(K):B[:,k]=normal_precision(rng,wtw/sig[k]+priorprec,wf.T@adjusted[:,k]/sig[k])
        # Probit latent normal under observed sign constraints.
        eta=Z@c+loading*b[:,0];cut=ndtr(-eta)
        u=rng.random(n);quant=np.where(event,cut+(1-cut)*u,cut*u)
        z=eta+ndtri(np.clip(quant,1e-14,1-1e-14))
        Cz=np.column_stack([Z,b[:,0]])
        cc=normal_precision(rng,Cz.T@Cz+cprec,Cz.T@z);c=cc[:-1];loading=cc[-1]
        # Event likelihood updates the vascular random intercept: genuine joint fit.
        prec=np.linalg.inv(D)+np.diag(V/sig);prec[0,0]+=loading**2
        cov=np.linalg.inv(prec)
        rhs=(y-np.einsum('nvp,pk->nvk',W,B)).sum(axis=1)/sig
        rhs[:,0]+=loading*(z-Z@c)
        b=rhs@cov+rng.normal(size=(n,K))@np.linalg.cholesky(cov).T
        if independent:
            vv=(2+.5*(b*b).sum(axis=0))/rng.gamma(3+n/2,size=K);D=np.diag(vv)
        else:D=invwishart.rvs(df=K+2+n,scale=np.eye(K)*.5+b.T@b,random_state=rng)
        residual=y-np.einsum('nvp,pk->nvk',W,B)-b[:,None,:]
        sig=(2+.5*(residual**2).sum(axis=(0,1)))/rng.gamma(3+n*V/2,size=K)
        sm2=am_rate/rng.gamma(am_shape);a=am+np.sqrt(sm2)*(Lmchol@rng.normal(size=P))
        if it>=warmup:
            j=it-warmup
            for key,val in [('B',B),('D',D),('sd',np.sqrt(sig)),('a',a),('sm',np.sqrt(sm2)),('c',c),('loading',loading),('b',b)]:storage[key][j]=val
        if (it+1)%1000==0:print(f'prior={prior_scale}, independent={independent}, chain={chain+1}, iteration={it+1}/{warmup+draws}',flush=True)
    return storage

def sample(data,draws=2000,warmup=1500,prior_scale=1.,label='joint',independent=False):
    chains=[run_chain(data,ch,draws,warmup,prior_scale,independent) for ch in range(4)]
    result={k:np.stack([c[k] for c in chains]) for k in chains[0]}
    # b is saved separately so population parameters are convenient to inspect.
    np.savez_compressed(OUT/f'{label}_population_draws.npz',**{k:v for k,v in result.items() if k!='b'})
    np.savez_compressed(OUT/f'{label}_participant_draws.npz',b=result['b'])
    return result

def summary_row(name,values,truth=np.nan,**extra):
    a=np.asarray(values).ravel();lo,med,hi=np.quantile(a,[.025,.5,.975])
    return dict(parameter=name,mean=float(a.mean()),median=float(med),lower=float(lo),upper=float(hi),
                probability_positive=float((a>0).mean()),truth=truth,**extra)

def scalar_draws(s):
    variables={}
    for j,f in enumerate(FEATURES+['module']):
        for k,dom in enumerate(DOMAINS):variables[f'{dom}:{f}']=s['B'][:,:,j,k]
    for j,f in enumerate(FEATURES):variables['module:'+f]=s['a'][:,:,j]
    for j,f in enumerate(FEATURES[:12]+['module']):variables['complication:'+f]=s['c'][:,:,j]
    variables['complication:vascular_loading']=s['loading'];variables['module:residual_sd']=s['sm']
    for k,dom in enumerate(DOMAINS):
        variables[f'{dom}:residual_sd']=s['sd'][:,:,k]
        variables[f'{dom}:participant_sd']=np.sqrt(s['D'][:,:,k,k])
        for l in range(k):variables[f'correlation:{DOMAINS[l]}:{dom}']=s['D'][:,:,k,l]/np.sqrt(s['D'][:,:,k,k]*s['D'][:,:,l,l])
    return variables

def diagnostics(s,label):
    from arviz_stats import rhat, ess, mcse
    variables=scalar_draws(s)
    # All participant effects are monitored, not only selected population parameters.
    for i in range(s['b'].shape[2]):
        for k,dom in enumerate(DOMAINS):variables[f'participant_{i+1}:{dom}']=s['b'][:,:,i,k]
    names=list(variables);arr=np.stack(list(variables.values()),axis=0).astype(float)
    rh=np.asarray(rhat(arr,chain_axis=1,draw_axis=2,method='rank'))
    bulk=np.asarray(ess(arr,chain_axis=1,draw_axis=2,method='bulk'))
    tail=np.asarray(ess(arr,chain_axis=1,draw_axis=2,method='tail',prob=.05))
    mc=np.asarray(mcse(arr,chain_axis=1,draw_axis=2,method='mean'))
    res=pd.DataFrame({'parameter':names,'rhat':rh,'ess_bulk':bulk,'ess_tail':tail,'mcse_mean':mc})
    res.to_csv(OUT/f'{label}_diagnostics.csv',index=False)
    diag={'max_rhat':float(rh.max()),'min_bulk_ess':float(bulk.min()),'min_tail_ess':float(tail.min()),
          'parameters_checked':len(names),'rhat_over_1_01':int((rh>=1.01).sum()),
          'bulk_below_400':int((bulk<400).sum()),'tail_below_400':int((tail<400).sum()),
          'sampler':'Blocked Gibbs with exact conditionals; HMC divergences/treedepth/BFMI not applicable'}
    (OUT/f'{label}_diagnostic_summary.json').write_text(json.dumps(diag,indent=2)+'\n')
    print(label,diag,flush=True)
    return diag

def effects(s,data,truth,label='joint'):
    flat={k:v.reshape((-1,)+v.shape[2:]) for k,v in s.items() if k!='b'}
    B,a,c,D,sd,sm,loading=[flat[k] for k in ['B','a','c','D','sd','sm','loading']]
    gc=data['X'][:,0,2];gbar=gc.mean();A=a[:,1]+a[:,3]*gbar
    rows=[];rho_rows=[]
    for k,dom in enumerate(DOMAINS):
        for g in ['population',0,1,2]:
            gg=gbar if g=='population' else g-.6
            med_slope=a[:,1]+a[:,3]*gg
            nde=B[:,1,k]+B[:,3,k]*gg;nie=med_slope*B[:,-1,k]
            at=np.array(truth['a']);bt=np.array(truth['B'])
            for typ,val,tr in [('direct',nde,bt[1,k]+bt[3,k]*gg),('indirect',nie,(at[1]+at[3]*gg)*bt[-1,k]),('total',nde+nie,bt[1,k]+bt[3,k]*gg+(at[1]+at[3]*gg)*bt[-1,k])]:
                rows.append(summary_row(f'{dom}:{typ}',val,tr,outcome=dom,effect=typ,genotype=str(g),scale='standardised outcome units'))
        for rho in np.arange(-.6,.61,.1):
            # Exact linear-normal residual-correlation reparameterisation.
            kappa=rho*sd[:,k]/(sm*np.sqrt(1-rho*rho))
            nie=A*(B[:,-1,k]-kappa)
            rho_rows.append(summary_row(dom,nie,rho=float(round(rho,2)),outcome=dom,scale='standardised outcome units'))
    # Population-standardised natural effects for antenatal +1 SD exposure shift.
    # Average over empirical baseline covariates, integrate new-subject random effects
    # and mediator residuals analytically on the probit probability scale.
    Xm=data['X'][:,:3].mean(axis=1);Z=data['Z'].copy();Z[:,-1]=0
    natural={key:[] for key in ['direct','indirect','total']};truth_event={}
    rho_event={round(float(rho),2):[] for rho in np.arange(-.6,.61,.1)}
    # Use all population draws, with batches to bound memory.
    for start in range(0,len(a),200):
        sl=slice(start,start+200);aa=a[sl];cc=c[sl];muM=aa@Xm.T
        dM=aa[:,1,None]+aa[:,3,None]*gc[None,:]
        direct=cc[:,1,None]+cc[:,3,None]*gc[None,:]
        base=cc@Z.T+cc[:,-1,None]*muM
        denom=np.sqrt(1+loading[sl]**2*D[sl,0,0]+cc[:,-1]**2*sm[sl]**2/3)[:,None]
        r00=ndtr(base/denom).mean(axis=1)
        r10=ndtr((base+direct)/denom).mean(axis=1)
        r11=ndtr((base+direct+cc[:,-1,None]*dM)/denom).mean(axis=1)
        for key,val in [('direct',r10-r00),('indirect',r11-r10),('total',r11-r00)]:natural[key].extend(val)
        for rho in rho_event:
            kap=rho/(sm[sl]/np.sqrt(3)*np.sqrt(1-rho*rho))
            # Total is invariant; changing residual correlation reallocates pathways.
            r10rho=ndtr((base+direct+kap[:,None]*dM)/denom).mean(axis=1)
            rho_event[rho].extend(r11-r10rho)
    at=np.array(truth['a']);ct=np.array(truth['clinical_c']);dt=np.array(truth['D'])
    muM=Xm@at;dM=at[1]+at[3]*gc;direct=ct[1]+ct[3]*gc
    base=Z@ct+ct[-1]*muM;denom=np.sqrt(1+truth['loading']**2*dt[0,0]+ct[-1]**2*truth['mediator_sd']**2/3)
    r00=ndtr(base/denom).mean();r10=ndtr((base+direct)/denom).mean();r11=ndtr((base+direct+ct[-1]*dM)/denom).mean()
    truth_event=dict(direct=r10-r00,indirect=r11-r10,total=r11-r00)
    for typ,val in natural.items():rows.append(summary_row(f'complication:{typ}',val,truth_event[typ],outcome='complication',effect=typ,genotype='population',scale='risk difference'))
    for rho,val in rho_event.items():rho_rows.append(summary_row('complication',val,rho=rho,outcome='complication',scale='risk difference'))
    pd.DataFrame(rows).to_csv(OUT/f'{label}_causal_effects.csv',index=False)
    pd.DataFrame(rho_rows).to_csv(OUT/f'{label}_mediation_sensitivity.csv',index=False)
    return pd.DataFrame(rows),pd.DataFrame(rho_rows)

def posterior_checks(s,data):
    rng=np.random.default_rng(SEED+99);rows=[];prior_rows=[]
    X,W,Z,y,med,event=[data[k] for k in ['X','W','Z','y','mediator','event']]
    n,V,K=y.shape
    def stats(y,e,m):
        ans={'complication_rate':e.mean(),'module_mean':m.mean(),'module_sd':m.std()}
        for k,dom in enumerate(DOMAINS):
            ans[f'{dom}_mean']=y[:,:,k].mean();ans[f'{dom}_sd']=y[:,:,k].std()
            ans[f'{dom}_antenatal_postpartum_change']=(y[:,3,k]-y[:,:3,k].mean(axis=1)).mean()
        for j,k in [(0,1),(0,2),(1,2)]:ans[f'correlation_{DOMAINS[j]}_{DOMAINS[k]}']=np.corrcoef(y[:,:,j].mean(axis=1),y[:,:,k].mean(axis=1))[0,1]
        return ans
    observed=stats(y,event,med);reps={k:[] for k in observed};priors={k:[] for k in observed}
    for t in range(400):
        ch=t%4;j=int(rng.integers(s['B'].shape[1]))
        # New-subject replication at observed covariates, preserving joint structure.
        b=rng.multivariate_normal(np.zeros(3),s['D'][ch,j],n)
        mm=np.einsum('nvp,p->nv',X,s['a'][ch,j])+rng.normal(0,s['sm'][ch,j],(n,V))
        ww=np.concatenate([X,mm[:,:,None]],axis=2)
        yy=np.einsum('nvp,pk->nvk',ww,s['B'][ch,j])+b[:,None,:]+rng.normal(size=y.shape)*s['sd'][ch,j]
        zz=Z.copy();zz[:,-1]=mm[:,:3].mean(axis=1)
        ee=rng.binomial(1,ndtr(zz@s['c'][ch,j]+s['loading'][ch,j]*b[:,0]))
        for key,val in stats(yy,ee,mm).items():reps[key].append(val)
        # Actual prior predictive simulation under the base model priors.
        ss=2/rng.gamma(3,size=3);dd=invwishart.rvs(df=5,scale=np.eye(3)*.5,random_state=rng)
        bp=rng.multivariate_normal(np.zeros(3),dd,n);beta=rng.normal(size=(P+1,3));beta[0]*=2
        sma=np.sqrt(2/rng.gamma(3));aa=rng.normal(0,sma,P)
        mp=np.einsum('nvp,p->nv',X,aa)+rng.normal(0,sma,(n,V))
        wp=np.concatenate([X,mp[:,:,None]],axis=2)
        yp=np.einsum('nvp,pk->nvk',wp,beta)+bp[:,None,:]+rng.normal(size=y.shape)*np.sqrt(ss)
        cp=rng.normal(size=Z.shape[1]);cp[0]*=2;zp=Z.copy();zp[:,-1]=mp[:,:3].mean(axis=1)
        ep=rng.binomial(1,ndtr(zp@cp+rng.normal()*bp[:,0]))
        for key,val in stats(yp,ep,mp).items():priors[key].append(val)
    for key,obs in observed.items():
        rr=np.array(reps[key]);pp=np.array(priors[key])
        rows.append(dict(statistic=key,observed=obs,lower=np.quantile(rr,.025),median=np.median(rr),upper=np.quantile(rr,.975),predictive_tail_probability=np.mean(rr>=obs)))
        prior_rows.append(dict(statistic=key,lower=np.quantile(pp,.025),median=np.median(pp),upper=np.quantile(pp,.975)))
    pd.DataFrame(rows).to_csv(OUT/'posterior_predictive_checks.csv',index=False)
    pd.DataFrame(prior_rows).to_csv(OUT/'prior_predictive_checks.csv',index=False)
    return pd.DataFrame(rows)

def transcriptomic_scan(data):
    """Exploratory post-joint scan; mediator module was fixed before this scan."""
    rng=np.random.default_rng(SEED+88);X=data['X'];n,V,p=X.shape;m=500
    effects=np.r_[np.repeat(.3,25),np.zeros(m-25)]
    subject=rng.normal(0,.6,(n,1,m));cell=rng.normal(size=(n,V));batch=rng.integers(0,2,(n,V))
    expr=subject+X[:,:,1,None]*effects+rng.normal(size=(n,V,m))+.15*cell[:,:,None]+.1*batch[:,:,None]
    # Demeaning removes all invariant covariates; adjust visits and seasonal harmonics.
    xx=np.concatenate([X[:,:,[1,3,10,11,12,13,14]],cell[:,:,None],batch[:,:,None]],axis=2)
    xx-=xx.mean(axis=1,keepdims=True);yy=expr-expr.mean(axis=1,keepdims=True)
    xf=xx.reshape(n*V,-1);yf=yy.reshape(n*V,m);inv=np.linalg.inv(xf.T@xf);b=inv@xf.T@yf
    res=(yf-xf@b).reshape(n,V,m);score=np.einsum('nvk,nvm->nkm',xx,res);influence=np.einsum('kl,nlm->nkm',inv,score)
    se=np.sqrt((influence**2).sum(axis=0)*n/(n-1)*(n*V-1)/(n*V-n-xf.shape[1]))
    pv=2*student_t.sf(np.abs(b[0]/se[0]),df=n-1);order=np.argsort(pv);q=np.empty(m)
    q[order]=np.minimum.accumulate((pv[order]*m/np.arange(1,m+1))[::-1])[::-1];q=np.minimum(q,1)
    r=pd.DataFrame({'feature':[f'EXT_GENE_{i+1:03}' for i in range(m)],'truth':effects,'beta':b[0],'se':se[0],'p':pv,'q':q})
    r.to_csv(OUT/'transcriptomic_scan.csv',index=False)
    # Separate synthetic matrix, not the source of the prespecified mediator module.
    np.savez_compressed(OUT/'synthetic_expression_scan.npz',expression=expr.astype(np.float32),cell=cell,batch=batch)
    return {'true_positives':int(((q<.05)&(effects!=0)).sum()),'false_positives':int(((q<.05)&(effects==0)).sum()),'injected':25,'tested':500}

def figures(s,data,eff,sens,ppc):
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
    d=pd.read_csv(OUT/'synthetic_longitudinal.csv')
    fig,ax=plt.subplots(1,2,figsize=(10,3.5))
    for site,g in d.groupby('site'):
        v=g.groupby('visit').temperature_C.mean();ax[0].plot(v.index,v,label=f'Synthetic site {site+1}',marker='o')
    ax[0].set(xlabel='Sample visit (4 = postpartum)',ylabel='Mean preceding 7-day Tmax (C)',xticks=[1,2,3,4]);ax[0].legend(fontsize=8)
    for k,dom in enumerate(DOMAINS):
        yy=data['y'][:,:,k];ax[1].errorbar(np.arange(1,5),yy.mean(axis=0),yerr=1.96*yy.std(axis=0)/np.sqrt(len(yy)),marker='o',label=dom)
    ax[1].set(xlabel='Outcome visit, one week after sample',ylabel='Synthetic standardised outcome',xticks=[1,2,3,4]);ax[1].legend(fontsize=8)
    fig.suptitle('SIMULATED spatiotemporal exposure and joint outcome trajectories');fig.tight_layout();fig.savefig(OUT/'05_joint_trajectories.png',dpi=180);plt.close(fig)
    fig,ax=plt.subplots(2,2,figsize=(10,6))
    for k,dom in enumerate(DOMAINS+['complication']):
        z=eff[(eff.outcome==dom)&(eff.genotype=='population')];a=ax.ravel()[k]
        for j,(_,r) in enumerate(z.iterrows()):
            a.errorbar(r['mean'],j,xerr=[[r['mean']-r.lower],[r.upper-r['mean']]],fmt='o',capsize=3);a.scatter(r.truth,j,c='black',marker='x')
        a.set(yticks=range(3),yticklabels=z.effect,title=dom,xlabel='Risk difference' if dom=='complication' else 'Standardised units');a.axvline(0,c='grey',ls=':')
    fig.suptitle('Bayesian effects for +1 SD heat; 95% credible intervals; x = truth');fig.tight_layout();fig.savefig(OUT/'06_joint_effects.png',dpi=180);plt.close(fig)
    fig,ax=plt.subplots(2,2,figsize=(10,6))
    for k,dom in enumerate(DOMAINS+['complication']):
        z=sens[sens.outcome==dom];a=ax.ravel()[k];a.plot(z.rho,z['mean']);a.fill_between(z.rho,z.lower,z.upper,alpha=.2)
        a.axhline(0,c='gray',ls=':');a.axvline(0,c='gray',ls=':');a.set(title=dom,xlabel='Assumed residual correlation rho',ylabel='Indirect effect')
    fig.suptitle('Mediation sensitivity: one outcome at a time; rho is not estimated');fig.tight_layout();fig.savefig(OUT/'07_mediation_sensitivity.png',dpi=180);plt.close(fig)
    fig,ax=plt.subplots(3,1,figsize=(10,6))
    vars=[s['a'][:,:,1],s['B'][:,:,-1,1],s['loading']]
    for a,v,name in zip(ax,vars,['Heat to module','Module to depression','Vascular loading to complication']):
        for ch in range(4):a.plot(v[ch],alpha=.7,lw=.5,label=f'Chain {ch+1}')
        a.set(ylabel=name,xlabel='Retained draw')
    ax[0].legend(ncol=4,fontsize=7);fig.suptitle('Post-warmup Gibbs traces');fig.tight_layout();fig.savefig(OUT/'08_traceplots.png',dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,5));chosen=ppc[~ppc.statistic.str.endswith('_sd')]
    for j,(_,r) in enumerate(chosen.iterrows()):ax.plot([r.lower,r.upper],[j,j],c='#176b87',lw=3);ax.scatter(r.observed,j,c='#c26927',marker='x',zorder=3)
    ax.set(yticks=range(len(chosen)),yticklabels=chosen.statistic,xlabel='Statistic',title='New-subject posterior predictive checks: intervals; x = observed')
    fig.tight_layout();fig.savefig(OUT/'09_predictive_checks.png',dpi=180);plt.close(fig)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--draws',type=int,default=2000);ap.add_argument('--warmup',type=int,default=1500)
    ap.add_argument('--n',type=int,default=600);ap.add_argument('--reuse',action='store_true');ap.add_argument('--skip-sensitivity-fit',action='store_true');args=ap.parse_args()
    start=time.time();data,truth=create_data(args.n)
    def get(label,prior=1.,independent=False):
        if args.reuse and (OUT/f'{label}_population_draws.npz').exists():
            s=dict(np.load(OUT/f'{label}_population_draws.npz'));s['b']=np.load(OUT/f'{label}_participant_draws.npz')['b'];return s
        return sample(data,args.draws,args.warmup,prior,label,independent)
    joint=get('joint');diag=diagnostics(joint,'joint');eff,sens=effects(joint,data,truth)
    pop=scalar_draws(joint);rows=[summary_row(k,v) for k,v in pop.items()]
    pd.DataFrame(rows).to_csv(OUT/'joint_parameter_summary.csv',index=False)
    ppc=posterior_checks(joint,data);scan=transcriptomic_scan(data);figures(joint,data,eff,sens,ppc)
    comparisons=[]
    if not args.skip_sensitivity_fit:
        # Same data, wider regression priors. Keeps covariance prior fixed.
        prior=get('wide_prior',2.);diagnostics(prior,'wide_prior');ep,_=effects(prior,data,truth,'wide_prior')
        for _,r in ep[ep.genotype=='population'].iterrows():
            rr=eff[(eff.parameter==r.parameter)&(eff.genotype=='population')].iloc[0]
            comparisons.append(dict(parameter=r.parameter,base_mean=rr['mean'],wide_mean=r['mean'],base_lower=rr.lower,base_upper=rr.upper,wide_lower=r.lower,wide_upper=r.upper))
        pd.DataFrame(comparisons).to_csv(OUT/'prior_sensitivity_comparison.csv',index=False)
        del prior
    run={'seed':SEED,'n':args.n,'visits':4,'chains':4,'draws_per_chain':joint['B'].shape[1],
         'warmup_per_chain':args.warmup,'clinical_events':truth['complications'],'diagnostics':diag,
         'scan':scan,'seconds':time.time()-start,'python':platform.python_version(),
         'status':'SIMULATED - no empirical or clinical conclusions'}
    (OUT/'run_summary.json').write_text(json.dumps(run,indent=2)+'\n');print(json.dumps(run,indent=2),flush=True)

if __name__=='__main__':main()
