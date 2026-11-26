"""Synthetic design exercise. No real people, genotypes, or weather observations."""
from pathlib import Path
import json, os, platform, sys
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.stats import norm, t
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'results'
OUT.mkdir(exist_ok=True)
SEED = 20261023
R = np.random.default_rng(SEED)
CFG = dict(seed=SEED, recruited=1500, retention=.90, molecular_n=300,
           genes=500, true_heat_genes=25, repetitions=300,
           warming_C_per_year=.04, heat_OR=1.25, interaction_OR=1.15,
           allele_frequency_reference=.30, expression_effect=.30,
           statement='ALL DATA AND EFFECT SIZES ARE SYNTHETIC ASSUMPTIONS')

def fit_logit(X, y):
    b = np.zeros(X.shape[1])
    for iteration in range(60):
        p = expit(X @ b)
        w = np.maximum(p * (1-p), 1e-8)
        h = (X.T * w) @ X
        step = np.linalg.solve(h, X.T @ (y-p))
        b += step
        if np.max(np.abs(step)) < 1e-8:
            p = expit(X @ b)
            h = (X.T * np.maximum(p*(1-p), 1e-8)) @ X
            se = np.sqrt(np.diag(np.linalg.inv(h)))
            return b, se, 2*norm.sf(np.abs(b/se))
    raise RuntimeError('Logistic model did not converge')

def weather():
    dates = pd.date_range('2010-01-01', '2030-12-31')
    years = (dates - dates[0]).days.to_numpy()/365.25
    season = 5*np.cos(2*np.pi*(dates.dayofyear.to_numpy()-20)/365.25)
    common = np.zeros(len(dates))
    innovations = R.normal(0, 1.7, len(dates))
    for j in range(1, len(dates)):
        common[j] = .65*common[j-1]+innovations[j]
    records=[]
    for s, base in enumerate([24., 27., 29.]):
        vals = base + season + .04*years + common + R.normal(0,.7,len(dates))
        for k in range(6):
            offset = (k-2.5)*.35
            records.append(pd.DataFrame(dict(date=dates, site=s, neighbourhood=k,
                tmax_C=vals+offset, year=dates.year, month=dates.month)))
    df=pd.concat(records, ignore_index=True)
    df.to_csv(OUT/'synthetic_daily_weather.csv.gz',index=False,compression={'method':'gzip','mtime':0})
    annual=df.groupby(['site','year']).tmax_C.mean().reset_index()
    annual.to_csv(OUT/'annual_weather.csv',index=False)
    trend=[]
    for s, d in annual.groupby('site'):
        X=np.column_stack([np.ones(len(d)), d.year.to_numpy()-2010])
        y=d.tmax_C.to_numpy(); inv=np.linalg.inv(X.T@X); b=inv@X.T@y
        e=y-X@b; hat=np.sum((X@inv)*X,axis=1)
        cov=inv@(X.T*(e/(1-hat))**2)@X@inv
        se=np.sqrt(cov[1,1]); q=t.ppf(.975,len(y)-2)
        trend.append(dict(site=s, slope=b[1], lower=b[1]-q*se, upper=b[1]+q*se))
    pd.DataFrame(trend).to_csv(OUT/'trend_estimates.csv',index=False)
    fig,ax=plt.subplots(1,2,figsize=(10,3.7))
    for s,d in annual.groupby('site'): ax[0].plot(d.year,d.tmax_C,label=f'Site {s+1}')
    ax[0].set(xlabel='Synthetic calendar year',ylabel='Annual mean daily maximum (C)',title='Assumed warming: 0.04 C/year')
    ax[0].legend(frameon=False)
    m=df[df.year.between(2027,2029)].groupby(['site','month']).tmax_C.mean().unstack()
    im=ax[1].imshow(m,aspect='auto',cmap='YlOrRd'); ax[1].set(xticks=range(12),xticklabels=range(1,13),yticks=range(3),yticklabels=['Site 1','Site 2','Site 3'],xlabel='Month',title='Synthetic spatial / seasonal contrast')
    fig.colorbar(im,ax=ax[1],label='Tmax (C)'); fig.suptitle('SIMULATED WEATHER - not observed South African trends',fontsize=11)
    fig.tight_layout();fig.savefig(OUT/'01_spatiotemporal.png',dpi=180);plt.close(fig)
    return df

def pool_from_weather(w):
    n=18000
    site=R.integers(0,3,n); hood=R.integers(0,6,n)
    start=R.integers(0,1095,n)
    conception=pd.Timestamp('2027-01-01')+pd.to_timedelta(start,unit='D')
    lookup={key:g.sort_values('date').tmax_C.to_numpy() for key,g in w.groupby(['site','neighbourhood'])}
    day=(conception-pd.Timestamp('2010-01-01')).days.to_numpy()
    exposure=np.empty(n)
    for i in range(n): exposure[i]=lookup[(site[i],hood[i])][day[i]+98:day[i]+140].mean()
    # Standardised relative to the synthetic source population, not each retained sample.
    z=(exposure-exposure.mean())/exposure.std()
    pool=pd.DataFrame(dict(site=site,neighbourhood=hood,conception=conception,
                           day_index=day,heat_C=exposure,heat_z=z,
                           season_sin=np.sin(2*np.pi*start/365.25),season_cos=np.cos(2*np.pi*start/365.25)))
    return pool,lookup

def cohort(n, rng, pool, int_or=1.15, noise=.4, heat_or=1.25):
    d=pool.iloc[rng.integers(0,len(pool),n)].reset_index(drop=True).copy()
    d['pc']=rng.normal(.3*(d.site.to_numpy()-1),1,n)
    # Synthetic allele; ancestry proxy is observed exactly in this favourable scenario.
    freq=expit(np.log(.3/.7)+.20*d.pc.to_numpy())
    d['G']=rng.binomial(2,freq); d['Gc']=d.G-.6
    d['age_z']=rng.normal(0,1,n);d['deprivation']=rng.normal(.3*d.site.to_numpy(),1,n)
    h=d.heat_z.to_numpy();g=d.Gc.to_numpy()
    eta=(-2.25+np.log(heat_or)*h+.10*g+np.log(int_or)*h*g+
         .15*d.age_z.to_numpy()+.20*d.deprivation.to_numpy()+.10*d.pc.to_numpy()+
         .08*d.site.to_numpy()+.10*d.season_sin.to_numpy())
    d['risk']=expit(eta);d['preterm']=rng.binomial(1,d.risk)
    d['heat_observed']=h+rng.normal(0,noise,n)
    d['retained']=rng.random(n)<.9
    d['synthetic_id']=[f'SYN-{i+1:05}' for i in range(n)]
    return d

def design(d, heat='heat_observed'):
    h=d[heat].to_numpy(); g=d.Gc.to_numpy()
    return np.column_stack([np.ones(len(d)),h,g,h*g,d.age_z,d.deprivation,d.pc,
                            d.site==1,d.site==2,d.season_sin,d.season_cos])

def main_analysis(d):
    rows=[];keep=d[d.retained]
    for name,col in [('Oracle exposure','heat_z'),('Assigned exposure','heat_observed')]:
        b,se,p=fit_logit(design(keep,col),keep.preterm.to_numpy())
        for j,term,truth in [(1,'Heat',1.25),(3,'Heat x allele',1.15)]:
            rows.append(dict(model=name,term=term,OR=np.exp(b[j]),lower=np.exp(b[j]-1.96*se[j]),upper=np.exp(b[j]+1.96*se[j]),p=p[j],assumed_OR=truth))
    res=pd.DataFrame(rows);res.to_csv(OUT/'clinical_estimates.csv',index=False)
    fig,ax=plt.subplots(figsize=(8,3.5))
    for i,row in res.iterrows():
        ax.errorbar(row.OR,i,xerr=[[row.OR-row.lower],[row.upper-row.OR]],fmt='o',color='#176b87',capsize=4)
        ax.scatter(row.assumed_OR,i,marker='x',color='#c26927',zorder=3)
    ax.set(yticks=range(len(res)),yticklabels=res.model+' | '+res.term,xlabel='Odds ratio with 95% Wald interval',title='SIMULATED clinical estimates; crosses = generating values')
    ax.axvline(1,color='gray',ls=':');fig.tight_layout();fig.savefig(OUT/'02_clinical.png',dpi=180);plt.close(fig)
    return res

def molecular(d,lookup):
    sub=d[d.retained].sample(n=300,random_state=SEED).reset_index(drop=True)
    n=len(sub); v=3; weeks=np.array([14,20,26]); m=500
    exposure=np.empty((n,v))
    for i,row in sub.iterrows():
        vals=lookup[(row.site,row.neighbourhood)]
        for j,wk in enumerate(weeks): exposure[i,j]=vals[int(row.day_index)+wk*7-7:int(row.day_index)+wk*7].mean()
    heat=(exposure-exposure.mean())/exposure.std()
    gest=np.tile((weeks-20)/6,(n,1));cells=R.normal(size=(n,v));batch=R.integers(0,2,(n,v))
    g=sub.Gc.to_numpy()[:,None]
    # Normalised expression surrogates, NOT RNA-seq counts or biological genes.
    effects=np.r_[np.repeat(.30,25),np.zeros(m-25)]
    expr=R.normal(0,.7,(n,1,m))+heat[:,:,None]*effects+R.normal(0,1,(n,v,m))
    expr+=.12*gest[:,:,None]+.18*cells[:,:,None]+.10*batch[:,:,None]
    X=np.stack([heat,heat*g,gest,cells,batch],axis=2)
    X-=X.mean(axis=1,keepdims=True);Y=expr-expr.mean(axis=1,keepdims=True)
    xf=X.reshape(-1,5);yf=Y.reshape(-1,m);inv=np.linalg.inv(xf.T@xf);b=inv@xf.T@yf
    e=(yf-xf@b).reshape(n,v,m)
    scores=np.einsum('nvk,nvm->nkm',X,e)
    influence=np.einsum('kl,nlm->nkm',inv,scores)
    # Cluster CR1 adjustment includes absorbed subject fixed effects.
    correction=n/(n-1)*(n*v-1)/(n*v-n-5)
    se=np.sqrt((influence**2).sum(axis=0)*correction)
    p=2*t.sf(np.abs(b[0]/se[0]),df=n-1)
    order=np.argsort(p);q=np.empty(m);q[order]=np.minimum.accumulate((p[order]*m/np.arange(1,m+1))[::-1])[::-1];q=np.minimum(q,1)
    result=pd.DataFrame(dict(feature=[f'SYN_GENE_{j+1:03}' for j in range(m)],true_beta=effects,beta=b[0],se=se[0],p=p,q=q))
    result.to_csv(OUT/'molecular_estimates.csv',index=False)
    long=pd.DataFrame(dict(synthetic_id=np.repeat(sub.synthetic_id.to_numpy(),v),week=np.tile(weeks,n),heat_z=heat.ravel(),cell_proxy=cells.ravel(),batch=batch.ravel()))
    long.to_csv(OUT/'synthetic_molecular_visits.csv',index=False)
    np.savez_compressed(OUT/'synthetic_expression.npz',expression=expr,features=result.feature.to_numpy(dtype=str))
    fig,ax=plt.subplots(figsize=(8,3.5));ax.scatter(result.beta,-np.log10(np.maximum(result.p,1e-300)),c=np.where(effects>0,'#c26927','#176b87'),alpha=.65,s=14)
    ax.set(xlabel='Within-person heat coefficient',ylabel='-log10(p)',title='SIMULATED expression: orange = 25 injected effects; blue = null')
    fig.tight_layout();fig.savefig(OUT/'03_molecular.png',dpi=180);plt.close(fig)
    return dict(tested=m,true_positive=int(((q<.05)&(effects>0)).sum()),false_positive=int(((q<.05)&(effects==0)).sum()),molecular_n=n,visits=n*v)

def power(pool):
    rows=[];rng=np.random.default_rng(SEED+100)
    scenarios=[('No interaction',1.25,1.),('Modest interaction',1.25,1.15),('Larger interaction',1.25,1.30),('Global heat null',1.,1.)]
    for n in [600,1500,3000]:
        for label,ho,io in scenarios:
            for noise in [0.,.6]:
                vals=[];failures=0
                for rep in range(CFG['repetitions']):
                    d=cohort(n,rng,pool,io,noise,ho);d=d[d.retained]
                    try:
                        b,se,p=fit_logit(design(d),d.preterm.to_numpy())
                        vals.append([p[1]<.05,p[3]<.05,p[3]<5e-8,b[1],b[3],b[3]-1.96*se[3]<=np.log(io)<=b[3]+1.96*se[3],len(d)])
                    except (RuntimeError,np.linalg.LinAlgError):failures+=1
                a=np.asarray(vals);rate=a[:,1].mean()
                rows.append(dict(n=n,scenario=label,heat_OR=ho,interaction_OR=io,noise_sd=noise,
                    successful=len(a),failures=failures,heat_rejection=a[:,0].mean(),interaction_rejection=rate,
                    mcse=np.sqrt(rate*(1-rate)/len(a)),genomewide_rejection=a[:,2].mean(),
                    heat_log_bias=a[:,3].mean()-np.log(ho),interaction_log_bias=a[:,4].mean()-np.log(io),
                    interaction_coverage=a[:,5].mean(),mean_retained=a[:,6].mean()))
                print(f'Completed n={n}, {label}, noise={noise}',flush=True)
    res=pd.DataFrame(rows);res.to_csv(OUT/'design_operating_characteristics.csv',index=False)
    fig,ax=plt.subplots(1,2,figsize=(10,3.6),sharey=True)
    for j,noise in enumerate([0.,.6]):
        for label,_,_ in scenarios[:3]:
            z=res[(res.scenario==label)&(res.noise_sd==noise)]
            ax[j].errorbar(z.n,z.interaction_rejection,yerr=1.96*z.mcse,marker='o',label=label,capsize=3)
        ax[j].axhline(.8,color='gray',ls=':');ax[j].set(title=f'Exposure error SD = {noise}',xlabel='Recruited pregnancies',ylim=(0,1))
    ax[0].set_ylabel('Interaction rejection probability (alpha .05)');ax[1].legend(fontsize=8,frameon=False)
    fig.suptitle('SIMULATION: 300 repetitions per scenario; bars = +/- 1.96 MCSE',fontsize=10)
    fig.tight_layout();fig.savefig(OUT/'04_design.png',dpi=180);plt.close(fig)
    return res

def main():
    (ROOT/'simulation_config.json').write_text(json.dumps(CFG,indent=2)+'\n')
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
    w=weather();pool,lookup=pool_from_weather(w)
    d=cohort(1500,R,pool);d.to_csv(OUT/'synthetic_cohort.csv',index=False)
    main_analysis(d);mol=molecular(d,lookup);power(pool)
    summary=dict(recruited=len(d),retained=int(d.retained.sum()),preterm_events=int(d.loc[d.retained,'preterm'].sum()),
        event_rate=float(d.loc[d.retained,'preterm'].mean()),**mol)
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    import scipy, reportlab
    versions={k:v for k,v in [('python',platform.python_version()),('numpy',np.__version__),('pandas',pd.__version__),('scipy',scipy.__version__),('matplotlib',matplotlib.__version__),('reportlab',reportlab.Version)]}
    (OUT/'software_versions.json').write_text(json.dumps(versions,indent=2)+'\n')
    assert d.synthetic_id.is_unique and d.G.between(0,2).all()
    assert len(w)==len(pd.date_range('2010-01-01','2030-12-31'))*18
    assert (OUT/'design_operating_characteristics.csv').exists()
    assert np.isfinite(pd.read_csv(OUT/'clinical_estimates.csv').OR).all()
    print(json.dumps(summary,indent=2))

if __name__=='__main__': main()
