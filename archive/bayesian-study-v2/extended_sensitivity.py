"""Assay-error stress test: naive model refit after adding module measurement error."""
import json, argparse
import numpy as np
import pandas as pd
from extend_simulation import OUT, SEED, sample, diagnostics, effects

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--reuse',action='store_true');args=ap.parse_args()
    data=dict(np.load(OUT/'analysis_inputs.npz'))
    truth=json.loads((OUT/'generating_truth.json').read_text())
    rng=np.random.default_rng(SEED+777)
    sd=.5
    noisy=data['mediator']+rng.normal(0,sd,data['mediator'].shape)
    data['mediator']=noisy;data['W']=data['W'].copy();data['W'][:,:,-1]=noisy
    data['Z']=data['Z'].copy();data['Z'][:,-1]=noisy[:,:3].mean(axis=1)
    np.savez_compressed(OUT/'noisy_module_inputs.npz',**data)
    if args.reuse and (OUT/'noisy_module_population_draws.npz').exists():
        s=dict(np.load(OUT/'noisy_module_population_draws.npz'));s['b']=np.load(OUT/'noisy_module_participant_draws.npz')['b']
    else:s=sample(data,2000,1500,label='noisy_module')
    diagnostics(s,'noisy_module');eff,_=effects(s,data,truth,'noisy_module')
    base=pd.read_csv(OUT/'joint_causal_effects.csv')
    a=base[base.genotype=='population'].merge(eff[eff.genotype=='population'],on=['outcome','effect'],suffixes=('_base','_noisy'))
    a[['outcome','effect','mean_base','lower_base','upper_base','mean_noisy','lower_noisy','upper_noisy','truth_base']].to_csv(OUT/'assay_error_comparison.csv',index=False)
    (OUT/'assay_error_design.json').write_text(json.dumps({'added_noise_sd':sd,
      'conditional_reliability':truth['mediator_sd']**2/(truth['mediator_sd']**2+sd**2),
      'interpretation':'Sensitivity stress test, not an error-corrected model; same outcomes and weather, corrupted module only.'},indent=2)+'\n')

if __name__=='__main__':main()
