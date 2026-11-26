"""Verify the current manuscript against the preserved scientific release."""
from pathlib import Path
import hashlib,json,re
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parent;E=ROOT/'results/extended'
manifest=json.loads((ROOT/'provenance/original_results_sha256.json').read_text())
canonical_matches=[]
for name,sha in manifest['files'].items():
 h=hashlib.sha256()
 with (ROOT/name).open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 if h.hexdigest()!=sha:
  expected=manifest.get('canonical_text_sha256',{}).get(name)
  assert expected and hashlib.sha256((ROOT/name).read_bytes().replace(b'\r\n',b'\n')).hexdigest()==expected,'Original result changed: '+name
  canonical_matches.append(name)
s=json.loads((E/'run_summary.json').read_text());clinical=pd.read_csv(E/'synthetic_clinical.csv');long=pd.read_csv(E/'synthetic_longitudinal.csv')
assert s['n']==600 and s['visits']==4 and s['clinical_events']==86
assert len(long)==2400 and len(clinical)==600 and clinical.complication.sum()==86
assert (long.outcome_day-long.sample_day).eq(7).all()
assert (long[long.visit<=3].merge(clinical,on='id').outcome_day<long[long.visit<=3].merge(clinical,on='id').ascertainment_day).all()
effects=pd.read_csv(E/'joint_causal_effects.csv');population=effects[effects.genotype=='population']
for outcome,g in population.groupby('outcome'):
 t=g.set_index('effect')['mean'];assert np.isclose(t['total'],t['direct']+t['indirect'],atol=1e-10)
for k in ['joint','wide_prior','noisy_module']:
 d=pd.read_csv(E/f'{k}_diagnostics.csv');assert len(d)==1887 and (d.rhat<1.01).all() and (d.ess_bulk>400).all() and (d.ess_tail>400).all()
pp=pd.read_csv(E/'posterior_predictive_checks.csv');assert len(pp)==15 and ((pp.observed>=pp.lower)&(pp.observed<=pp.upper)).all()
scan=pd.read_csv(E/'transcriptomic_scan.csv');selected=scan.q<.05;assert ((scan.truth>0)&selected).sum()==22 and ((scan.truth==0)&selected).sum()==1
readme=(ROOT/'README.md').read_text(encoding='utf-8');source=(ROOT/'report/Heat_Pregnancy_Manuscript.qmd').read_text(encoding='utf-8')
for name in ['Abstract','Introduction','Methods','Results','Discussion','Conclusion','References']:assert f'## {name}' in readme and f'## {name}' in source
images=re.findall(r'!\[[^]]*\]\(([^)]+)\)',readme);assert len(images)==5 and all((ROOT/p).is_file() for p in images)
assert not re.search(r'(DESIGN_TABLE|EFFECT_TABLE|ASSETS/|TRAJECTORY_FIGURE)',readme)
result=dict(status='PASS',original_result_files_unchanged=len(manifest['files']),participants=600,visits=2400,complications=86,figures_in_readme=5,manuscript_structure=True,all_three_diagnostic_sets_pass=True,posterior_predictive_summaries_checked=15,spatial_analysis_added=False,scope='Saved-result manuscript revision; no posterior resampling')
result['text_files_matching_after_line_ending_normalisation']=canonical_matches
(ROOT/'provenance/manuscript_validation.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
