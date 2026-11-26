"""Record the unchanged legacy scientific release before manuscript revision."""
from pathlib import Path
import hashlib,json,shutil,subprocess
ROOT=Path(__file__).resolve().parent
prov=ROOT/'provenance';prov.mkdir(exist_ok=True)
p=prov/'original_results_sha256.json'
def digest(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
if not p.exists():
 paths=sorted(x for x in (ROOT/'results').rglob('*') if x.is_file())
 p.write_text(json.dumps(dict(original_commit='01e49b7',scope='All original scientific results; manuscript reorganisation only; no new spatial experiment',files={x.relative_to(ROOT).as_posix():digest(x) for x in paths}),indent=2))
 original=ROOT/'_working/original';original.mkdir(parents=True,exist_ok=True)
 shutil.copy2(ROOT/'README.md',original/'README_before_manuscript.md')
 shutil.copy2(ROOT/'output/pdf/heat_pregnancy_simulation_report.pdf',original/'heat_pregnancy_simulation_report_before_manuscript.pdf')
data=json.loads(p.read_text());failed=[name for name,sha in data['files'].items() if digest(ROOT/name)!=sha]
assert not failed,failed
if 'canonical_text_sha256' not in data:
 data['canonical_text_sha256']={name:hashlib.sha256((ROOT/name).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for name in data['files'] if Path(name).suffix in ['.csv','.json']}
 data['text_portability']='Canonical text hashes ignore only CRLF versus LF conversion by Git; binary files require exact bytes'
 p.write_text(json.dumps(data,indent=2))
print('Original result files preserved:',len(data['files']))
