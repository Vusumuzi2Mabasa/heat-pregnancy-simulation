"""Validate the current report and frozen scientific sources."""
from pathlib import Path
import hashlib,json,re
from zipfile import ZipFile
from xml.etree import ElementTree as ET
ROOT=Path(__file__).resolve().parent
NS={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','m':'http://schemas.openxmlformats.org/officeDocument/2006/math'}
manifest=json.loads((ROOT/'provenance/original_results_sha256.json').read_text())
for name,expected in manifest['files'].items():
    b=(ROOT/name).read_bytes()
    if hashlib.sha256(b).hexdigest()!=expected:
        assert hashlib.sha256(b.replace(b'\r\n',b'\n')).hexdigest()==manifest['canonical_text_sha256'].get(name),name
validation=json.loads((ROOT/'results/epidemiology/validation.json').read_text())
assert validation['status']=='passed'
readme=(ROOT/'README.md').read_text(encoding='utf-8')
assert '@@' not in readme
assert len(re.findall(r'!\[.*?\]\(',readme))==9
for target in re.findall(r'\]\(([^)]+)\)',readme):
    if '://' not in target and not target.startswith('#'):
        assert (ROOT/target).exists(),target
with ZipFile(ROOT/'output/docx/heat_pregnancy_simulation_report.docx') as z:
    doc=ET.fromstring(z.read('word/document.xml'));text=' '.join(doc.itertext())
    tables=len(doc.findall('.//w:tbl',NS));equations=len(doc.findall('.//m:oMath',NS));drawings=len(doc.findall('.//w:drawing',NS))
    assert tables==8,(tables,'native Word tables')
    assert drawings==9,(drawings,'figures')
    assert equations>=10,equations
    assert 'Vusumuzi Mabasa' in text and 'synthetic' in text and '@@' not in text
    # flextable explicitly emits invisible paragraph borders inside table cells.
    assert all(edge.get('{'+NS['w']+'}val') in ('none','nil') for border in doc.findall('.//w:pBdr',NS) for edge in border)
    assert 'Vusumuzi Mabasa' in z.read('docProps/core.xml').decode()
result=dict(status='passed',source_files_preserved=len(manifest['files']),native_word_tables=tables,word_figures=drawings,native_equations=equations,README_figures=9,visual_quality_assurance='See epidemiology_quality_assurance.json after render review')
(ROOT/'provenance/epidemiology_manuscript_validation.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result))
