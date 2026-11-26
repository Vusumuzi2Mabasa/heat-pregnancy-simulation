"""Build the consolidated research manuscript from unchanged saved results."""
from pathlib import Path
import argparse,os,shutil,subprocess,sys
ROOT=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--quarto',default=shutil.which('quarto'));p.add_argument('--pdf',action='store_true');a=p.parse_args()
if not a.quarto:sys.exit('Quarto is required; add it to PATH or supply --quarto /path/to/quarto')
os.environ['R_ENVIRON_USER']='NUL' if os.name=='nt' else '/dev/null';os.environ['R_PROFILE_USER']=os.environ['R_ENVIRON_USER']
def run(*cmd):subprocess.run([str(c) for c in cmd],cwd=ROOT,check=True)
run(sys.executable,ROOT/'write_manuscript.py');run(sys.executable,ROOT/'validate_manuscript.py')
subprocess.run([a.quarto,'render','Heat_Pregnancy_Manuscript.qmd','--to','docx','--output','heat_pregnancy_simulation_report.docx'],cwd=ROOT/'report',check=True)
target=ROOT/'output/docx/heat_pregnancy_simulation_report.docx';target.parent.mkdir(parents=True,exist_ok=True)
shutil.move(ROOT/'report/heat_pregnancy_simulation_report.docx',target)
run(sys.executable,ROOT/'style_manuscript.py',target)
if a.pdf:
 pdf=ROOT/'output/pdf/heat_pregnancy_simulation_report.pdf'
 if os.name=='nt':run(shutil.which('pwsh') or 'powershell','-NoProfile','-File',ROOT/'export_manuscript.ps1','-DocumentPath',target,'-PdfPath',pdf)
 elif shutil.which('soffice'):
  run(shutil.which('soffice'),'--headless','--convert-to','pdf','--outdir',pdf.parent,target)
 else:sys.exit('PDF export needs Microsoft Word on Windows or LibreOffice; Word manuscript was saved')
print('Saved manuscript:',target)
