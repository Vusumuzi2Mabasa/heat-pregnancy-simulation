"""Compatibility entry point for the current R/Quarto association manuscript.
The earlier builder is preserved in archive/bayesian-study-v2.
"""
from pathlib import Path
import shutil
import subprocess
ROOT = Path(__file__).resolve().parent
rscript = shutil.which('Rscript') or 'C:/Program Files/R/R-4.5.1/bin/Rscript.exe'
subprocess.run([rscript, '--vanilla', str(ROOT/'run_epidemiology.R'), '--report'], cwd=ROOT, check=True)
if shutil.which('pwsh'):
    subprocess.run(['pwsh', '-NoProfile', '-File', str(ROOT/'export_manuscript.ps1'),
        '-DocumentPath', str(ROOT/'output/docx/heat_pregnancy_simulation_report.docx'),
        '-PdfPath', str(ROOT/'output/pdf/heat_pregnancy_simulation_report.pdf')],check=True)
