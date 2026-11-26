"""Apply manuscript typography after Quarto rendering using the bundled runtime."""
import sys
from docx import Document
from docx.shared import Inches,Pt,RGBColor
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
d=Document(sys.argv[1])
for sec in d.sections:
 sec.page_width=Inches(8.5);sec.page_height=Inches(11);sec.top_margin=sec.bottom_margin=Inches(.8);sec.left_margin=sec.right_margin=Inches(.85)
 sec.header_distance=sec.footer_distance=Inches(.35)
 p=sec.footer.paragraphs[0];p.clear();p.alignment=2;r=p.add_run('Vusumuzi Mabasa  |  Simulation study  |  ');r.font.size=Pt(8)
 field=OxmlElement('w:fldSimple');field.set(qn('w:instr'),'PAGE');p._p.append(field)
for s in d.styles:
 if s.type==1:
  s.font.name='Calibri';s.font.color.rgb=RGBColor(0,0,0)
  if s.name in ['Normal','Body Text']:s.font.size=Pt(11);s.paragraph_format.space_after=Pt(6)
d.styles['Title'].font.size=Pt(23)
for p in d.paragraphs:
 if p.style.name=='Title':
  pp=p._p.get_or_add_pPr()
  for b in list(pp.findall(qn('w:pBdr'))):pp.remove(b)
 if p.style.name.startswith('Heading'):p.paragraph_format.keep_with_next=True
for tab in d.tables:
 borders=OxmlElement('w:tblBorders')
 for n in ['top','left','bottom','right','insideH','insideV']:
  b=OxmlElement('w:'+n);b.set(qn('w:val'),'single');b.set(qn('w:sz'),'4');b.set(qn('w:color'),'D9D9D9');borders.append(b)
 tab._tbl.tblPr.append(borders)
 for i,row in enumerate(tab.rows):
  if i==0:row._tr.get_or_add_trPr().append(OxmlElement('w:tblHeader'))
  for cell in row.cells:
   shade=OxmlElement('w:shd');shade.set(qn('w:fill'),'E6EEF2' if i==0 else 'FFFFFF');cell._tc.get_or_add_tcPr().append(shade)
   for p in cell.paragraphs:
    p.paragraph_format.space_after=Pt(3);p.paragraph_format.space_before=Pt(3)
    for r in p.runs:r.font.size=Pt(9)
d.core_properties.author='Vusumuzi Mabasa';d.core_properties.subject='Public weather and simulated pregnancy outcomes';d.save(sys.argv[1])
