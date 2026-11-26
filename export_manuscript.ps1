param([Parameter(Mandatory=$true)][string]$DocumentPath,[Parameter(Mandatory=$true)][string]$PdfPath)
$ErrorActionPreference='Stop'
$wordApp=New-Object -ComObject Word.Application
$wordApp.Visible=$false
$wordApp.DisplayAlerts=0
try {
 $manuscriptDoc=$wordApp.Documents.Open($DocumentPath,$false,$true)
 try { $manuscriptDoc.ExportAsFixedFormat($PdfPath,17) } finally { $manuscriptDoc.Close(0) }
} finally { $wordApp.Quit() }
