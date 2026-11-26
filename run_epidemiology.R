# Rscript --vanilla run_epidemiology.R [--report]
if (.Platform$OS.type == 'windows') invisible(Sys.setlocale('LC_CTYPE', 'English_United States.utf8'))
args <- commandArgs(trailingOnly=TRUE)
entry <- grep('^--file=', commandArgs(), value=TRUE)
root <- normalizePath(dirname(sub('^--file=', '', entry[1])), winslash='/', mustWork=TRUE)
setwd(root)
required <- c('lme4','lmerTest','gtsummary','flextable','officer','sandwich','car','jsonlite','broom.mixed','ggplot2','emmeans','lpSolve','dplyr','knitr','digest','rmarkdown')
missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if(length(missing)) stop('Install required R packages: ', paste(missing,collapse=', '))
out <- file.path(root,'results/epidemiology'); figs <- file.path(root,'figures/epidemiology'); tabs <- file.path(out,'tables')
for(p in c(out,figs,tabs)) dir.create(p, recursive=TRUE, showWarnings=FALSE)
for(f in c('01_preprocess.R','02_associations.R','03_secondary.R','04_tables_figures.R','05_manuscript.R','validate_epidemiology.R')) {
 message('Running ',f); source(file.path(root,'R',f), local=TRUE, encoding='UTF-8')
}
if('--report' %in% args) {
 q <- Sys.which('quarto')
 if(q=='') q <- 'C:/Program Files/RStudio/resources/app/bin/quarto/bin/quarto.exe'
 if(!file.exists(q)) stop('Quarto is required for --report')
 Sys.setenv(R_ENVIRON_USER=if(.Platform$OS.type=='windows') 'NUL' else '/dev/null', R_PROFILE_USER=if(.Platform$OS.type=='windows') 'NUL' else '/dev/null', QUARTO_R=file.path(R.home('bin'),'Rscript.exe'))
 if(.Platform$OS.type!='windows') Sys.setenv(QUARTO_R=file.path(R.home('bin'),'Rscript'))
 code <- system2(q,c('render',shQuote(file.path(root,'report/Heat_Pregnancy_Manuscript.qmd'))))
 if(code!=0) stop('Quarto rendering failed')
 file.copy(file.path(root,'report/Heat_Pregnancy_Manuscript.docx'),file.path(root,'output/docx/heat_pregnancy_simulation_report.docx'),overwrite=TRUE)
}
message('Applied epidemiology analysis and validation complete')
