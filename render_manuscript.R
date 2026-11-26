# Regenerate the narrative from saved results without refitting models.
if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','English_United States.utf8'))
entry<-grep('^--file=',commandArgs(),value=TRUE)
root<-normalizePath(dirname(sub('^--file=','',entry[1])),winslash='/');setwd(root)
out<-file.path(root,'results/epidemiology');figs<-file.path(root,'figures/epidemiology');tabs<-file.path(out,'tables')
read_result<-function(n)read.csv(file.path(out,paste0(n,'.csv')))
heat_rows<-read_result('heat_mean_differences');clinical_terms<-read_result('clinical_odds_ratios')
mediation<-read_result('mediation_decomposition');within_rows<-read_result('within_between_associations');sensitivity<-read_result('depression_sensitivity')
diag_rows<-read_result('mixed_model_diagnostics');functional<-read_result('functional_form_tests');vif_table<-read_result('collinearity');cs<-read_result('clinical_influence_sensitivity')
gslopes<-read_result('genotype_specific_slopes');interaction_tests<-read_result('interaction_tests')
audit<-jsonlite::fromJSON(file.path(out,'bootstrap_audit.json'));B<-audit$requested;valid<-c(rep(TRUE,audit$valid),rep(FALSE,B-audit$valid))
fp<-function(p)ifelse(is.na(p),'',ifelse(p<.001,'<0.001',sprintf('%.3f',p)))
tables<-setNames(vector('list',8),c('baseline','visits','primary','heat','clinical','sensitivity','mediation','diagnostics'))
for(n in names(tables)){
 p<-file.path(tabs,paste0(n,'_flextable.rds'));ft<-readRDS(p);ft<-flextable::paginate(ft,init=TRUE,hdr_ftr=TRUE);saveRDS(ft,p)
}
source('R/05_manuscript.R',encoding='UTF-8')
q<-Sys.which('quarto');if(q=='')q<-'C:/Program Files/RStudio/resources/app/bin/quarto/bin/quarto.exe'
Sys.setenv(R_ENVIRON_USER=if(.Platform$OS.type=='windows')'NUL' else '/dev/null',R_PROFILE_USER=if(.Platform$OS.type=='windows')'NUL' else '/dev/null',QUARTO_R=file.path(R.home('bin'),if(.Platform$OS.type=='windows')'Rscript.exe' else 'Rscript'))
stopifnot(system2(q,c('render',shQuote(file.path(root,'report/Heat_Pregnancy_Manuscript.qmd'))))==0)
stopifnot(file.copy('report/Heat_Pregnancy_Manuscript.docx','output/docx/heat_pregnancy_simulation_report.docx',overwrite=TRUE))
