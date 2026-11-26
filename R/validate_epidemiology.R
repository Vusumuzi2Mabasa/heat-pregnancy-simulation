# Fail loudly on analysis inconsistencies; these are software checks, not study outcomes.
stopifnot(nrow(d)==2400,length(unique(d$id))==600,nrow(preg)==600,sum(preg$complication)==86)
stopifnot(!'vascular_intercept_truth' %in% names(d),!'vascular_intercept_truth' %in% names(preg))
stopifnot(all(vapply(models,function(x)!lme4::isSingular(x),logical(1))),all(is.finite(heat_rows$estimate)),all(heat_rows$conf.low<heat_rows$conf.high))
stopifnot(all(diag(sandwich::vcovHC(clinical_model,type='HC3'))>0),all(fitted(clinical_model)>0 & fitted(clinical_model)<1))
stopifnot(all(separation_check(cbind(1,c(-2,-1,1,2)),c(0,0,1,1))))
stopifnot(!any(separation_check(cbind(1,c(-1,-1,1,1)),c(0,1,0,1))))
qs<-separation_check(cbind(1,c(-1,0,0,1)),c(0,0,1,1));stopifnot(!qs['complete'],qs['quasi_or_complete'])
stopifnot(functional$df[1]==2,interaction_tests$df[1]==1)
# Independent direct calculation of the primary fixed-effect coefficient and t CI.
z<-coef(summary(primary))['heat5',];direct_interval<-z['Estimate']+c(-1,1)*qt(.975,z['df'])*z['Std. Error']
stopifnot(max(abs(direct_interval-c(h$conf.low,h$conf.high)))<1e-10)
# Entire cluster blocks and unique bootstrap draw labels were used; audit all statuses.
stopifnot(nrow(boot_table)==config$mediation_bootstrap_replicates,sum(valid)>=.95*nrow(boot_table),!anyNA(status))
stopifnot(!any(grepl('@@',readLines('README.md',warn=FALSE))),length(list.files(figs,pattern='\\.png$'))==9)
# Preserve the legacy input and output bytes, tolerating only Git text line endings.
manifest<-jsonlite::fromJSON('provenance/original_results_sha256.json')
sha256<-function(p,canonical=FALSE){b<-readBin(p,'raw',n=file.info(p)$size);if(canonical)b<-charToRaw(gsub('\r\n','\n',rawToChar(b),fixed=TRUE));digest::digest(b,algo='sha256',serialize=FALSE)}
preserved<-vapply(names(manifest$files),function(p){expected<-manifest$files[[p]];observed<-sha256(p);if(observed==expected)return(TRUE);canonical<-manifest$canonical_text_sha256[[p]];!is.null(canonical) && sha256(p,TRUE)==canonical},logical(1))
stopifnot(all(preserved))
jsonlite::write_json(list(status='passed',women=600,visits=2400,events=86,original_scientific_files_checked=length(preserved),original_files_preserved=all(preserved),primary_interval_reproduced=TRUE,logistic_separation_toy_checks='complete quasi and overlap passed',exposure_reconstruction='80 windows checked',figures=9,tables=length(tables),bootstrap_valid=sum(valid)),file.path(out,'validation.json'),pretty=TRUE,auto_unbox=TRUE)
message('All analysis checks passed; ',length(preserved),' original scientific files preserved')
