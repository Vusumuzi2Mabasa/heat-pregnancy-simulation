library(lmerTest)
adjust <- 'factor(visit) + age + deprivation + ancestry + baseline_distress + site + season_sin + season_cos + g_centered'
fit_lmm <- function(y,data=d,extra='',REML=TRUE) lmerTest::lmer(as.formula(paste(y,'~ heat5 +',adjust,extra,'+ (1|id)')),data=data,REML=REML)
tidy_lmm <- function(x,conf.int=TRUE,conf.level=.95,exponentiate=FALSE,...) {
 a<-as.data.frame(coef(summary(x)));q<-qt((1+conf.level)/2,a$df)
 data.frame(term=rownames(a),estimate=a$Estimate,std.error=a$`Std. Error`,df=a$df,statistic=a$`t value`,p.value=a$`Pr(>|t|)`,conf.low=a$Estimate-q*a$`Std. Error`,conf.high=a$Estimate+q*a$`Std. Error`,row.names=NULL)
}
models <- setNames(lapply(c('depression_z','vascular_z','anxiety_z'),fit_lmm),c('depression','vascular','anxiety'))
primary <- models$depression
heat_rows <- do.call(rbind,lapply(names(models),function(n){a<-tidy_lmm(models[[n]]);a<-a[a$term=='heat5',];a$outcome<-n;a$analysis<-'Fully adjusted';a}))
heat_rows$secondary_Holm_p <- NA_real_
heat_rows$secondary_Holm_p[heat_rows$outcome!='depression'] <- p.adjust(heat_rows$p.value[heat_rows$outcome!='depression'],method='holm')
minimal <- lmerTest::lmer(depression_z~heat5+factor(visit)+(1|id),data=d)
ant_fit <- fit_lmm('depression_z',antenatal)
get_heat <- function(m,label){x<-tidy_lmm(m);x<-x[x$term=='heat5',];x$analysis<-label;x}
sensitivity <- rbind(get_heat(minimal,'Visit adjusted only'),get_heat(primary,'Fully adjusted'),get_heat(ant_fit,'Antenatal visits only'))
# Within-person sensitivity separates time-varying heat from person-average heat.
d$mean_heat5 <- ave(d$heat5,d$id,FUN=mean);d$within_heat5 <- d$heat5-d$mean_heat5
within_fit <- lmerTest::lmer(as.formula(paste('depression_z ~ mean_heat5 + within_heat5 +',adjust,'+ (1|id)')),data=d)
within_rows <- tidy_lmm(within_fit);within_rows <- within_rows[within_rows$term %in% c('mean_heat5','within_heat5'),]
write.csv(within_rows,file.path(out,'within_between_associations.csv'),row.names=FALSE)
write.csv(heat_rows,file.path(out,'heat_mean_differences.csv'),row.names=FALSE);write.csv(sensitivity,file.path(out,'depression_sensitivity.csv'),row.names=FALSE)
write.csv(tidy_lmm(primary),file.path(out,'depression_coefficients.csv'),row.names=FALSE)
clinical_model <- glm(complication~antenatal_heat5+age+deprivation+ancestry+baseline_distress+site+antenatal_season_sin+antenatal_season_cos,data=preg,family=binomial())
robust_tidy <- function(x,conf.int=TRUE,conf.level=.95,exponentiate=FALSE,...) {
 b<-coef(x);se<-sqrt(diag(sandwich::vcovHC(x,type='HC3')));q<-qnorm((1+conf.level)/2)
 ans<-data.frame(term=names(b),estimate=b,std.error=se,statistic=b/se,p.value=2*pnorm(-abs(b/se)),conf.low=b-q*se,conf.high=b+q*se,row.names=NULL)
 if(exponentiate) ans[c('estimate','conf.low','conf.high')]<-lapply(ans[c('estimate','conf.low','conf.high')],exp)
 ans
}
clinical_terms <- robust_tidy(clinical_model,exponentiate=TRUE)
write.csv(clinical_terms,file.path(out,'clinical_odds_ratios.csv'),row.names=FALSE)
separation_check <- function(X,y) {
 A<-X*(2*y-1);B<-cbind(A,-A);n<-nrow(B);p<-ncol(B)
 complete<-lpSolve::lp('min',rep(0,p),B,rep('>=',n),rep(1,n))$status==0
 quasi<-lpSolve::lp('min',rep(0,p),rbind(B,colSums(B)),rep('>=',n+1),c(rep(0,n),1))$status==0
 c(complete=complete,quasi_or_complete=quasi)
}
sep <- separation_check(model.matrix(clinical_model),preg$complication)
stopifnot(!any(sep),clinical_model$converged,clinical_model$rank==length(coef(clinical_model)))
# Compare fixed effects with ML, not REML.
linear_ml <- update(primary,REML=FALSE)
spline_heat <- update(linear_ml,.~.-heat5+splines::ns(heat5,df=3))
spline_age <- update(linear_ml,.~.-age+splines::ns(age,df=3))
lr <- function(a,b) {z<-anova(a,b);c(statistic=z$Chisq[2],df=attr(logLik(b),'df')-attr(logLik(a),'df'),p=z$`Pr(>Chisq)`[2])}
functional <- rbind(data.frame(term='Heat',t(lr(linear_ml,spline_heat))),data.frame(term='Age surrogate',t(lr(linear_ml,spline_age))))
write.csv(functional,file.path(out,'functional_form_tests.csv'),row.names=FALSE)
v <- car::vif(primary)
vif_table <- if(is.matrix(v)) data.frame(term=rownames(v),v,row.names=NULL) else data.frame(term=names(v),GVIF=v,Df=1,adjusted_GVIF=sqrt(v))
write.csv(vif_table,file.path(out,'collinearity.csv'),row.names=FALSE)
diag_rows <- do.call(rbind,lapply(names(models),function(n){m<-models[[n]];u<-as.numeric(lme4::VarCorr(m)$id[1,1]);data.frame(outcome=n,singular=lme4::isSingular(m),participant_variance=u,residual_variance=sigma(m)^2,ICC=u/(u+sigma(m)^2),convergence_messages=paste(m@optinfo$conv$lme4$messages,collapse='; '))}))
write.csv(diag_rows,file.path(out,'mixed_model_diagnostics.csv'),row.names=FALSE)
residuals_table <- data.frame(id=d$id,visit=d$visit,fitted=fitted(primary),residual=residuals(primary),temperature_C=d$temperature_C)
write.csv(residuals_table,file.path(out,'depression_residuals.csv'),row.names=FALSE)
clinical_influence <- data.frame(id=preg$id,cooks_distance=cooks.distance(clinical_model),leverage=hatvalues(clinical_model),fitted=fitted(clinical_model),observed=preg$complication)
write.csv(clinical_influence,file.path(out,'clinical_influence.csv'),row.names=FALSE)
largest<-order(clinical_influence$cooks_distance,decreasing=TRUE)[1:5]
clinical_sensitivity<-update(clinical_model,data=preg[-largest,]);cs<-robust_tidy(clinical_sensitivity,exponentiate=TRUE);cs<-cs[cs$term=='antenatal_heat5',]
cs$n<-nrow(preg)-5;cs$events<-sum(preg$complication[-largest]);write.csv(cs,file.path(out,'clinical_influence_sensitivity.csv'),row.names=FALSE)
saveRDS(list(continuous=models,clinical=clinical_model,within=within_fit,antenatal=ant_fit),file.path(out,'models.rds'),compress='xz')
diagnostics<-list(clinical_converged=clinical_model$converged,clinical_rank=clinical_model$rank,clinical_coefficients=length(coef(clinical_model)),clinical_separation=as.list(sep),clinical_cooks_above_4n=sum(clinical_influence$cooks_distance>4/nrow(preg)),all_primary_records_retained=TRUE,ML_used_for_fixed_effect_tests=TRUE)
jsonlite::write_json(diagnostics,file.path(out,'model_diagnostics.json'),pretty=TRUE,auto_unbox=TRUE)
