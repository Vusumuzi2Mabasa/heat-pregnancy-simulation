# Two prespecified exploratory modifiers of the primary continuous association.
interaction_g <- update(primary,.~.+heat5:g_centered)
interaction_phase <- update(primary,.~.+heat5:postpartum)
interaction_tests <- rbind(data.frame(test='Heat by genotype',t(lr(linear_ml,update(interaction_g,REML=FALSE)))),data.frame(test='Heat by postpartum phase',t(lr(linear_ml,update(interaction_phase,REML=FALSE)))))
interaction_tests$Holm_p <- p.adjust(interaction_tests$p,method='holm')
write.csv(interaction_tests,file.path(out,'interaction_tests.csv'),row.names=FALSE)
gslopes <- as.data.frame(emmeans::emtrends(interaction_g,~g_centered,var='heat5',at=list(g_centered=c(-.6,.4,1.4)),lmer.df='satterthwaite'))
names(gslopes)[names(gslopes)=='heat5.trend']<-'estimate';gslopes$G<-0:2
write.csv(gslopes,file.path(out,'genotype_specific_slopes.csv'),row.names=FALSE)
phase_slopes <- as.data.frame(emmeans::emtrends(interaction_phase,~postpartum,var='heat5',at=list(postpartum=0:1),lmer.df='satterthwaite'))
write.csv(phase_slopes,file.path(out,'phase_specific_slopes.csv'),row.names=FALSE)
# Visit-specific model means are predictions at mean continuous covariates, with
# equal weighting over the three fictional sites; not observed outcome changes.
visit_means <- as.data.frame(emmeans::emmeans(primary,~visit,lmer.df='satterthwaite'))
write.csv(visit_means,file.path(out,'adjusted_visit_means.csv'),row.names=FALSE)

# Exploratory linear path decomposition. This is not a longitudinal causal
# g-formula: exposure-induced time-varying confounding is not modelled.
med_formula <- as.formula(paste('module ~ heat5*g_centered +',adjust))
outcome_formula <- as.formula(paste('depression_z ~ heat5*g_centered + module +',adjust,'+ (1|id)'))
path_fit <- function(z) {
 a<-lm(med_formula,data=z)
 b<-lme4::lmer(outcome_formula,data=z,REML=TRUE)
 mean_g<-mean(z$g_centered)
 ah<-coef(a)['heat5']+coef(a)['heat5:g_centered']*mean_g
 bh<-lme4::fixef(b)['module']
 ch<-lme4::fixef(b)['heat5']+lme4::fixef(b)['heat5:g_centered']*mean_g
 c(indirect=unname(ah*bh),direct=unname(ch),total=unname(ch+ah*bh),a=unname(ah),b=unname(bh),sigma_M=sigma(a),sigma_Y=sigma(b),singular=as.integer(lme4::isSingular(b)),converged=as.integer(is.null(b@optinfo$conv$lme4$messages)))
}
path_point <- path_fit(d)
set.seed(config$bootstrap_seed)
B <- config$mediation_bootstrap_replicates
boot <- matrix(NA_real_,nrow=B,ncol=length(path_point),dimnames=list(NULL,names(path_point)))
status <- rep('failed',B); ids<-levels(d$id);blocks<-split(d,d$id)
for(i in seq_len(B)) {
 draw<-sample(seq_along(ids),length(ids),replace=TRUE)
 z<-do.call(rbind,lapply(seq_along(draw),function(j){x<-blocks[[draw[j]]];x$id<-paste0('B',j);x}))
 z$id<-factor(z$id)
 val<-tryCatch(suppressMessages(suppressWarnings(path_fit(z))),error=function(e)e)
 if(!inherits(val,'error')){boot[i,]<-val;status[i]<-if(val['converged']==1)'valid' else 'nonconverged'}
 if(i%%50==0)message('Participant bootstrap ',i,' / ',B)
}
boot_table <- data.frame(replicate=seq_len(B),status=status,boot,check.names=FALSE)
write.csv(boot_table,file.path(out,'mediation_cluster_bootstrap.csv'),row.names=FALSE)
valid<-status=='valid' & apply(is.finite(boot),1,all)
stopifnot(mean(valid)>=.95)
ci<-apply(boot[valid,c('indirect','direct','total'),drop=FALSE],2,quantile,probs=c(.025,.975))
mediation <- data.frame(component=c('indirect','direct','total'),estimate=path_point[c('indirect','direct','total')],conf.low=ci[1,],conf.high=ci[2,],valid_bootstraps=sum(valid),requested_bootstraps=B,row.names=NULL)
write.csv(mediation,file.path(out,'mediation_decomposition.csv'),row.names=FALSE)
rho<-seq(-.6,.6,.05)
rho_sens<-do.call(rbind,lapply(rho,function(r){kap<-r*path_point['sigma_Y']/(path_point['sigma_M']*sqrt(1-r*r));v<-boot[valid,'a']*(boot[valid,'b']-r*boot[valid,'sigma_Y']/(boot[valid,'sigma_M']*sqrt(1-r*r)));q<-quantile(v,c(.025,.975));data.frame(rho=r,estimate=unname(path_point['a']*(path_point['b']-kap)),conf.low=q[1],conf.high=q[2])}))
write.csv(rho_sens,file.path(out,'mediation_residual_correlation_sensitivity.csv'),row.names=FALSE)
jsonlite::write_json(list(seed=config$bootstrap_seed,requested=B,valid=sum(valid),failed=sum(status=='failed'),nonconverged=sum(status=='nonconverged'),singular=sum(boot[valid,'singular']),resampling_unit='Woman with all four visits; duplicate draws receive distinct new identifiers',interpretation='Exploratory linear path decomposition; not proof of a causal mechanism'),file.path(out,'bootstrap_audit.json'),pretty=TRUE,auto_unbox=TRUE)
