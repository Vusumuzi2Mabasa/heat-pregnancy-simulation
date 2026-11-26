library(ggplot2);library(gtsummary);library(flextable)
fp <- function(p)ifelse(is.na(p),'',ifelse(p<.001,'<0.001',sprintf('%.3f',p)))
ci_text <- function(e,l,u)sprintf('%.2f (%.2f, %.2f)',e,l,u)
table_data <- base[c('complication','site','G','age','deprivation','ancestry','baseline_distress','temperature_C','depression_z')]
table_data$complication <- factor(table_data$complication,0:1,c('No complication','Complication'));table_data$G<-factor(table_data$G)
tbl1 <- tbl_summary(table_data,by=complication,missing='no',statistic=all_continuous()~'{mean} ({sd})',digits=all_continuous()~2,label=list(site~'Fictional site',G~'Fictional allele dosage',age~'Age surrogate',deprivation~'Deprivation surrogate',ancestry~'Ancestry score',baseline_distress~'Baseline distress',temperature_C~'First visit heat in degrees C',depression_z~'First depression outcome')) |> add_overall() |> bold_labels()
tbl_primary <- tbl_regression(primary,tidy_fun=tidy_lmm,intercept=FALSE,label=list(heat5~'Heat per 5 degrees C',`factor(visit)`~'Sample visit',age~'Age surrogate',deprivation~'Deprivation surrogate',ancestry~'Ancestry score',baseline_distress~'Baseline distress',site~'Fictional site',season_sin~'Season sine',season_cos~'Season cosine',g_centered~'Genotype per allele'),estimate_fun=label_style_sigfig(digits=3),pvalue_fun=label_style_pvalue(digits=3)) |> modify_header(estimate='**Mean difference**') |> bold_labels()
tbl_clinical <- tbl_regression(clinical_model,exponentiate=TRUE,tidy_fun=robust_tidy,label=list(antenatal_heat5~'Antenatal heat per 5 degrees C',age~'Age surrogate',deprivation~'Deprivation surrogate',ancestry~'Ancestry score',baseline_distress~'Baseline distress',site~'Fictional site',antenatal_season_sin~'Mean season sine',antenatal_season_cos~'Mean season cosine'),pvalue_fun=label_style_pvalue(digits=3)) |> bold_labels()
visits <- do.call(rbind,lapply(levels(d$visit),function(v){z<-d[d$visit==v,];data.frame(Visit=v,N=nrow(z),Heat=sprintf('%.2f (%.2f)',mean(z$temperature_C),sd(z$temperature_C)),Depression=sprintf('%.2f (%.2f)',mean(z$depression_z),sd(z$depression_z)),Vascular=sprintf('%.2f (%.2f)',mean(z$vascular_z),sd(z$vascular_z)),Anxiety=sprintf('%.2f (%.2f)',mean(z$anxiety_z),sd(z$anxiety_z)))}))
heat_display <- data.frame(Outcome=heat_rows$outcome,`Adjusted mean difference (95% CI)`=ci_text(heat_rows$estimate,heat_rows$conf.low,heat_rows$conf.high),P=fp(heat_rows$p.value),`Secondary Holm P`=fp(heat_rows$secondary_Holm_p),check.names=FALSE)
sens_display <- data.frame(Analysis=sensitivity$analysis,`Mean difference (95% CI)`=ci_text(sensitivity$estimate,sensitivity$conf.low,sensitivity$conf.high),P=fp(sensitivity$p.value),check.names=FALSE)
secondary_display <- data.frame(Component=mediation$component,`Mean difference (95% bootstrap CI)`=ci_text(mediation$estimate,mediation$conf.low,mediation$conf.high),check.names=FALSE)
diagnostic_display <- data.frame(Check=c('Depression random intercept fit','Clinical logistic fit','Complete or quasi separation','Heat functional form','Age functional form','Participant bootstrap','Influence sensitivity'),Finding=c(sprintf('ICC %.2f; singular: %s',diag_rows$ICC[1],diag_rows$singular[1]),sprintf('Converged; design rank %d of %d',clinical_model$rank,length(coef(clinical_model))),'Absent by linear programming',paste('ML spline comparison p',fp(functional$p[1])),paste('ML spline comparison p',fp(functional$p[2])),sprintf('%d of %d converged resamples',sum(valid),B),sprintf('Omitting five largest Cook distances: clinical heat OR %s',ci_text(cs$estimate,cs$conf.low,cs$conf.high)) ))
make_ft <- function(x){
 ft<-if(inherits(x,'gtsummary'))as_flex_table(x) else flextable(x)
 ft<-ft |> font(fontname='Calibri',part='all') |> fontsize(size=9,part='all') |> padding(padding=3,part='all') |> bg(bg='#E6EEF2',part='header') |> color(color='black',part='all') |> autofit() |> set_table_properties(layout='autofit',width=1)
 border<-officer::fp_border(color='#D9D9D9',width=.5)
 ft<-border_outer(ft,border=border);ft<-border_inner_h(ft,border=border);ft<-border_inner_v(ft,border=border);paginate(ft,init=TRUE,hdr_ftr=TRUE)
}
tables <- list(baseline=tbl1,visits=visits,primary=tbl_primary,heat=heat_display,clinical=tbl_clinical,sensitivity=sens_display,mediation=secondary_display,diagnostics=diagnostic_display)
write_md <- function(x,name){
 z<-if(inherits(x,'gtsummary'))as_tibble(x,col_labels=TRUE) else x
 cell<-function(v){v<-as.character(v);v[is.na(v)]<-'';gsub('\n','<br>',v,fixed=TRUE)}
 z[]<-lapply(z,cell)
 lines<-c(paste0('| ',paste(cell(names(z)),collapse=' | '),' |'),paste0('| ',paste(rep('---',ncol(z)),collapse=' | '),' |'),apply(z,1,function(r)paste0('| ',paste(r,collapse=' | '),' |')))
 writeLines(lines,file.path(tabs,paste0(name,'.md')),useBytes=TRUE);write.csv(z,file.path(tabs,paste0(name,'_display.csv')),row.names=FALSE)
}
for(n in names(tables)){saveRDS(make_ft(tables[[n]]),file.path(tabs,paste0(n,'_flextable.rds')));write_md(tables[[n]],n)}
saveRDS(tbl1$table_body,file.path(tabs,'baseline_gtsummary_body.rds'));saveRDS(tbl_primary$table_body,file.path(tabs,'primary_gtsummary_body.rds'));saveRDS(tbl_clinical$table_body,file.path(tabs,'clinical_gtsummary_body.rds'))
theme_set(theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),plot.title=element_text(face='bold'),plot.subtitle=element_text(size=10)))
save_plot<-function(name,p,width=8,height=4.5){ggsave(file.path(figs,paste0(name,'.png')),p,width=width,height=height,dpi=300,bg='white');ggsave(file.path(figs,paste0(name,'.pdf')),p,width=width,height=height,device=cairo_pdf)}
p1<-ggplot(d,aes(visit,temperature_C,fill=site))+geom_boxplot(position=position_dodge(.8),outlier.size=.3)+scale_x_discrete(labels=c('Week 14','Week 20','Week 26','Postpartum'))+scale_fill_manual(values=c('#167D9A','#B56528','#667C54'))+labs(x='Sample timing',y='Preceding seven day mean Tmax (degrees C)',fill='Fictional site',title='Heat exposure across the four sample visits',subtitle='All weather and participant values are simulated')
save_plot('01_exposure',p1,8,4)
traj<-aggregate(depression_z~id+visit,d,mean)
p2<-ggplot(traj,aes(visit,depression_z,group=id))+geom_line(alpha=.035,colour='#558499')+geom_pointrange(data=visit_means,aes(x=visit,y=emmean,ymin=lower.CL,ymax=upper.CL,group=1),inherit.aes=FALSE,colour='#A34B29',linewidth=.6)+geom_line(data=visit_means,aes(visit,emmean,group=1),inherit.aes=FALSE,colour='#A34B29',linewidth=.7)+scale_x_discrete(labels=c('Week 15','Week 21','Week 27','7 weeks postpartum'))+labs(x='Outcome timing one week after each sample',y='Continuous depression surrogate',title='Repeated depression outcomes and adjusted visit means',subtitle='Thin lines: observed synthetic participants; orange: adjusted mean and 95% CI')
save_plot('02_trajectories',p2,8,4.2)
f<-heat_rows;f$outcome<-factor(f$outcome,levels=c('anxiety','vascular','depression'),labels=c('Anxiety secondary','Vascular secondary','Depression primary'))
p3<-ggplot(f,aes(estimate,outcome))+geom_vline(xintercept=0,linetype=2,colour='grey50')+geom_errorbar(aes(xmin=conf.low,xmax=conf.high),orientation='y',width=.15,colour='#167D9A')+geom_point(size=2.5,colour='#167D9A')+labs(x='Adjusted mean difference per 5 degrees C with 95% CI',y=NULL,title='Heat associations with continuous health surrogates',subtitle='Participant random intercept models; no expression module adjustment')
save_plot('03_heat_associations',p3,8,3.3)
f<-clinical_terms[clinical_terms$term %in% c('antenatal_heat5','age','deprivation','baseline_distress'),];f$label<-c(antenatal_heat5='Antenatal heat per 5 degrees C',age='Age surrogate per unit',deprivation='Deprivation surrogate per unit',baseline_distress='Baseline distress per unit')[f$term]
p4<-ggplot(f,aes(estimate,reorder(label,estimate)))+geom_vline(xintercept=1,linetype=2,colour='grey50')+geom_errorbar(aes(xmin=conf.low,xmax=conf.high),orientation='y',width=.15,colour='#167D9A')+geom_point(size=2.5,colour='#167D9A')+scale_x_log10()+labs(x='Adjusted odds ratio with HC3 95% CI',y=NULL,title='Associations with the generated pregnancy complication',subtitle='600 pregnancies and 86 events; antenatal exposure windows only')
save_plot('04_clinical_associations',p4,8,3.5)
qq<-data.frame(theoretical=qnorm(ppoints(nrow(d))),observed=sort(residuals(primary)/sigma(primary)))
p5<-ggplot(qq,aes(theoretical,observed))+geom_point(size=.6,alpha=.45,colour='#167D9A')+geom_abline(slope=1,intercept=0,linetype=2)+labs(x='Theoretical standard normal quantile',y='Ordered scaled conditional residual',title='Depression model residual normality assessment',subtitle='Conditional residuals after accounting for participant random intercepts')
save_plot('05_residual_qq',p5,7,4)
residuals_table$phase<-ifelse(as.integer(residuals_table$visit)==4,'Postpartum','Antenatal')
p6<-ggplot(residuals_table,aes(fitted,residual,colour=phase))+geom_point(size=.65,alpha=.3)+geom_hline(yintercept=0,linetype=2)+scale_colour_manual(values=c('#167D9A','#B56528'))+labs(x='Conditional fitted depression surrogate',y='Conditional residual',colour='Phase',title='Depression model residual variance assessment',subtitle='Antenatal and postpartum observations shown separately')
save_plot('06_residual_fitted',p6,7,4)
p7<-ggplot(gslopes,aes(factor(G),estimate))+geom_pointrange(aes(ymin=lower.CL,ymax=upper.CL),colour='#167D9A')+geom_hline(yintercept=0,linetype=2)+labs(x='Fictional allele dosage',y='Depression mean difference per 5 degrees C',title='Exploratory genetic modification of the heat association',subtitle='Separate slopes from one interaction model; 95% confidence intervals')
save_plot('07_genotype_interaction',p7,7,4)
p8<-ggplot(rho_sens,aes(rho,estimate))+geom_ribbon(aes(ymin=conf.low,ymax=conf.high),fill='#D6E7EE')+geom_line(colour='#167D9A',linewidth=.8)+geom_hline(yintercept=0,linetype=2)+geom_vline(xintercept=0,linetype=3)+labs(x='Assumed mediator and outcome residual correlation',y='Exploratory indirect mean difference per 5 degrees C',title='Sensitivity of the expression pathway decomposition',subtitle='Participant bootstrap percentile intervals; assumptions determine causal interpretation')
save_plot('08_mediation_sensitivity',p8,8,4)
# Random-intercept normality is a separate assumption from residual normality.
u<-lme4::ranef(primary)$id[,1];uq<-data.frame(theoretical=qnorm(ppoints(length(u))),observed=sort(as.numeric(scale(u))))
save_plot('09_random_intercept_qq',ggplot(uq,aes(theoretical,observed))+geom_point(size=.7,colour='#167D9A')+geom_abline(slope=1,intercept=0,linetype=2)+labs(x='Theoretical normal quantile',y='Scaled conditional random intercept',title='Participant random intercept distribution',subtitle='Conditional estimates are shrunken; diagnostic interpretation is approximate'),7,4)
packages<-c(required,'rmarkdown')
write.csv(data.frame(package=packages,version=vapply(packages,function(x)as.character(packageVersion(x)),character(1))),file.path(out,'R_packages.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
