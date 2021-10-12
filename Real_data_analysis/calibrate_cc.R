library(survival)
library(survey)
cohort.study = survival::nwtco
cohort.study$age = as.numeric(scale(cohort.study$age/12))
cohort.study$age1 = as.numeric(cohort.study$age/12)
#cohort.study$age1 = as.numeric(scale(cohort.study$age1))
cohort.study$age2 = if_else(cohort.study$age1 >= 1, cohort.study$age1-1, 0 )
#cohort.study$age2 = as.numeric(scale(cohort.study$age2))
cohort.study$histol = factor(cohort.study$histol)
#cohort.study$instit <- factor(cohort.study$instit)

cohort.study$stage.collapse = if_else(cohort.study$stage==4 | cohort.study$stage==3, 1, 0)
cohort.study$stage.collapse = factor(cohort.study$stage.collapse)
cohort.study$stage <- factor(cohort.study$stage)
beta_true = glm(rel~histol*stage + stage*age + histol*age, data = cohort.study, family = quasibinomial())
#cohort.study$histol <- factor(cohort.study$histol)
#cohort.study$instit <- factor(cohort.study$instit)
#cohort.study$stage <- factor(cohort.study$stage)
results.npmle.survey = matrix(NA, 1000,25)



for(i in 1:1000)
{
  calibrate.fail = 0
  glm.fail = 0
  cohort.study = survival::nwtco
  #cohort.study$histol <- factor(cohort.study$histol)
  cohort.study$instit <- factor(cohort.study$instit)
  cohort.study$stage <- factor(cohort.study$stage)
  cohort.study$age  = as.numeric(cohort.study$age/12)
  cohort.study$age.stratify = if_else(cohort.study$age >= median(cohort.study$age), 1, 0 )
  cohort.strata.table <- table(cohort.study$rel, cohort.study$instit)
  #cohort.study$age1 = as.numeric(cohort.study$age/12)
  #cohort.study$age1 = as.numeric(scale(cohort.study$age1))
  #cohort.study$age2 = if_else(cohort.study$age >= 1, cohort.study$age1-1, 0 )
  
  
  sampling.fraction.cases.strata1 = 1
  sampling.fraction.cases.strata2 = 1
  sampling.fraction.controls.strata1 = sum(cohort.study$rel)/sum(1 - cohort.study$rel)
  sampling.fraction.controls.strata2 = sum(cohort.study$rel)/sum(1 - cohort.study$rel)
  
  cohort.study$sampling.fraction = NA
  cohort.study$sampling.fraction[which(cohort.study$rel == 1 & cohort.study$instit == 1)] = sampling.fraction.cases.strata1
  cohort.study$sampling.fraction[which(cohort.study$rel == 1 & cohort.study$instit == 2)] = sampling.fraction.cases.strata2
  cohort.study$sampling.fraction[which(cohort.study$rel == 0 & cohort.study$instit == 1)] = sampling.fraction.controls.strata1
  cohort.study$sampling.fraction[which(cohort.study$rel == 0 & cohort.study$instit == 2)] = sampling.fraction.controls.strata2
  set.seed(i)
  R.cases.strata1 = rbinom(length(which(cohort.study$rel == 1 & cohort.study$instit == 1)), 1, sampling.fraction.cases.strata1)
  set.seed(i)
  R.cases.strata2 = rbinom(length(which(cohort.study$rel == 1 & cohort.study$instit == 2)), 1, sampling.fraction.cases.strata2)
  set.seed(i)
  R.controls.strata1 = rbinom(length(which(cohort.study$rel == 0 & cohort.study$instit == 1)), 1, sampling.fraction.controls.strata1)
  set.seed(i)
  R.controls.strata2 = rbinom(length(which(cohort.study$rel == 0 & cohort.study$instit == 2)), 1, sampling.fraction.controls.strata2)
  
  CC.cases.strata1 = cbind(cohort.study[which(cohort.study$rel == 1 & cohort.study$instit == 1),], R=R.cases.strata1) 
  CC.cases.strata2 = cbind(cohort.study[which(cohort.study$rel == 1 & cohort.study$instit == 2),], R=R.cases.strata2) 
  CC.controls.strata1 = cbind(cohort.study[which(cohort.study$rel == 0 & cohort.study$instit == 1),], R=R.controls.strata1) 
  CC.controls.strata2 = cbind(cohort.study[which(cohort.study$rel == 0 & cohort.study$instit == 2),], R=R.controls.strata2)
  CC.study = rbind(CC.cases.strata1,CC.cases.strata2, CC.controls.strata1, CC.controls.strata2)
  #CC.study = CC.study[CC.study$R==1, ]
  
  cohort.study = CC.study
  #cohort.study$histol[which(cohort.study$R == 0)] = NA
  
  
  
  
  # cohort.study = nwtco
  # cohort.study$age = as.numeric(cohort.study$age/12)
  # cohort.study$age1 = as.numeric(cohort.study$age/12)
  # cohort.study$age2 = if_else(cohort.study$age1 >= 1, cohort.study$age1-1, 0 )
  # cohort.study$age.cat = if_else(cohort.study$age>=1, 1,0)
  # 
  # cohort.study$stage.collapse = if_else(cohort.study$stage==4 | cohort.study$stage==3, 1, 0)
  # cohort.study$stage.collapse = factor(cohort.study$stage.collapse)
  # 
  # #cohort.study$age.cat = cut(cohort.study$age, breaks = c(min(cohort.study$age),median(cohort.study$age), max(cohort.study$age)), include.lowest = T)
  # cohort.study$strata = 1
  # cohort.study$strata[which(cohort.study$rel==0 & cohort.study$age.cat==0 & cohort.study$instit == 1 & (cohort.study$stage==1|cohort.study$stage==2))] = 2
  # cohort.study$strata[which(cohort.study$rel==0 & cohort.study$age.cat==1 & cohort.study$instit == 1 & (cohort.study$stage==1|cohort.study$stage==2))] = 3
  # cohort.study$strata[which(cohort.study$rel==0 & cohort.study$age.cat==1 & cohort.study$instit == 1 & (cohort.study$stage==3|cohort.study$stage==4))] = 4
  # 
  # cohort.study$R = 0
  # cohort.study$R[which(cohort.study$strata == 1)] = 1
  # set.seed(i)
  # cohort.study$R[which(cohort.study$strata == 2)] = rbinom(length(which(cohort.study$strata == 2)), 1, 90/length(which(cohort.study$strata == 2)))
  # set.seed(i)
  # cohort.study$R[which(cohort.study$strata == 3)] = rbinom(length(which(cohort.study$strata == 3)), 1, 120/length(which(cohort.study$strata == 3)))
  # set.seed(i)
  # cohort.study$R[which(cohort.study$strata == 4)] = rbinom(length(which(cohort.study$strata == 4)), 1, 90/length(which(cohort.study$strata == 4)))
  # 
  # cohort.study$histol[which(cohort.study$R == 0)] = NA
  # 
  cohort.study$age.impute = if_else(cohort.study$age<=1, 0, 1)
  cohort.study$age.median = if_else(cohort.study$age<=median(cohort.study$age), 0, 1)
  cohort.study$stage.impute = if_else(cohort.study$stage==4 , 1, 0)
  # 
  #cohort.study$histol <- factor(cohort.study$histol)
  #cohort.study$instit <- factor(cohort.study$instit)
  cohort.study$stage <- factor(cohort.study$stage)
  ##cohort.study$stage.impute <- factor(cohort.study$stage.impute)
  cohort.study$age = as.numeric(cohort.study$age)
  cohort.study$study = as.factor(cohort.study$study)
  cohort.study$histol = if_else(cohort.study$histol == 1, 0, 1, missing = NULL)
  #est.wts.data = estWeights(cohort.study, formula=~instit*stage.impute+age.impute + study,
  #              subset=~I(!is.na(histol)))
  #est.wts.fit = svyglm(histol~instit*stage.impute + age.impute + study, design = est.wts.data, family=quasibinomial())
  cohort.study$in.subsample = TRUE
  cohort.study$in.subsample[which(cohort.study$R==0)] = FALSE
  #cohort.study$in.subsample = as.factor(cohort.study$in.subsample)
  imp.model = glm(histol ~ instit + age.impute + stage.impute*study, data = cohort.study, subset=in.subsample, family = binomial())
  #imp.model = glm(histol ~ instit*stage + age + study, data = cohort.study, subset=in.subsample, family = quasibinomial())
  #imphistol = predict.glm(imp.model, newdata = cohort.study, type = "response")
  
  cohort.study$imphist = predict(imp.model, newdata=cohort.study, type="response")
  cohort.study$imphist[which(cohort.study$in.subsample ==T)] = cohort.study$histol[which(cohort.study$in.subsample == T)]
  
  phaseI.fit = glm(rel~imphist*stage + stage*scale(age) + imphist*scale(age)  , data = cohort.study, family = quasibinomial())
  df_beta = dfbetas(phaseI.fit) 
  colnames(df_beta) = paste0("z", seq(1, ncol(df_beta), 1))
  #df_beta = resid(phaseI.fit, type="deviance") + 1
  cohort.study = as.data.frame(cbind(cohort.study, df_beta))
  #pop.totals = c(z1=sum(cohort.study$z1), z2=sum(cohort.study$z2), z3=sum(cohort.study$z3), z4=sum(cohort.study$z4), z5=sum(cohort.study$z5))
  two_phase_design = twophase(id = list(~seqno, ~seqno), subset=~in.subsample, data=cohort.study, strata = list (NULL, ~interaction(instit , rel)))
  
  result = tryCatch({
    calibrate.design = calibrate(two_phase_design, phase=2,formula =~ z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8 + z9 + z10 + z11 + z12 + z13 + instit*rel*stage,
                                 calfun=c("raking"))
  }, warning = function(warning_condition) {
    print("warning")
    calibrate.fail = 1
  }, error = function(error_condition) {
    print("error")
    calibrate.fail = 1
  })
  
  result2 = tryCatch({
    final.fit = svyglm(rel~histol*stage + stage*scale(age) + histol*scale(age), design = calibrate.design, family=quasibinomial())
  }, warning = function(warning_condition) {
    print("warning")
    glm.fail = 1
  }, error = function(error_condition) {
    print("error")
    glm.fail = 1
  })
  
  
  
  if(is.numeric(result) != TRUE & is.numeric(result2) != TRUE)
  {
    results.npmle.survey[i,] = c(final.fit$coefficients[-1], diag(vcov(final.fit))[-1], as.numeric(final.fit$converged))
  }else{
    results.npmle.survey[i,] = NA
  }
  print(i)
}
saveRDS(results.npmle.survey, "/Users/prosenjitkundu/Desktop/Code_biometrics_github/Real_data_analysis/calibrate_cc.rds")
