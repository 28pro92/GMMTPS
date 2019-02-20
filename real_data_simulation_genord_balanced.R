library(missreg3)
library(survival)
library(OrdNor)
instit_margin = as.numeric(prop.table(margin.table(table(nwtco$instit, nwtco$histol, nwtco$stage), 1)))
hist_margin = as.numeric(prop.table(margin.table(table(nwtco$instit, nwtco$histol, nwtco$stage), 2)))
stage_margin = as.numeric(prop.table(margin.table(table(nwtco$instit, nwtco$histol, nwtco$stage), 3)))
#set.seed(1)
# Sets the marginals.
# The values are cumulative so for the first variable the first marginal will be .1, the second is .2, the third is .3, and the fourth is .4
marginal = list(cumsum(instit_margin)[-2], cumsum(hist_margin)[-2],cumsum(stage_margin)[-4])



beta_real_data = as.numeric(glm(rel~ instit + histol * factor(stage) + scale(age/12), data = nwtco, family = binomial())$coefficients)
beta_real_data[1] = -3.6
beta_real_data[2] = beta_real_data[3]
rho_nwtco = cor(cbind(nwtco$instit, nwtco$histol, nwtco$stage,scale(nwtco$age/12)))
new_corr = rbind(cbind(rho_nwtco, c(0,0,0,0)), c(0,0,0,0,1))

#cohort.study = nwtco
#cohort.study$histol <- factor(cohort.study$histol)
#cohort.study$instit <- factor(cohort.study$instit)
#cohort.study$stage <- factor(cohort.study$stage)
#cohort.study$age = as.numeric(scale(cohort.study$age/12))

#cohort.study$age.cat = cut(cohort.study$age, breaks = c(min(cohort.study$age),median(cohort.study$age), max(cohort.study$age)), include.lowest = T)
mean = c()

results.npmle.missreg = matrix(NA, 1000,21)
for(i in 1:1000)
{
  set.seed(i)
  mm.cmat = cmat.star(plist = marginal, CorrMat = new_corr, no.ord = 3, no.norm = 2)
  set.seed(i)
  data = genOrdNor(n = 10000, plist = marginal, cmat.star = mm.cmat, mean.vec = c(0,0), sd.vec = c(1,1), no.ord = 3, no.norm = 2)[, -5]
  #corr[[i]] = cor(data)
  data = as.data.frame(data)
  colnames(data) = c("instit", "histol", "stage", "age")
  data$instit = as.factor(data$instit)
  data$histol = as.factor(data$histol)
  data$stage = as.factor(data$stage)
  data.model.matrix = model.matrix(~ instit + histol *stage + age, data = data)
  set.seed(i)
  outcome = rbinom(10000, 1, expit(beta_real_data, data.model.matrix))
  data_sim = cbind(outcome, data)
  colnames(data_sim)[1] = "rel"
  mean = c(mean, mean(data_sim$rel))
  cohort.study = data_sim
  cohort.study$age.cat = cut(cohort.study$age, breaks = c(min(cohort.study$age),median(cohort.study$age), max(cohort.study$age)), include.lowest = T)
  sampling.fraction.cases.strata1 = 1
  sampling.fraction.cases.strata2 = 1
  sampling.fraction.controls.strata1 = sum(cohort.study$rel == 1 & cohort.study$instit == 2)/sum(cohort.study$rel == 0 & cohort.study$instit == 1)
  sampling.fraction.controls.strata2 = sum(cohort.study$rel == 1 & cohort.study$instit == 1)/sum(cohort.study$rel == 0 & cohort.study$instit == 2)
  
  cohort.study$sampling.fraction = NA
  cohort.study$sampling.fraction[which(cohort.study$rel == 1)] = sampling.fraction.cases.strata1
  cohort.study$sampling.fraction[which(cohort.study$rel == 0)] = sampling.fraction.controls.strata2
  
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
  CC.study = CC.study[CC.study$R==1, ]
  
  data.phase.I.missreg = cohort.study
  #data.phase.I.missreg = cbind(aggregate(data.phase.I$rel, by = list(data.phase.I$instit, data.phase.I$stage, data.phase.I$age.cat), FUN=sum), (aggregate(data.phase.I$rel, by = list(data.phase.I$instit, data.phase.I$stage, data.phase.I$age.cat), FUN=length)$x - aggregate(data.phase.I$rel, by = list(data.phase.I$instit, data.phase.I$stage, data.phase.I$age.cat), FUN=sum)$x))
  nn1 = data.phase.I.missreg[,2]
  nn0 = data.phase.I.missreg[,3]
  #colnames(data.phase.I.missreg) = c("instit", "stage","age.cat","case","control")
  data.phase.I.missreg$histol = factor(NA)
  data.phase.I.missreg$age = NA
  #data.phase.I.missreg$V3 =NA
  data.phase.I.missreg$obstype = "strata"
  data.phase.I.missreg = data.phase.I.missreg[c("rel","instit", "stage", "age","age.cat", "histol","obstype")]
  
  data.phase.II.missreg = CC.study
  
  #data.phase.II.missreg = cbind(aggregate(data.phase.II$rel, by = list(data.phase.II$instit, data.phase.II$stage, data.phase.II$age, data.phase.II$histol), FUN=sum), (aggregate(data.phase.II$rel, by = list(data.phase.II$instit, data.phase.II$stage, data.phase.II$age, data.phase.II$histol), FUN=length)$x - aggregate(data.phase.II$rel, by = list(data.phase.II$instit, data.phase.II$stage, data.phase.II$age, data.phase.II$histol), FUN=sum)$x))
  #data.phase.II.missreg$case = data.phase.II.missreg$rel
  #data.phase.II.missreg$control = 1 - data.phase.II.missreg$rel
  #colnames(data.phase.II.missreg) = c("instit", "stage", "age", "histol", "case","control")
  data.phase.II.missreg$obstype = "retro"
  data.phase.II.missreg = data.phase.II.missreg[c("rel","instit","stage", "age", "age.cat","histol","obstype")]
  
  data.missreg = rbind(data.phase.I.missreg, data.phase.II.missreg)
  #fit.npmle <- bin2stg(cbind(case,control) ~ histol*stage + age, xstrata=c("instit", "stage", "age"), data=data.missreg, xs.includes=TRUE)
  
  possible_error = tryCatch({
    fit.npmle <- bin2stg(formula = rel ~ instit + histol * stage + age, xstrata=c("instit","stage"), data=data.missreg, xs.includes=TRUE)
  }, error=function(e)e)
  
  if(length(possible_error) > 2)
  {
    print(paste0("No Error in simulation:", i))
    if(fit.npmle$error == 0)
    {
      results.npmle.missreg[i, ] = c(fit.npmle$coefficients, diag(fit.npmle$cov), fit.npmle$error)
    } else {
      
      results.npmle.missreg[i,] = rep(NA, 21) 
    }
  }
  
}

saveRDS(results.npmle.missreg, file="missreg_realdata_simulation_balanced_w_age.rds")
