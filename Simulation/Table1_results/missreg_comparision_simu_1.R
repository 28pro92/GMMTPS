library(mvtnorm)
source("/Users/prosenjitkundu/Desktop/Code_biometrics_github/Simulation/Table2_results/missreg3/bin2stg.R")
source("/Users/prosenjitkundu/Desktop/Code_biometrics_github/Simulation/Table2_results/missreg3/binmods.R")
source("/Users/prosenjitkundu/Desktop/Code_biometrics_github/Simulation/Table2_results/missreg3/MLEfn.R")
source("/Users/prosenjitkundu/Desktop/Code_biometrics_github/Simulation/Table2_results/missreg3/MLInf.R")
source("/Users/prosenjitkundu/Desktop/Code_biometrics_github/Simulation/Table2_results/missreg3/utilities.R")
expit <- function(beta, X)
{
  return(as.numeric(1/(1 + exp(-as.matrix(X) %*% beta))))
}

results.npmle.missreg = matrix(NA, 1000,9)
beta.true = c(-3.0, 0, log(1.3), log(1.3),  log(1.3))
for(i in 1:1000)
{
  set.seed(i)
  
  X = rmvnorm(10000, c(0,0,0), matrix(c(1,0.3,0.6,0.3,1,0.1,0.6,0.1,1),3,3))
  X1.cont = X[,1]
  X[,1]= ifelse(X[,1]<=0, 0, 1)
  #X[which(X[,1] <= 0), 1] = 0
  #X[which(X[,1] > 0), 1] = 1
  X.design = cbind(1, X, X1.cont)
  p = expit(beta.true, X.design)
  # X.design = cbind(1, X)
  # p = expit(beta.true, X.design)
  set.seed(i)
  Y = rbinom(10000, 1, p)
  
  cohort.study = data.frame(cbind(Y, X.design))
  
  #fit.full.cohort = glm(Y ~ V3 + V4 + V5, family = binomial(), data = cohort.study)
  
  
  sampling.fraction.cases.strata1 = 1
  sampling.fraction.cases.strata2 = 1
  sampling.fraction.controls.strata1 = sum(cohort.study$Y)/sum(1 - cohort.study$Y)
  sampling.fraction.controls.strata2 = sum(cohort.study$Y)/sum(1 - cohort.study$Y)
  
  cohort.study$sampling.fraction = NA
  cohort.study$sampling.fraction[which(cohort.study$Y == 1)] = sampling.fraction.cases.strata1
  cohort.study$sampling.fraction[which(cohort.study$Y == 0)] = sampling.fraction.controls.strata2
  
  set.seed(i)
  R.cases.strata1 = rbinom(length(which(cohort.study$Y == 1 & cohort.study$V3 == 0)), 1, sampling.fraction.cases.strata1)
  set.seed(i)
  R.cases.strata2 = rbinom(length(which(cohort.study$Y == 1 & cohort.study$V3 == 1)), 1, sampling.fraction.cases.strata2)
  set.seed(i)
  R.controls.strata1 = rbinom(length(which(cohort.study$Y == 0 & cohort.study$V3 == 0)), 1, sampling.fraction.controls.strata1)
  set.seed(i)
  R.controls.strata2 = rbinom(length(which(cohort.study$Y == 0 & cohort.study$V3 == 1)), 1, sampling.fraction.controls.strata2)
  
  CC.cases.strata1 = cbind(cohort.study[which(cohort.study$Y == 1 & cohort.study$V3 == 0),], R=R.cases.strata1) 
  CC.cases.strata2 = cbind(cohort.study[which(cohort.study$Y == 1 & cohort.study$V3 == 1),], R=R.cases.strata2) 
  CC.controls.strata1 = cbind(cohort.study[which(cohort.study$Y == 0 & cohort.study$V3 == 0),], R=R.controls.strata1) 
  CC.controls.strata2 = cbind(cohort.study[which(cohort.study$Y == 0 & cohort.study$V3 == 1),], R=R.controls.strata2)
  CC.study = rbind(CC.cases.strata1,CC.cases.strata2, CC.controls.strata1, CC.controls.strata2)
  CC.study = CC.study[CC.study$R==1, ]
  
  data.phase.I = cohort.study
  data.phase.I.missreg = cbind(aggregate(data.phase.I$Y, by = list(data.phase.I$V3), FUN=sum), (aggregate(data.phase.I$Y, by = list(data.phase.I$V3), FUN=length)$x - aggregate(data.phase.I$Y, by = list(data.phase.I$V3), FUN=sum)$x))
  nn1 = data.phase.I.missreg[,2]
  nn0 = data.phase.I.missreg[,3]
  colnames(data.phase.I.missreg) = c("V3","case","control")
  data.phase.I.missreg$V4 = NA
  data.phase.I.missreg$V5 = NA
  #data.phase.I.missreg$V3 =NA
  data.phase.I.missreg$X1.cont =NA
  data.phase.I.missreg$obstype = "strata"
  data.phase.I.missreg = data.phase.I.missreg[c("V3","X1.cont","V4","V5","case","control","obstype")]
  
  data.phase.II.missreg = CC.study
  data.phase.II.missreg$case = data.phase.II.missreg$Y
  data.phase.II.missreg$control = 1 - data.phase.II.missreg$Y
  data.phase.II.missreg$obstype = "retro"
  data.phase.II.missreg = data.phase.II.missreg[c("V3","X1.cont","V4","V5","case","control","obstype")]
  
  data.missreg = rbind(data.phase.II.missreg, data.phase.I.missreg)
  fit.npmle <- bin2stg(cbind(case,control) ~ X1.cont + V4 + V5, xstrata=c("V3"), data=data.missreg, xs.includes=TRUE)
  if(fit.npmle$error == 0)
  {
    results.npmle.missreg[i, ] = c(fit.npmle$coefficients, diag(fit.npmle$cov), fit.npmle$error)
  } else {
    
    results.npmle.missreg[i,] = rep(NA, 9) 
  }
  
}
saveRDS(results.npmle.missreg, file="/Users/prosenjitkundu/Desktop/Code_biometrics_github/Simulation/Table1_results/missreg_case_control_table1.rds")
#saveRDS(results.npmle.missreg, file="missreg_GENMETA_without_measurement_error_simulation_1_case_control_Feb_18_2021.rds")
