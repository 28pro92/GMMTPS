library(OrdNor)
library(survival)
source("myoptim_realdata.R")

myexpit = function(x)
{
  return(as.numeric(1/(1 + exp(-x))))
}

expit <- function(beta, X)
{
  return(as.numeric(1/(1 + exp(-as.matrix(X) %*% beta))))
}

myexpit_d1 = function(x)
{
  return(as.numeric((1/(1 + exp(-x))) * (1/(1 + exp(x)))))
}
expit_d1 <- function(beta, X)
{
  return(as.numeric((1/(1 + exp(-as.matrix(X) %*% beta))) * (1/(1 + exp(as.matrix(X) %*% beta)))))
}



#---- Calculating Omega -----#

Omega = function(beta, study.list, phase_I_var, phase_II_var)
{
  
  cohort.study = study.list[[1]]
  nested.case.control.study = study.list[[2]]
  
  phase_I.indices = match(phase_I_var, colnames(study.list[[2]]))
  phase_II.indices = match(phase_II_var, colnames(study.list[[2]]))
  
  
  fit.R = glm(as.formula(paste("rel", "~", paste(phase_I_var[-1], collapse = "+"))), data = cohort.study, family = binomial())
  theta.R = fit.R$coefficients
  
  
  correction_factor_var_index = match(correction_factor_var, colnames(study.list[[2]]))
  
  beta_old_matrix = replicate(dim(study.list[[2]])[1], beta)
  beta_old_matrix[1, ] = beta_old_matrix[1,] + log(study.list[[2]][, correction_factor_var_index])
  theta_CC = t(beta_old_matrix)
  
  no.of.cases.CC.inst1 = length(which(nested.case.control.study$rel == 1 & nested.case.control.study$instit2 == 0))
  no.of.cases.CC.inst2 = length(which(nested.case.control.study$rel == 1 & nested.case.control.study$instit2 == 1))
  no.of.controls.CC.inst1 = length(which(nested.case.control.study$rel == 0 & nested.case.control.study$instit2 == 0))
  no.of.controls.CC.inst2 = length(which(nested.case.control.study$rel == 0 & nested.case.control.study$instit2 == 1))
  
  no.of.cases.R.inst1 = length(which(cohort.study$rel == 1 & cohort.study$instit2 == 0))
  no.of.cases.R.inst2 = length(which(cohort.study$rel == 1 & cohort.study$instit2 == 1))
  no.of.controls.R.inst1 = length(which(cohort.study$rel == 0 & cohort.study$instit2 == 0))
  no.of.controls.R.inst2 = length(which(cohort.study$rel == 0 & cohort.study$instit2 == 1))
  
  cases.CC.inst1 = nested.case.control.study[which(nested.case.control.study$rel == 1 & nested.case.control.study$instit2 == 0), ]
  cases.CC.inst2 = nested.case.control.study[which(nested.case.control.study$rel == 1 & nested.case.control.study$instit2 == 1), ]
  controls.CC.inst1 = nested.case.control.study[which(nested.case.control.study$rel == 0 & nested.case.control.study$instit2 == 0), ]
  controls.CC.inst2 = nested.case.control.study[which(nested.case.control.study$rel == 0 & nested.case.control.study$instit2 == 1),]
  
  theta.cases.CC.inst1 = theta_CC[which(nested.case.control.study$rel == 1 & nested.case.control.study$instit2 == 0), ]
  theta.cases.CC.inst2 = theta_CC[which(nested.case.control.study$rel == 1 & nested.case.control.study$instit2 == 1), ]
  theta.controls.CC.inst1 = theta_CC[which(nested.case.control.study$rel == 0 & nested.case.control.study$instit2 == 0), ]
  theta.controls.CC.inst2 = theta_CC[which(nested.case.control.study$rel == 0 & nested.case.control.study$instit2 == 1), ]
  
  Y.nested.case.control.study = nested.case.control.study$rel
  Y.cases.CC.inst1 = cases.CC.inst1$rel
  Y.cases.CC.inst2 = cases.CC.inst2$rel
  Y.controls.CC.inst1 = controls.CC.inst1$rel
  Y.controls.CC.inst2 = controls.CC.inst2$rel
  
  cases.CC.inst1 = as.matrix(cases.CC.inst1)
  cases.CC.inst2 = as.matrix(cases.CC.inst2)
  controls.CC.inst1 = as.matrix(controls.CC.inst1)
  controls.CC.inst2 = as.matrix(controls.CC.inst2)
  nested.case.control.study = as.matrix(nested.case.control.study)
  
  
  
  Omega_11.cases.inst1 = (t(cases.CC.inst1[, phase_I.indices ]) %*% diag((expit(beta, cases.CC.inst1[, phase_II.indices]) - expit(theta.R, cases.CC.inst1[, phase_I.indices ]))^2*((no.of.cases.R.inst1/no.of.cases.CC.inst1)/cases.CC.inst1[,9])) %*% cases.CC.inst1[, phase_I.indices ])/study.list[[3]]
  Omega_11.cases.inst2 = (t(cases.CC.inst2[, phase_I.indices ]) %*% diag((expit(beta, cases.CC.inst2[, phase_II.indices]) - expit(theta.R, cases.CC.inst2[, phase_I.indices ]))^2*((no.of.cases.R.inst2/no.of.cases.CC.inst2)/cases.CC.inst2[,9])) %*% cases.CC.inst2[, phase_I.indices ])/study.list[[3]]
  Omega_11.controls.inst1 = (t(controls.CC.inst1[,phase_I.indices ]) %*% diag((expit(beta, controls.CC.inst1[, phase_II.indices]) - expit(theta.R, controls.CC.inst1[, phase_I.indices ]))^2*((no.of.controls.R.inst1/no.of.controls.CC.inst1)/controls.CC.inst1[,9])) %*% controls.CC.inst1[, phase_I.indices])/study.list[[3]]
  Omega_11.controls.inst2 = (t(controls.CC.inst2[,phase_I.indices ]) %*% diag((expit(beta, controls.CC.inst2[, phase_II.indices]) - expit(theta.R, controls.CC.inst2[, phase_I.indices ]))^2*((no.of.controls.R.inst2/no.of.controls.CC.inst2)/controls.CC.inst2[,9])) %*% controls.CC.inst2[, phase_I.indices ])/study.list[[3]]
  Omega_11 = Omega_11.cases.inst1 + Omega_11.cases.inst2 + Omega_11.controls.inst1 + Omega_11.controls.inst2
  
  
  #Omega_22 = (t(nested.case.control.study[, 3:6]) %*% diag((Y.nested.case.control.study - expit(theta.CC, nested.case.control.study[, 3:6]))^2) %*% nested.case.control.study[, 3:6])/study.list[[3]]
  Omega_22.cases.inst1 = (t(cases.CC.inst1[,  phase_II.indices]) %*% diag((Y.cases.CC.inst1 - myexpit(rowSums(theta.cases.CC.inst1 * cases.CC.inst1[, phase_II.indices])))^2 * cases.CC.inst1[,9]) %*% cases.CC.inst1[,  phase_II.indices])  * ((no.of.cases.R.inst1/study.list[[3]])/no.of.cases.CC.inst1)
  Omega_22.cases.inst2 = (t(cases.CC.inst2[,  phase_II.indices]) %*% diag((Y.cases.CC.inst2 - myexpit(rowSums(theta.cases.CC.inst2 * cases.CC.inst2[, phase_II.indices])))^2  * cases.CC.inst2[,9]) %*% cases.CC.inst2[,  phase_II.indices]) * ((no.of.cases.R.inst2/study.list[[3]])/no.of.cases.CC.inst2)
  Omega_22.controls.inst1 = (t(controls.CC.inst1[,  phase_II.indices]) %*% diag((Y.controls.CC.inst1 - myexpit(rowSums(theta.controls.CC.inst1 * controls.CC.inst1[, phase_II.indices])))^2 * controls.CC.inst1[,9]) %*% controls.CC.inst1[,  phase_II.indices])  * ((no.of.controls.R.inst1/study.list[[3]])/no.of.controls.CC.inst1)
  Omega_22.controls.inst2 = (t(controls.CC.inst2[,  phase_II.indices]) %*% diag((Y.controls.CC.inst2 - myexpit(rowSums(theta.controls.CC.inst2 * controls.CC.inst2[, phase_II.indices])))^2 * controls.CC.inst2[,9]) %*% controls.CC.inst2[, phase_II.indices])  * ((no.of.controls.R.inst2/study.list[[3]])/no.of.controls.CC.inst2)
  Omega_22 = Omega_22.cases.inst1 + Omega_22.cases.inst2 + Omega_22.controls.inst1 + Omega_22.controls.inst2
  
  Omega_33.cases.inst1 = t(cases.CC.inst1[, phase_I.indices ]) %*% diag((Y.cases.CC.inst1 - expit(theta.R, cases.CC.inst1[, phase_I.indices ]))^2) %*% cases.CC.inst1[, phase_I.indices ] * ((no.of.cases.R.inst1/study.list[[3]])/no.of.cases.CC.inst1)
  Omega_33.cases.inst2 = t(cases.CC.inst2[, phase_I.indices ]) %*% diag((Y.cases.CC.inst2 - expit(theta.R, cases.CC.inst2[, phase_I.indices ]))^2) %*% cases.CC.inst2[, phase_I.indices ] * ((no.of.cases.R.inst2/study.list[[3]])/no.of.cases.CC.inst2)
  Omega_33.controls.inst1 = t(controls.CC.inst1[, phase_I.indices ]) %*% diag((Y.controls.CC.inst1 - expit(theta.R, controls.CC.inst1[, phase_I.indices ]))^2) %*% controls.CC.inst1[, phase_I.indices ] * ((no.of.controls.R.inst1/study.list[[3]])/no.of.controls.CC.inst1)
  Omega_33.controls.inst2 = t(controls.CC.inst2[, phase_I.indices ]) %*% diag((Y.controls.CC.inst2 - expit(theta.R, controls.CC.inst2[, phase_I.indices ]))^2) %*% controls.CC.inst2[, phase_I.indices ] * ((no.of.controls.R.inst2/study.list[[3]])/no.of.controls.CC.inst2)
  Omega_33 = Omega_33.cases.inst1 + Omega_33.cases.inst2 + Omega_33.controls.inst1 + Omega_33.controls.inst2
  
  Omega_12.cases.inst1 = t(cases.CC.inst1[, phase_I.indices]) %*% diag((expit(beta, cases.CC.inst1[,  phase_II.indices]) - expit(theta.R, cases.CC.inst1[, phase_I.indices ]))*(Y.cases.CC.inst1 - myexpit(rowSums(theta.cases.CC.inst1 * cases.CC.inst1[, phase_II.indices])))) %*% cases.CC.inst1[,  phase_II.indices] * ((no.of.cases.R.inst1/study.list[[3]])/no.of.cases.CC.inst1)
  Omega_12.cases.inst2 = t(cases.CC.inst2[, phase_I.indices]) %*% diag((expit(beta, cases.CC.inst2[,  phase_II.indices]) - expit(theta.R, cases.CC.inst2[, phase_I.indices ]))*(Y.cases.CC.inst2 - myexpit(rowSums(theta.cases.CC.inst2 * cases.CC.inst2[, phase_II.indices])))) %*% cases.CC.inst2[,  phase_II.indices] * ((no.of.cases.R.inst2/study.list[[3]])/no.of.cases.CC.inst2)
  Omega_12.controls.inst1 = t(controls.CC.inst1[, phase_I.indices ]) %*% diag((expit(beta, controls.CC.inst1[,  phase_II.indices]) - expit(theta.R, controls.CC.inst1[,phase_I.indices ]))*(Y.controls.CC.inst1 - myexpit(rowSums(theta.controls.CC.inst1 * controls.CC.inst1[, phase_II.indices])))) %*% controls.CC.inst1[,  phase_II.indices] * ((no.of.controls.R.inst1/study.list[[3]])/no.of.controls.CC.inst1)
  Omega_12.controls.inst2= t(controls.CC.inst2[, phase_I.indices ]) %*% diag((expit(beta, controls.CC.inst2[,  phase_II.indices]) - expit(theta.R, controls.CC.inst2[,phase_I.indices ]))*(Y.controls.CC.inst2 - myexpit(rowSums(theta.controls.CC.inst2 * controls.CC.inst2[, phase_II.indices])))) %*% controls.CC.inst2[,  phase_II.indices] * ((no.of.controls.R.inst2/study.list[[3]])/no.of.controls.CC.inst2)
  Omega_12 = Omega_12.cases.inst1 + Omega_12.cases.inst2 + Omega_12.controls.inst1 + Omega_12.controls.inst2
  
  #Omega_23 = (t(nested.case.control.study[, 3:6]) %*% diag((Y.nested.case.control.study - expit(theta.CC, nested.case.control.study[, 3:6]))*(Y.nested.case.control.study - expit(theta.R, nested.case.control.study[, 3:4]))) %*% nested.case.control.study[, 3:4])/study.list[[3]]
  Omega_23.cases.inst1 = (t(cases.CC.inst1[,  phase_II.indices]) %*% diag((Y.cases.CC.inst1 - myexpit(rowSums(theta.cases.CC.inst1 * cases.CC.inst1[, phase_II.indices])))*(Y.cases.CC.inst1 - expit(theta.R, cases.CC.inst1[, phase_I.indices ])) * cases.CC.inst1[,9]) %*% cases.CC.inst1[, phase_I.indices ]) * ((no.of.cases.R.inst1/study.list[[3]])/no.of.cases.CC.inst1)
  Omega_23.cases.inst2 = (t(cases.CC.inst2[,  phase_II.indices]) %*% diag((Y.cases.CC.inst2 - myexpit(rowSums(theta.cases.CC.inst2 * cases.CC.inst2[, phase_II.indices])))*(Y.cases.CC.inst2 - expit(theta.R, cases.CC.inst2[, phase_I.indices ])) * cases.CC.inst2[,9]) %*% cases.CC.inst2[, phase_I.indices ]) * ((no.of.cases.R.inst2/study.list[[3]])/no.of.cases.CC.inst2)
  Omega_23.controls.inst1 = (t(controls.CC.inst1[,  phase_II.indices]) %*% diag((Y.controls.CC.inst1 - myexpit(rowSums(theta.controls.CC.inst1 * controls.CC.inst1[, phase_II.indices])))*(Y.controls.CC.inst1 - expit(theta.R, controls.CC.inst1[, phase_I.indices ])) * controls.CC.inst1[,9]) %*% controls.CC.inst1[, phase_I.indices ]) * ((no.of.controls.R.inst1/study.list[[3]])/no.of.controls.CC.inst1)
  Omega_23.controls.inst2 = (t(controls.CC.inst2[,  phase_II.indices]) %*% diag((Y.controls.CC.inst2 - myexpit(rowSums(theta.controls.CC.inst2 * controls.CC.inst2[, phase_II.indices])))*(Y.controls.CC.inst2 - expit(theta.R, controls.CC.inst2[, phase_I.indices ]))* controls.CC.inst2[,9]) %*% controls.CC.inst2[, phase_I.indices ]) * ((no.of.controls.R.inst2/study.list[[3]])/no.of.controls.CC.inst2)
  Omega_23 = Omega_23.cases.inst1 + Omega_23.cases.inst2 + Omega_23.controls.inst1 + Omega_23.controls.inst2
  
  Omega_13.cases.inst1 = t(cases.CC.inst1[, phase_I.indices ]) %*% diag((expit(beta, cases.CC.inst1[,  phase_II.indices]) - expit(theta.R, cases.CC.inst1[, phase_I.indices ]))*(Y.cases.CC.inst1 - expit(theta.R, cases.CC.inst1[, phase_I.indices ]))) %*% cases.CC.inst1[, phase_I.indices ] * ((no.of.cases.R.inst1/study.list[[3]])/no.of.cases.CC.inst1)
  Omega_13.cases.inst2 = t(cases.CC.inst2[, phase_I.indices ]) %*% diag((expit(beta, cases.CC.inst2[,  phase_II.indices]) - expit(theta.R, cases.CC.inst2[, phase_I.indices ]))*(Y.cases.CC.inst2 - expit(theta.R, cases.CC.inst2[, phase_I.indices ]))) %*% cases.CC.inst2[, phase_I.indices ] * ((no.of.cases.R.inst2/study.list[[3]])/no.of.cases.CC.inst2)
  Omega_13.controls.inst1 = t(controls.CC.inst1[, phase_I.indices ]) %*% diag((expit(beta, controls.CC.inst1[,  phase_II.indices]) - expit(theta.R, controls.CC.inst1[, phase_I.indices ]))*(Y.controls.CC.inst1 - expit(theta.R, controls.CC.inst1[, phase_I.indices ]))) %*% controls.CC.inst1[, phase_I.indices ] * ((no.of.controls.R.inst1/study.list[[3]])/no.of.controls.CC.inst1)
  Omega_13.controls.inst2 = t(controls.CC.inst2[, phase_I.indices ]) %*% diag((expit(beta, controls.CC.inst2[,  phase_II.indices]) - expit(theta.R, controls.CC.inst2[, phase_I.indices ]))*(Y.controls.CC.inst2 - expit(theta.R, controls.CC.inst2[, phase_I.indices ]))) %*% controls.CC.inst2[, phase_I.indices ] * ((no.of.controls.R.inst2/study.list[[3]])/no.of.controls.CC.inst2)
  Omega_13 = Omega_13.cases.inst1 + Omega_13.cases.inst2 + Omega_13.controls.inst1 + Omega_13.controls.inst2
  
  Omega.est = rbind(cbind(Omega_11, Omega_12, Omega_13), cbind(t(Omega_12), Omega_22, Omega_23), cbind(t(Omega_13), t(Omega_23), Omega_33))
  
  Omega.est
  
  
  
}




#---- Calculating Delta -----#
Delta = function(beta, study.list, phase_I_var, phase_II_var,lambda.in.Delta)
{
  
  #--Change here when you increase number of risk factors in both phases---#
  no.of.covariates.inc.intercept.phase.I = length(phase_I_var)
  no.of.covariates.inc.intercept.phase.II = length(phase_II_var)
  
  
  lambda.hat.matrix = (1/lambda.in.Delta) * diag(no.of.covariates.inc.intercept.phase.II)
  Identity_I = diag(no.of.covariates.inc.intercept.phase.I)
  Identity_II = diag(no.of.covariates.inc.intercept.phase.II)
  Zero_I_II = matrix(0, no.of.covariates.inc.intercept.phase.II, no.of.covariates.inc.intercept.phase.I)
  
  #V1I = Robust.info.mtrx(beta, study.list)/study.list[[3]]
  
  Delta.est = cbind(adiag(Identity_I, lambda.hat.matrix), rbind(-Identity_I, Zero_I_II))
  
  return(Delta.est)
}




#---- Calculating Gamma ---#


Gamma = function(beta, study.list, phase_I_var, phase_II_var)
{
  cohort.study = study.list[[1]]
  nested.case.control.study = study.list[[2]]
  
  phase_I.indices = match(phase_I_var, colnames(study.list[[2]]))
  phase_II.indices = match(phase_II_var, colnames(study.list[[2]]))
  
  #fit.R = glm(cohort.study$rel ~ cohort.study$V3 + cohort.study$stage2 + cohort.study$stage3 + cohort.study$stage4 + cohort.study$age, family = binomial())
  fit.R = glm(as.formula(paste("rel", "~", paste(phase_I_var[-1], collapse = "+"))), data = cohort.study, family = binomial())
  theta.R = fit.R$coefficients
  #status.converged = fit.R$converged
  
  correction_factor_var_index = match(correction_factor_var, colnames(study.list[[2]]))
  
  beta_old_matrix = replicate(dim(study.list[[2]])[1], beta)
  beta_old_matrix[1, ] = beta_old_matrix[1,] + log(study.list[[2]][, correction_factor_var_index])
  theta_CC = t(beta_old_matrix)
  
  no.of.cases.CC.inst1 = length(which(nested.case.control.study$rel == 1 & nested.case.control.study$instit2 == 0))
  no.of.cases.CC.inst2 = length(which(nested.case.control.study$rel == 1 & nested.case.control.study$instit2 == 1))
  no.of.controls.CC.inst1 = length(which(nested.case.control.study$rel == 0 & nested.case.control.study$instit2 == 0))
  no.of.controls.CC.inst2 = length(which(nested.case.control.study$rel == 0 & nested.case.control.study$instit2 == 1))
  
  no.of.cases.R.inst1 = length(which(cohort.study$rel == 1 & cohort.study$instit2 == 0))
  no.of.cases.R.inst2 = length(which(cohort.study$rel == 1 & cohort.study$instit2 == 1))
  no.of.controls.R.inst1 = length(which(cohort.study$rel == 0 & cohort.study$instit2 == 0))
  no.of.controls.R.inst2 = length(which(cohort.study$rel == 0 & cohort.study$instit2 == 1))
  
  cases.CC.inst1 = nested.case.control.study[which(nested.case.control.study$rel == 1 & nested.case.control.study$instit2 == 0), ]
  cases.CC.inst2 = nested.case.control.study[which(nested.case.control.study$rel == 1 & nested.case.control.study$instit2 == 1), ]
  controls.CC.inst1 = nested.case.control.study[which(nested.case.control.study$rel == 0 & nested.case.control.study$instit2 == 0), ]
  controls.CC.inst2 = nested.case.control.study[which(nested.case.control.study$rel == 0 & nested.case.control.study$instit2 == 1),]
  
  theta.cases.CC.inst1 = theta_CC[which(nested.case.control.study$rel == 1 & nested.case.control.study$instit2 == 0), ]
  theta.cases.CC.inst2 = theta_CC[which(nested.case.control.study$rel == 1 & nested.case.control.study$instit2 == 1), ]
  theta.controls.CC.inst1 = theta_CC[which(nested.case.control.study$rel == 0 & nested.case.control.study$instit2 == 0), ]
  theta.controls.CC.inst2 = theta_CC[which(nested.case.control.study$rel == 0 & nested.case.control.study$instit2 == 1), ]
  
  
  
  cases.CC.inst1 = as.matrix(cases.CC.inst1)
  cases.CC.inst2 = as.matrix(cases.CC.inst2)
  controls.CC.inst1 = as.matrix(controls.CC.inst1)
  controls.CC.inst2 = as.matrix(controls.CC.inst2)
  nested.case.control.study = as.matrix(nested.case.control.study)
  
  #print(head(cases.CC[, 3:6]))
  #print(class(cases.CC[, 3:6]))
  #print(expit_d1(beta, cases.CC[, 3:6]))
  #print(((no.of.cases.R/study.list[[3]])/no.of.cases.CC))
  
  
  Gamma.1.cases.inst1 = (t(cases.CC.inst1[, phase_I.indices]) %*% diag(expit_d1(beta, cases.CC.inst1[,  phase_II.indices])) %*% cases.CC.inst1[,  phase_II.indices]) * ((no.of.cases.R.inst1/study.list[[3]])/no.of.cases.CC.inst1)
  Gamma.1.cases.inst2 = (t(cases.CC.inst2[, phase_I.indices]) %*% diag(expit_d1(beta, cases.CC.inst2[,  phase_II.indices])) %*% cases.CC.inst2[,  phase_II.indices]) * ((no.of.cases.R.inst2/study.list[[3]])/no.of.cases.CC.inst2)
  Gamma.1.controls.inst1 = (t(controls.CC.inst1[, phase_I.indices]) %*% diag(expit_d1(beta, controls.CC.inst1[,  phase_II.indices])) %*% controls.CC.inst1[,  phase_II.indices]) * ((no.of.controls.R.inst1/study.list[[3]])/no.of.controls.CC.inst1)
  Gamma.1.controls.inst2 = (t(controls.CC.inst2[, phase_I.indices]) %*% diag(expit_d1(beta, controls.CC.inst2[,  phase_II.indices])) %*% controls.CC.inst2[,  phase_II.indices]) * ((no.of.controls.R.inst2/study.list[[3]])/no.of.controls.CC.inst2)
  Gamma.1 = Gamma.1.cases.inst1 + Gamma.1.cases.inst2 + Gamma.1.controls.inst1 + Gamma.1.controls.inst2
  
  #Gamma.2 = -(t(nested.case.control.study[, 3:6]) %*% diag(expit_d1(theta.CC, nested.case.control.study[, 3:6])) %*% nested.case.control.study[, 3:6])/study.list[[4]]
  
  
  
  Gamma.2.cases.inst1 = (t(cases.CC.inst1[,  phase_II.indices]) %*% diag(myexpit_d1(rowSums(theta.cases.CC.inst1 * cases.CC.inst1[, phase_II.indices])) * cases.CC.inst1[,9]) %*% cases.CC.inst1[,  phase_II.indices]) * ((no.of.cases.R.inst1/study.list[[3]])/no.of.cases.CC.inst1)
  Gamma.2.cases.inst2 = (t(cases.CC.inst2[,  phase_II.indices]) %*% diag(myexpit_d1(rowSums(theta.cases.CC.inst2 * cases.CC.inst2[, phase_II.indices])) * cases.CC.inst2[,9]) %*% cases.CC.inst2[,  phase_II.indices]) * ((no.of.cases.R.inst2/study.list[[3]])/no.of.cases.CC.inst2)
  Gamma.2.controls.inst1 = (t(controls.CC.inst1[,  phase_II.indices]) %*% diag(myexpit_d1(rowSums(theta.controls.CC.inst1 * controls.CC.inst1[, phase_II.indices])) * controls.CC.inst1[,9]) %*% controls.CC.inst1[,  phase_II.indices]) * ((no.of.controls.R.inst1/study.list[[3]])/no.of.controls.CC.inst1)
  Gamma.2.controls.inst2 = (t(controls.CC.inst2[, phase_II.indices]) %*% diag(myexpit_d1(rowSums(theta.controls.CC.inst2 * controls.CC.inst2[, phase_II.indices])) * controls.CC.inst2[,9]) %*% controls.CC.inst2[, phase_II.indices]) * ((no.of.controls.R.inst2/study.list[[3]])/no.of.controls.CC.inst2)
  Gamma.2 = -(1/lambda.in.Delta)*(Gamma.2.cases.inst1 + Gamma.2.cases.inst2 + Gamma.2.controls.inst1 + Gamma.2.controls.inst2)
  
  
  Gamma.est = rbind(Gamma.1, Gamma.2)
  
  Gamma.est
  
}


#----------Calculating C_opt-----#
C_opt <- function(beta, study.list, phase_I_var, phase_II_var, lambda.in.Delta)
{
  if(det(Delta(beta, study.list, phase_I_var, phase_II_var, lambda.in.Delta) %*% Omega(beta, study.list, phase_I_var, phase_II_var) %*% t(Delta(beta, study.list, phase_I_var, phase_II_var, lambda.in.Delta ))) == 0)
  {
    return(NA)
  }else{
    return(solve(Delta(beta, study.list, phase_I_var, phase_II_var, lambda.in.Delta) %*% Omega(beta, study.list, phase_I_var, phase_II_var) %*% t(Delta(beta, study.list, phase_I_var, phase_II_var, lambda.in.Delta)), tol=1e-60))
  }
  
}

instit_margin = as.numeric(prop.table(margin.table(table(nwtco$instit, nwtco$histol, nwtco$stage), 1)))
hist_margin = as.numeric(prop.table(margin.table(table(nwtco$instit, nwtco$histol, nwtco$stage), 2)))
stage_margin = as.numeric(prop.table(margin.table(table(nwtco$instit, nwtco$histol, nwtco$stage), 3)))
#set.seed(1)
# Sets the marginals.
# The values are cumulative so for the first variable the first marginal will be .1, the second is .2, the third is .3, and the fourth is .4
marginal = list(cumsum(instit_margin)[-2], cumsum(hist_margin)[-2],cumsum(stage_margin)[-4])
# Checks the lower and upper bounds of the correlation coefficients.
corrcheck(marginal)

beta_real_data = as.numeric(glm(rel~factor(histol)*factor(stage) + factor(histol)*scale(age/12) + factor(stage)*scale(age/12) , data = nwtco, family = binomial())$coefficients)
beta_real_data[1] = -3.6
#beta_real_data[2] = beta_real_data[3]
#beta_real_data = beta_real_data[-11]
#beta_real_data = beta_real_data[-2]
corr = list()
# rho_inst_histol =  polycor::polychor(as.factor(nwtco$instit), as.factor(nwtco$histol))
# rho_inst_stage =  polycor::polychor(as.factor(nwtco$instit), as.factor(nwtco$stage))
# rho_inst_age = polycor::polyserial(nwtco$age, as.factor(nwtco$instit))
# rho_histol_stage = polycor::polychor(as.factor(nwtco$histol), as.factor(nwtco$stage))
# rho_histol_age = polycor::polyserial(nwtco$age, as.factor(nwtco$histol))
# rho_stage_age = polycor::polyserial(nwtco$age, as.factor(nwtco$stage))
# rho_nwtco_matrix = matrix(c(1,rho_inst_histol,rho_inst_stage,rho_inst_age,rho_inst_histol,1,rho_histol_stage,rho_histol_age, rho_inst_stage, rho_histol_stage, 1, rho_stage_age, rho_inst_age, rho_histol_age, rho_stage_age,1), nrow = 4, ncol = 4)
# 
# rho_nwtco_matrix[1,2] = 0.92
# rho_nwtco_matrix[2,1] = 0.92

rho_nwtco_matrix = cor(cbind(nwtco$instit, nwtco$histol, nwtco$stage,scale(nwtco$age/12)))
new_corr = rbind(cbind(rho_nwtco_matrix, c(0,0,0,0)), c(0,0,0,0,1))

#rownames(rho_nwtco) = colnames(rho_nwtco) = c("O1","O2","O3", "C1")
mean.Y = c()
sim.results = matrix(NA, 1000, 27)
for(i in 1:1000)
{
  #data.list = corrvar(n = 10000, k_cat = 3, k_cont = 1, means = 0, vars = 1, method = "Polynomial",skews = 0, skurts = 0, fifths = 0, sixths = 0,
  #                marginal = marginal, rho = rho_nwtco, quiet = TRUE, seed=i, errorloop = TRUE, epsilon = 0.01)
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
  data.model.matrix = model.matrix(~ instit * histol * stage * age, data = data)
  true.model.variables = c("(Intercept)", "histol2","stage2","stage3","stage4","age","histol2:stage2", "histol2:stage3", "histol2:stage4", "histol2:age","stage2:age", "stage3:age", "stage4:age")
  data.model.matrix = data.model.matrix[, match(true.model.variables,colnames(data.model.matrix))]
  set.seed(i)
  outcome = rbinom(10000, 1, expit(beta_real_data, data.model.matrix))
  data_sim = cbind(outcome, data)
  colnames(data_sim)[1] = "rel"
  cohort.study = as.data.frame(data_sim)
  #mean.Y = c(mean.Y, mean(cohort.study$rel))
  #cohort.study$histol <- factor(cohort.study$histol)
  #cohort.study$instit <- factor(cohort.study$instit)
  #cohort.study$stage <- factor(cohort.study$stage)
  #cohort.study$age  = as.numeric(scale(cohort.study$age/12))
  #cohort.strata.table <- table(cohort.study$rel, cohort.study$instit)
  
  sampling.fraction.cases.strata1 = 1
  sampling.fraction.cases.strata2 = 1
  sampling.fraction.controls.strata1 = sum(cohort.study$rel == 1 & cohort.study$instit == 2)/sum(cohort.study$rel == 0 & cohort.study$instit == 1)
  sampling.fraction.controls.strata2 = sum(cohort.study$rel == 1 & cohort.study$instit == 1)/sum(cohort.study$rel == 0 & cohort.study$instit == 2)
  
  
  cohort.study$sampling.fraction = NA
  cohort.study$sampling.fraction[which(cohort.study$rel == 1 & cohort.study$instit == 1)] = sampling.fraction.cases.strata1
  cohort.study$sampling.fraction[which(cohort.study$rel == 1 & cohort.study$instit == 2)] = sampling.fraction.cases.strata2
  cohort.study$sampling.fraction[which(cohort.study$rel == 0 & cohort.study$instit == 1)] = sampling.fraction.controls.strata1
  cohort.study$sampling.fraction[which(cohort.study$rel == 0 & cohort.study$instit == 2)] = sampling.fraction.controls.strata2
  
  divison = function(a,b){b[a==1]/b[a==0]}
  start.time = Sys.time()
  cohort.study.design.matrix <- data.frame(model.matrix(~ rel + instit + histol + stage + age + instit*stage*age + histol*stage*age, data = cohort.study))
  #phase_I_var = c("X.Intercept.","instit2", "stage2", "stage3", "stage4", "age")
  phase_I_var = c("X.Intercept.","instit2", "stage2", "stage3", "stage4", "age","instit2.stage2", "instit2.stage3", "instit2.stage4", "instit2.age", "stage2.age", "stage3.age", "stage4.age")
  phase_II_var = c("X.Intercept.", "histol2","stage2","stage3","stage4","age","histol2.stage2", "histol2.stage3", "histol2.stage4", "histol2.age", "stage2.age", "stage3.age", "stage4.age")
  fit.full.cohort = glm(as.formula(paste("rel", "~", paste(phase_II_var[-1], collapse = "+"))), family = binomial(), data = cohort.study.design.matrix)
  Q.beta.initial.cohort = as.numeric(fit.full.cohort$coefficients)
  #Q.beta.initial.cohort =beta_real_data
  
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
  
  
  #CC.study = CC.study %>% mutate(correction.factor = CC.study[CC.study$rel == 1, ]$sampling.fraction[1]/CC.study[CC.study$rel == 0, ]$sampling.fraction[1])
  
  correction_factor_var = "correction.factor"
  
  #CC.study <- CC.study[order(CC.study$seqno), ]
  CC.size <- dim(CC.study)[1]
  cohort.size <- dim(cohort.study)[1]
  
  disease_var = "rel"
  samp_fraction_var = "sampling.fraction"
  sampling_fraction_est_var = "sampling.fraction.est"
  strata_var = "instit"
  
  table.disease.strata.samp.frac = aggregate(as.formula(paste(samp_fraction_var, "~", paste(c(disease_var, strata_var), collapse = "+"))), data=CC.study, FUN="unique")
  strata_by_samp_frac = table.disease.strata.samp.frac %>% group_by(instit) %>% summarise(correction.factor = divison(rel,sampling.fraction))
  
  CC.study = full_join(CC.study, strata_by_samp_frac, by = strata_var)
  
  CC.study$NDS = NA
  CC.study$nDS = NA
  CC.study$NDS[which(CC.study$rel == 1 & CC.study$instit == 1)] = length(R.cases.strata1)
  CC.study$NDS[which(CC.study$rel == 1 & CC.study$instit == 2)] = length(R.cases.strata2)
  CC.study$NDS[which(CC.study$rel == 0 & CC.study$instit == 1)] = length(R.controls.strata1)
  CC.study$NDS[which(CC.study$rel == 0 & CC.study$instit == 2)] = length(R.controls.strata2)
  CC.study$nDS[which(CC.study$rel == 1 & CC.study$instit == 1)] = sum(R.cases.strata1)
  CC.study$nDS[which(CC.study$rel == 1 & CC.study$instit == 2)] = sum(R.cases.strata2)
  CC.study$nDS[which(CC.study$rel == 0 & CC.study$instit == 1)] = sum(R.controls.strata1)
  CC.study$nDS[which(CC.study$rel == 0 & CC.study$instit == 2)] = sum(R.controls.strata2)
  CC.study$sampling.fraction.est = CC.study$nDS/CC.study$NDS
  
  
  
  CC.study.design.matrix <- data.frame(model.matrix(~ rel + instit + histol + stage + age + instit*stage*age + histol*stage*age + sampling.fraction + correction.factor + sampling.fraction.est, data = CC.study))
  #cohort.study <- cbind(cohort.study, ones.cohort.study)
  #CC.study <- cbind(CC.study, ones.CC.study)
  
  
  
  Q.study.list <- list(cohort.study.design.matrix, CC.study.design.matrix, cohort.size, CC.size,  cohort.study, CC.study)
  
  lambda.in.Delta = ((sampling.fraction.cases.strata1 * length(R.cases.strata1)) + (sampling.fraction.cases.strata2 * length(R.cases.strata2)) + (sampling.fraction.controls.strata1 * length(R.controls.strata1))  + (sampling.fraction.controls.strata2 * length(R.controls.strata2)))/Q.study.list[[3]] 
  
  
  
  #Q.beta.initial.cohort = append(Q.beta.initial.cohort, 0, after = 1)
  #Q.beta.initial.CC = as.numeric(fit.CC$coefficients)
  #Q.beta.initial.CC = append(Q.beta.initial.CC, 0, after = 1)
  fit.R = glm(as.formula(paste("rel", "~", paste(phase_I_var[-1], collapse = "+"))), data = cohort.study.design.matrix, family = binomial())
  theta.R = fit.R$coefficients
  
  sampling_fraction_var = c("sampling.fraction.est")
  no_of_iter_outer = 0
  eps = 1e-06
  outer_iter = 0
  C = C_opt(Q.beta.initial.cohort, Q.study.list, phase_I_var, phase_II_var,lambda.in.Delta)
  #C = diag(length(phase_I_var) + length(phase_II_var))
  res.sim = myoptim(Q.study.list, C, Q.beta.initial.cohort, eps, no_of_iter_outer, phase_I_var, phase_II_var, theta.R, correction_factor_var, sampling_fraction_est_var)
  if(res.sim$Status == 0)
  {
    sim.results[i, ] = rep(NA, 27)
  }else{
    beta.initial = as.numeric(res.sim$beta_optim)
    C = C_opt(as.numeric(res.sim$beta_optim), Q.study.list, phase_I_var, phase_II_var,lambda.in.Delta)
    repeat{
      res.sim = myoptim(Q.study.list, C, beta.initial, eps, no_of_iter_outer, phase_I_var, phase_II_var, theta.R, correction_factor_var, sampling_fraction_est_var)
      #res.sim = optimx(Q.beta.initial.cohort, Q, study.list = Q.study.list, C = C, theta_R = theta.R, phase_I_var=phase_I_var, phase_II_var=phase_II_var, correction_factor_var = correction_factor_var, sampling_fraction_est_var = sampling_fraction_est_var, method = "BFGS")
      #beta.next <- as.numeric(res.sim)[1:4]
      beta.next = as.numeric(res.sim$beta_optim)
      
      status.convergence <- res.sim$Status
      if(status.convergence == 0)
      {
        break;
      }
      if(sqrt(sum(beta.next - beta.initial)^2) <= eps)
      {
        break;
      }
      
      beta.initial <- beta.next
      
      #print(C)
    }
    if(status.convergence == 1)
    {
      C = C_opt(beta.next, Q.study.list, phase_I_var, phase_II_var,lambda.in.Delta)
      #C = C_opt(beta.next, Q.study.list, lambda.in.Delta)
      #Var.est.formula.1 = solve((t(Gamma(beta.next, Q.study.list)) %*% C %*% Gamma(beta.next, Q.study.list)), tol=1e-30) %*% t(Gamma(beta.next, Q.study.list)) %*% C %*% Delta(Q.study.list)
      #Var.est.formula = (Var.est.formula.1 %*% Omega(beta.next, Q.study.list) %*% t(Var.est.formula.1))/Q.study.list[[4]]
      #sim.results[i,] = c(as.numeric(beta.next), as.numeric(diag(Var.est.formula)), as.numeric(Q.beta.initial.CC), as.numeric(diag(vcov(fit.CC))), status.convergence)
      #Var.est.formula = solve(t(Gamma(beta.next, Q.study.list)) %*% C %*% Gamma(beta.next, Q.study.list), tol = 1e-60 )/Q.study.list[[3]]
      Var.est.formula = solve(t(Gamma(beta.next, Q.study.list, phase_I_var, phase_II_var)) %*% C %*% Gamma(beta.next, Q.study.list, phase_I_var, phase_II_var), tol = 1e-60 )/Q.study.list[[3]]
      sim.results[i,] = c(as.numeric(beta.next), as.numeric(diag(Var.est.formula)), status.convergence)
    } else {
      sim.results[i, ] = rep(NA, 27)
    }
    #probability.cases[i,] = c(mean(cohort.study$Y), mean(CC.study$Y))
    no_of_iter_outer = res.sim$no.of.iter.outer
    
  }
  print(paste0("sim.number.",i))
}

saveRDS(sim.results, file="/Users/prosenjitkundu/Desktop/Code_biometrics_github/Simulation/Table2_results/genmeta_realdata_simulation_balanced_w_additional_interaction.rds")
CI.L = sim.results[,1:10] - 1.96*sqrt(sim.results[,11:20])
CI.U = sim.results[,1:10] + 1.96*sqrt(sim.results[,11:20])
length(which(CI.L[,2] <= beta_real_data[2] & CI.U[,2] >= beta_real_data[2]))
length(which(CI.L[,3] <= beta_real_data[3] & CI.U[,3] >= beta_real_data[3]))
length(which(CI.L[,4] <= beta_real_data[4] & CI.U[,4] >= beta_real_data[4]))
length(which(CI.L[,5] <= beta_real_data[5] & CI.U[,5] >= beta_real_data[5]))
length(which(CI.L[,6] <= beta_real_data[6] & CI.U[,6] >= beta_real_data[6]))
length(which(CI.L[,7] <= beta_real_data[7] & CI.U[,7] >= beta_real_data[7]))
length(which(CI.L[,8] <= beta_real_data[8] & CI.U[,8] >= beta_real_data[8]))
length(which(CI.L[,9] <= beta_real_data[9] & CI.U[,9] >= beta_real_data[9]))
length(which(CI.L[,10] <= beta_real_data[10] & CI.U[,10] >= beta_real_data[10]))
apply(CI.U-CI.L, 2, mean)



