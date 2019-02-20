#----- Where there is measurement error and the design is balanced---#
library(mvtnorm)
library(magic)
library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)
library(expm)
sourceCpp("matrix_multiply.cpp")
source("myoptim.R")

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



# ---- Calculating qudratic form ------ #
Q <- function(beta, study.list, C)
{
  cohort.study = study.list[[1]]
  #cohort.study = study.list[[1]]
  nested.case.control.study = study.list[[2]]
  #nested.case.control.study = study.list[[2]]
  
  
  #fit.R = glm(cohort.study$Y ~ cohort.study$V3 + cohort.study$stage2 + cohort.study$stage3 + cohort.study$stage4 + cohort.study$age, family = binomial())
  fit.R = glm(Y ~ V3, data = cohort.study, family = binomial())
  theta.R = as.numeric(fit.R$coefficients)
  #status.converged = fit.R$converged
  
  no.of.cases.CC.inst1 = length(which(nested.case.control.study$Y == 1 & nested.case.control.study$V3 == 0))
  no.of.cases.CC.inst2 = length(which(nested.case.control.study$Y == 1 & nested.case.control.study$V3 == 1))
  no.of.controls.CC.inst1 = length(which(nested.case.control.study$Y == 0 & nested.case.control.study$V3 == 0))
  no.of.controls.CC.inst2 = length(which(nested.case.control.study$Y == 0 & nested.case.control.study$V3 == 1))
  
  no.of.cases.R.inst1 = length(which(cohort.study$Y == 1 & cohort.study$V3 == 0))
  no.of.cases.R.inst2 = length(which(cohort.study$Y == 1 & cohort.study$V3 == 1))
  no.of.controls.R.inst1 = length(which(cohort.study$Y == 0 & cohort.study$V3 == 0))
  no.of.controls.R.inst2 = length(which(cohort.study$Y == 0 & cohort.study$V3 == 1))
  
  cases.CC.inst1 = nested.case.control.study[which(nested.case.control.study$Y == 1 & nested.case.control.study$V3 == 0), ]
  cases.CC.inst2 = nested.case.control.study[which(nested.case.control.study$Y == 1 & nested.case.control.study$V3 == 1), ]
  controls.CC.inst1 = nested.case.control.study[which(nested.case.control.study$Y == 0 & nested.case.control.study$V3 == 0), ]
  controls.CC.inst2 = nested.case.control.study[which(nested.case.control.study$Y == 0 & nested.case.control.study$V3 == 1),]
  
  Y.nested.case.control.study = nested.case.control.study$Y
  Y.cases.CC.inst1 = cases.CC.inst1$Y
  Y.cases.CC.inst2 = cases.CC.inst2$Y
  Y.controls.CC.inst1 = controls.CC.inst1$Y
  Y.controls.CC.inst2 = controls.CC.inst2$Y
  
  cases.CC.inst1 = as.matrix(cases.CC.inst1)
  cases.CC.inst2 = as.matrix(cases.CC.inst2)
  controls.CC.inst1 = as.matrix(controls.CC.inst1)
  controls.CC.inst2 = as.matrix(controls.CC.inst2)
  #nested.case.control.study = as.matrix(nested.case.control.study)
  
  U1.cases.inst1 = t(cases.CC.inst1[, c(2,3)]) %*% ((expit(beta, cases.CC.inst1[,c(2,6,4,5)]) - expit(theta.R, cases.CC.inst1[, c(2,3)]))*(no.of.cases.R.inst1/no.of.cases.CC.inst1))
  U1.cases.inst2 = t(cases.CC.inst2[, c(2,3)]) %*% ((expit(beta, cases.CC.inst2[,c(2,6,4,5)]) - expit(theta.R, cases.CC.inst2[, c(2,3)]))*(no.of.cases.R.inst2/no.of.cases.CC.inst2))
  U1.controls.inst1 = t(controls.CC.inst1[, c(2,3)]) %*% ((expit(beta, controls.CC.inst1[,c(2,6,4,5)]) - expit(theta.R, controls.CC.inst1[, c(2,3)]))*(no.of.controls.R.inst1/no.of.controls.CC.inst1))
  U1.controls.inst2 = t(controls.CC.inst2[, c(2,3)]) %*% ((expit(beta, controls.CC.inst2[,c(2,6,4,5)]) - expit(theta.R, controls.CC.inst2[, c(2,3)]))*(no.of.controls.R.inst2/no.of.controls.CC.inst2))
  U1 = (U1.cases.inst1 + U1.cases.inst2 + U1.controls.inst1 + U1.controls.inst2)/study.list[[3]]
  
  p1 = study.list[[2]][study.list[[2]]$Y == 1, ]$sampling.fraction[1]
  p2 = study.list[[2]][study.list[[2]]$Y == 0, ]$sampling.fraction[1]
  correction.factor = c(log(p1/p2), rep(0, (length(beta) - 1)))
  theta.CC = theta.cases = theta.controls = beta + correction.factor
  
  U2 = (t(nested.case.control.study[, c(2,6,4,5)]) %*% (nested.case.control.study$Y - expit(theta.CC, nested.case.control.study[, c(2,6,4,5)])))/study.list[[4]]
  
  U = c(U1,U2)
  #C = diag(length(U))
  
  as.numeric(t(U) %*% C %*% U)
  
}


#----- Calculating SE ------#

#---- Calculating Omega -----#

Omega = function(beta, study.list, phase_I_var, phase_II_var)
{
  cohort.study = study.list[[1]]
  nested.case.control.study = study.list[[2]]
  
  phase_I_index = match(phase_I_var, colnames(study.list[[2]]))
  phase_II_index = match(phase_II_var, colnames(study.list[[2]]))
  
  fit.R = glm(Y ~ X1.cont, data = cohort.study, family = binomial())
  theta.R = as.numeric(fit.R$coefficients)
  
  correction_factor_var_index = match(correction_factor_var, colnames(study.list[[2]]))
  
  beta_old_matrix = replicate(dim(study.list[[2]])[1], beta)
  beta_old_matrix[1, ] = beta_old_matrix[1,] + log(study.list[[2]][, correction_factor_var_index])
  theta_CC = t(beta_old_matrix)
  
  no.of.cases.CC.inst1 = length(which(nested.case.control.study$Y == 1 & nested.case.control.study$V3 == 0))
  no.of.cases.CC.inst2 = length(which(nested.case.control.study$Y == 1 & nested.case.control.study$V3 == 1))
  no.of.controls.CC.inst1 = length(which(nested.case.control.study$Y == 0 & nested.case.control.study$V3 == 0))
  no.of.controls.CC.inst2 = length(which(nested.case.control.study$Y == 0 & nested.case.control.study$V3 == 1))
  
  no.of.cases.R.inst1 = length(which(cohort.study$Y == 1 & cohort.study$V3 == 0))
  no.of.cases.R.inst2 = length(which(cohort.study$Y == 1 & cohort.study$V3 == 1))
  no.of.controls.R.inst1 = length(which(cohort.study$Y == 0 & cohort.study$V3 == 0))
  no.of.controls.R.inst2 = length(which(cohort.study$Y == 0 & cohort.study$V3 == 1))
  
  cases.CC.inst1 = nested.case.control.study[which(nested.case.control.study$Y == 1 & nested.case.control.study$V3 == 0), ]
  cases.CC.inst2 = nested.case.control.study[which(nested.case.control.study$Y == 1 & nested.case.control.study$V3 == 1), ]
  controls.CC.inst1 = nested.case.control.study[which(nested.case.control.study$Y == 0 & nested.case.control.study$V3 == 0), ]
  controls.CC.inst2 = nested.case.control.study[which(nested.case.control.study$Y == 0 & nested.case.control.study$V3 == 1),]
  
  theta.cases.CC.inst1 = theta_CC[which(nested.case.control.study$Y == 1 & nested.case.control.study$V3 == 0), ]
  theta.cases.CC.inst2 = theta_CC[which(nested.case.control.study$Y == 1 & nested.case.control.study$V3 == 1), ]
  theta.controls.CC.inst1 = theta_CC[which(nested.case.control.study$Y == 0 & nested.case.control.study$V3 == 0), ]
  theta.controls.CC.inst2 = theta_CC[which(nested.case.control.study$Y == 0 & nested.case.control.study$V3 == 1), ]
  
  Y.nested.case.control.study = nested.case.control.study$Y
  Y.cases.CC.inst1 = cases.CC.inst1$Y
  Y.cases.CC.inst2 = cases.CC.inst2$Y
  Y.controls.CC.inst1 = controls.CC.inst1$Y
  Y.controls.CC.inst2 = controls.CC.inst2$Y
  
  cases.CC.inst1 = as.matrix(cases.CC.inst1)
  cases.CC.inst2 = as.matrix(cases.CC.inst2)
  controls.CC.inst1 = as.matrix(controls.CC.inst1)
  controls.CC.inst2 = as.matrix(controls.CC.inst2)
  nested.case.control.study = as.matrix(nested.case.control.study)
  
  #correction.factor = c(log((no.of.cases.CC * no.of.controls.R)/(no.of.cases.R * no.of.controls.CC)), rep(0, (length(beta) - 1)))
  #p1 = study.list[[2]][study.list[[2]]$Y == 1, ]$sampling.fraction[1]
  #p2 = study.list[[2]][study.list[[2]]$Y == 0, ]$sampling.fraction[1]
  #correction.factor = c(log(p1/p2), rep(0, (length(beta) - 1)))
  #theta.CC = theta.cases = theta.controls = beta + correction.factor
  
  
  
  Omega_11.cases.inst1 = (t(cases.CC.inst1[,  phase_I_index]) %*% diag((expit(beta, cases.CC.inst1[,phase_II_index]) - expit(theta.R, cases.CC.inst1[,  phase_I_index]))^2*((no.of.cases.R.inst1/no.of.cases.CC.inst1)/cases.CC.inst1[,7])) %*% cases.CC.inst1[,  phase_I_index])/study.list[[3]]
  Omega_11.cases.inst2 = (t(cases.CC.inst2[,  phase_I_index]) %*% diag((expit(beta, cases.CC.inst2[,phase_II_index]) - expit(theta.R, cases.CC.inst2[,  phase_I_index]))^2*((no.of.cases.R.inst2/no.of.cases.CC.inst2)/cases.CC.inst2[,7])) %*% cases.CC.inst2[,  phase_I_index])/study.list[[3]]
  Omega_11.controls.inst1 = (t(controls.CC.inst1[, phase_I_index]) %*% diag((expit(beta, controls.CC.inst1[,phase_II_index]) - expit(theta.R, controls.CC.inst1[,  phase_I_index]))^2*((no.of.controls.R.inst1/no.of.controls.CC.inst1)/controls.CC.inst1[,7])) %*% controls.CC.inst1[,  phase_I_index])/study.list[[3]]
  Omega_11.controls.inst2 = (t(controls.CC.inst2[, phase_I_index]) %*% diag((expit(beta, controls.CC.inst2[,phase_II_index]) - expit(theta.R, controls.CC.inst2[,  phase_I_index]))^2*((no.of.controls.R.inst2/no.of.controls.CC.inst2)/controls.CC.inst2[,7])) %*% controls.CC.inst2[,  phase_I_index])/study.list[[3]]
  Omega_11 = Omega_11.cases.inst1 + Omega_11.cases.inst2 + Omega_11.controls.inst1 + Omega_11.controls.inst2
  
  
  #Omega_22 = (t(nested.case.control.study[, 3:6]) %*% diag((Y.nested.case.control.study - expit(theta.CC, nested.case.control.study[, 3:6]))^2) %*% nested.case.control.study[, 3:6])/study.list[[3]]
  Omega_22.cases.inst1 = (t(cases.CC.inst1[, phase_II_index]) %*% diag((Y.cases.CC.inst1 - myexpit(rowSums(theta.cases.CC.inst1 * cases.CC.inst1[,phase_II_index])))^2 * cases.CC.inst1[,7]) %*% cases.CC.inst1[, phase_II_index])  * ((no.of.cases.R.inst1/study.list[[3]])/no.of.cases.CC.inst1)
  Omega_22.cases.inst2 = (t(cases.CC.inst2[, phase_II_index]) %*% diag((Y.cases.CC.inst2 - myexpit(rowSums(theta.cases.CC.inst2 * cases.CC.inst2[,phase_II_index])))^2  * cases.CC.inst2[,7]) %*% cases.CC.inst2[, phase_II_index]) * ((no.of.cases.R.inst2/study.list[[3]])/no.of.cases.CC.inst2)
  Omega_22.controls.inst1 = (t(controls.CC.inst1[, phase_II_index]) %*% diag((Y.controls.CC.inst1 - myexpit(rowSums(theta.controls.CC.inst1 * controls.CC.inst1[,phase_II_index])))^2 * controls.CC.inst1[,7]) %*% controls.CC.inst1[,phase_II_index])  * ((no.of.controls.R.inst1/study.list[[3]])/no.of.controls.CC.inst1)
  Omega_22.controls.inst2 = (t(controls.CC.inst2[, phase_II_index]) %*% diag((Y.controls.CC.inst2 - myexpit(rowSums(theta.controls.CC.inst2 * controls.CC.inst2[,phase_II_index])))^2 * controls.CC.inst2[,7]) %*% controls.CC.inst2[, phase_II_index])  * ((no.of.controls.R.inst2/study.list[[3]])/no.of.controls.CC.inst2)
  Omega_22 = Omega_22.cases.inst1 + Omega_22.cases.inst2 + Omega_22.controls.inst1 + Omega_22.controls.inst2
  
  Omega_33.cases.inst1 = t(cases.CC.inst1[,  phase_I_index]) %*% diag((Y.cases.CC.inst1 - expit(theta.R, cases.CC.inst1[, phase_I_index]))^2) %*% cases.CC.inst1[,  phase_I_index] * ((no.of.cases.R.inst1/study.list[[3]])/no.of.cases.CC.inst1)
  Omega_33.cases.inst2 = t(cases.CC.inst2[, phase_I_index]) %*% diag((Y.cases.CC.inst2 - expit(theta.R, cases.CC.inst2[,  phase_I_index]))^2) %*% cases.CC.inst2[,  phase_I_index] * ((no.of.cases.R.inst2/study.list[[3]])/no.of.cases.CC.inst2)
  Omega_33.controls.inst1 = t(controls.CC.inst1[,  phase_I_index]) %*% diag((Y.controls.CC.inst1 - expit(theta.R, controls.CC.inst1[,  phase_I_index]))^2) %*% controls.CC.inst1[,  phase_I_index] * ((no.of.controls.R.inst1/study.list[[3]])/no.of.controls.CC.inst1)
  Omega_33.controls.inst2 = t(controls.CC.inst2[,  phase_I_index]) %*% diag((Y.controls.CC.inst2 - expit(theta.R, controls.CC.inst2[,  phase_I_index]))^2) %*% controls.CC.inst2[,  phase_I_index] * ((no.of.controls.R.inst2/study.list[[3]])/no.of.controls.CC.inst2)
  Omega_33 = Omega_33.cases.inst1 + Omega_33.cases.inst2 + Omega_33.controls.inst1 + Omega_33.controls.inst2
  
  Omega_12.cases.inst1 = t(cases.CC.inst1[,  phase_I_index]) %*% diag((expit(beta, cases.CC.inst1[, phase_II_index]) - expit(theta.R, cases.CC.inst1[,  phase_I_index]))*(Y.cases.CC.inst1 - myexpit(rowSums(theta.cases.CC.inst1 * cases.CC.inst1[,phase_II_index])))) %*% cases.CC.inst1[, phase_II_index] * ((no.of.cases.R.inst1/study.list[[3]])/no.of.cases.CC.inst1)
  Omega_12.cases.inst2 = t(cases.CC.inst2[,  phase_I_index]) %*% diag((expit(beta, cases.CC.inst2[, phase_II_index]) - expit(theta.R, cases.CC.inst2[,  phase_I_index]))*(Y.cases.CC.inst2 - myexpit(rowSums(theta.cases.CC.inst2 * cases.CC.inst2[,phase_II_index])))) %*% cases.CC.inst2[, phase_II_index] * ((no.of.cases.R.inst2/study.list[[3]])/no.of.cases.CC.inst2)
  Omega_12.controls.inst1 = t(controls.CC.inst1[,  phase_I_index]) %*% diag((expit(beta, controls.CC.inst1[, phase_II_index]) - expit(theta.R, controls.CC.inst1[, phase_I_index]))*(Y.controls.CC.inst1 - myexpit(rowSums(theta.controls.CC.inst1 * controls.CC.inst1[,phase_II_index])))) %*% controls.CC.inst1[, phase_II_index] * ((no.of.controls.R.inst1/study.list[[3]])/no.of.controls.CC.inst1)
  Omega_12.controls.inst2= t(controls.CC.inst2[,  phase_I_index]) %*% diag((expit(beta, controls.CC.inst2[, phase_II_index]) - expit(theta.R, controls.CC.inst2[, phase_I_index]))*(Y.controls.CC.inst2 - myexpit(rowSums(theta.controls.CC.inst2 * controls.CC.inst2[,phase_II_index])))) %*% controls.CC.inst2[, phase_II_index] * ((no.of.controls.R.inst2/study.list[[3]])/no.of.controls.CC.inst2)
  Omega_12 = Omega_12.cases.inst1 + Omega_12.cases.inst2 + Omega_12.controls.inst1 + Omega_12.controls.inst2
  
  #Omega_23 = (t(nested.case.control.study[, 3:6]) %*% diag((Y.nested.case.control.study - expit(theta.CC, nested.case.control.study[, 3:6]))*(Y.nested.case.control.study - expit(theta.R, nested.case.control.study[, 3:4]))) %*% nested.case.control.study[, 3:4])/study.list[[3]]
  Omega_23.cases.inst1 = (t(cases.CC.inst1[, phase_II_index]) %*% diag((Y.cases.CC.inst1 - myexpit(rowSums(theta.cases.CC.inst1 * cases.CC.inst1[,phase_II_index])))*(Y.cases.CC.inst1 - expit(theta.R, cases.CC.inst1[,  phase_I_index])) * cases.CC.inst1[,7]) %*% cases.CC.inst1[,  phase_I_index]) * ((no.of.cases.R.inst1/study.list[[3]])/no.of.cases.CC.inst1)
  Omega_23.cases.inst2 = (t(cases.CC.inst2[, phase_II_index]) %*% diag((Y.cases.CC.inst2 - myexpit(rowSums(theta.cases.CC.inst2 * cases.CC.inst2[,phase_II_index])))*(Y.cases.CC.inst2 - expit(theta.R, cases.CC.inst2[,  phase_I_index])) * cases.CC.inst2[,7]) %*% cases.CC.inst2[,  phase_I_index]) * ((no.of.cases.R.inst2/study.list[[3]])/no.of.cases.CC.inst2)
  Omega_23.controls.inst1 = (t(controls.CC.inst1[, phase_II_index]) %*% diag((Y.controls.CC.inst1 - myexpit(rowSums(theta.controls.CC.inst1 * controls.CC.inst1[,phase_II_index])))*(Y.controls.CC.inst1 - expit(theta.R, controls.CC.inst1[,  phase_I_index])) * controls.CC.inst1[,7]) %*% controls.CC.inst1[,  phase_I_index]) * ((no.of.controls.R.inst1/study.list[[3]])/no.of.controls.CC.inst1)
  Omega_23.controls.inst2 = (t(controls.CC.inst2[, phase_II_index]) %*% diag((Y.controls.CC.inst2 - myexpit(rowSums(theta.controls.CC.inst2 * controls.CC.inst2[,phase_II_index])))*(Y.controls.CC.inst2 - expit(theta.R, controls.CC.inst2[,  phase_I_index]))* controls.CC.inst2[,7]) %*% controls.CC.inst2[,  phase_I_index]) * ((no.of.controls.R.inst2/study.list[[3]])/no.of.controls.CC.inst2)
  Omega_23 = Omega_23.cases.inst1 + Omega_23.cases.inst2 + Omega_23.controls.inst1 + Omega_23.controls.inst2
  
  Omega_13.cases.inst1 = t(cases.CC.inst1[,  phase_I_index]) %*% diag((expit(beta, cases.CC.inst1[, phase_II_index]) - expit(theta.R, cases.CC.inst1[,  phase_I_index]))*(Y.cases.CC.inst1 - expit(theta.R, cases.CC.inst1[,  phase_I_index]))) %*% cases.CC.inst1[,  phase_I_index] * ((no.of.cases.R.inst1/study.list[[3]])/no.of.cases.CC.inst1)
  Omega_13.cases.inst2 = t(cases.CC.inst2[,  phase_I_index]) %*% diag((expit(beta, cases.CC.inst2[, phase_II_index]) - expit(theta.R, cases.CC.inst2[,  phase_I_index]))*(Y.cases.CC.inst2 - expit(theta.R, cases.CC.inst2[,  phase_I_index]))) %*% cases.CC.inst2[,  phase_I_index] * ((no.of.cases.R.inst2/study.list[[3]])/no.of.cases.CC.inst2)
  Omega_13.controls.inst1 = t(controls.CC.inst1[,  phase_I_index]) %*% diag((expit(beta, controls.CC.inst1[, phase_II_index]) - expit(theta.R, controls.CC.inst1[,  phase_I_index]))*(Y.controls.CC.inst1 - expit(theta.R, controls.CC.inst1[,  phase_I_index]))) %*% controls.CC.inst1[,  phase_I_index] * ((no.of.controls.R.inst1/study.list[[3]])/no.of.controls.CC.inst1)
  Omega_13.controls.inst2 = t(controls.CC.inst2[,  phase_I_index]) %*% diag((expit(beta, controls.CC.inst2[, phase_II_index]) - expit(theta.R, controls.CC.inst2[,  phase_I_index]))*(Y.controls.CC.inst2 - expit(theta.R, controls.CC.inst2[,  phase_I_index]))) %*% controls.CC.inst2[,  phase_I_index] * ((no.of.controls.R.inst2/study.list[[3]])/no.of.controls.CC.inst2)
  Omega_13 = Omega_13.cases.inst1 + Omega_13.cases.inst2 + Omega_13.controls.inst1 + Omega_13.controls.inst2
  
  Omega.est = rbind(cbind(Omega_11, Omega_12, Omega_13), cbind(t(Omega_12), Omega_22, Omega_23), cbind(t(Omega_13), t(Omega_23), Omega_33))
  
  Omega.est
  
  
}


#---- Calculating Delta -----#
#---- Calculating Delta -----#
Delta = function(beta, study.list, lambda.in.Delta)
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

Gamma = function(beta, study.list, phase_I_var, phase_II_var, lambda.in.Delta)
{
  cohort.study = study.list[[1]]
  nested.case.control.study = study.list[[2]]
  
  phase_I_index = match(phase_I_var, colnames(study.list[[2]]))
  phase_II_index = match(phase_II_var, colnames(study.list[[2]]))
  
  #fit.R = glm(cohort.study$Y ~ cohort.study$V3 + cohort.study$stage2 + cohort.study$stage3 + cohort.study$stage4 + cohort.study$age, family = binomial())
  fit.R = glm(Y~X1.cont, data = cohort.study, family = binomial())
  theta.R = as.numeric(fit.R$coefficients)
  #status.converged = fit.R$converged
  
  correction_factor_var_index = match(correction_factor_var, colnames(study.list[[2]]))
  
  beta_old_matrix = replicate(dim(study.list[[2]])[1], beta)
  beta_old_matrix[1, ] = beta_old_matrix[1,] + log(study.list[[2]][, correction_factor_var_index])
  theta_CC = t(beta_old_matrix)
  
  no.of.cases.CC.inst1 = length(which(nested.case.control.study$Y == 1 & nested.case.control.study$V3 == 0))
  no.of.cases.CC.inst2 = length(which(nested.case.control.study$Y == 1 & nested.case.control.study$V3 == 1))
  no.of.controls.CC.inst1 = length(which(nested.case.control.study$Y == 0 & nested.case.control.study$V3 == 0))
  no.of.controls.CC.inst2 = length(which(nested.case.control.study$Y == 0 & nested.case.control.study$V3 == 1))
  
  no.of.cases.R.inst1 = length(which(cohort.study$Y == 1 & cohort.study$V3 == 0))
  no.of.cases.R.inst2 = length(which(cohort.study$Y == 1 & cohort.study$V3 == 1))
  no.of.controls.R.inst1 = length(which(cohort.study$Y == 0 & cohort.study$V3 == 0))
  no.of.controls.R.inst2 = length(which(cohort.study$Y == 0 & cohort.study$V3 == 1))
  
  cases.CC.inst1 = nested.case.control.study[which(nested.case.control.study$Y == 1 & nested.case.control.study$V3 == 0), ]
  cases.CC.inst2 = nested.case.control.study[which(nested.case.control.study$Y == 1 & nested.case.control.study$V3 == 1), ]
  controls.CC.inst1 = nested.case.control.study[which(nested.case.control.study$Y == 0 & nested.case.control.study$V3 == 0), ]
  controls.CC.inst2 = nested.case.control.study[which(nested.case.control.study$Y == 0 & nested.case.control.study$V3 == 1),]
  
  theta.cases.CC.inst1 = theta_CC[which(nested.case.control.study$Y == 1 & nested.case.control.study$V3 == 0), ]
  theta.cases.CC.inst2 = theta_CC[which(nested.case.control.study$Y == 1 & nested.case.control.study$V3 == 1), ]
  theta.controls.CC.inst1 = theta_CC[which(nested.case.control.study$Y == 0 & nested.case.control.study$V3 == 0), ]
  theta.controls.CC.inst2 = theta_CC[which(nested.case.control.study$Y == 0 & nested.case.control.study$V3 == 1), ]
  #Y.nested.case.control.study = nested.case.control.study$Y
  #Y.cases.CC.inst1 = cases.CC.inst1$Y
  #Y.cases.CC.inst2 = cases.CC.inst2$Y
  #Y.controls.CC.inst1 = controls.CC.inst1$Y
  #Y.controls.CC.inst2 = controls.CC.inst2$Y
  
  cases.CC.inst1 = as.matrix(cases.CC.inst1)
  cases.CC.inst2 = as.matrix(cases.CC.inst2)
  controls.CC.inst1 = as.matrix(controls.CC.inst1)
  controls.CC.inst2 = as.matrix(controls.CC.inst2)
  nested.case.control.study = as.matrix(nested.case.control.study)
  
  #correction.factor = c(log((no.of.cases.CC * no.of.controls.R)/(no.of.cases.R * no.of.controls.CC)), rep(0, (length(beta) - 1)))
  #p1 = study.list[[2]][study.list[[2]]$Y == 1, ]$sampling.fraction[1]
  #p2 = study.list[[2]][study.list[[2]]$Y == 0, ]$sampling.fraction[1]
  #correction.factor = c(log(p1/p2), rep(0, (length(beta) - 1)))
  #theta.CC = theta.cases = theta.controls = beta + correction.factor
  
  #print(head(cases.CC[, 3:6]))
  #print(class(cases.CC[, 3:6]))
  #print(expit_d1(beta, cases.CC[, 3:6]))
  #print(((no.of.cases.R/study.list[[3]])/no.of.cases.CC))
  
  
  Gamma.1.cases.inst1 = (t(cases.CC.inst1[, phase_I_index]) %*% diag(expit_d1(beta, cases.CC.inst1[, phase_II_index])) %*% cases.CC.inst1[, phase_II_index]) * ((no.of.cases.R.inst1/study.list[[3]])/no.of.cases.CC.inst1)
  Gamma.1.cases.inst2 = (t(cases.CC.inst2[, phase_I_index]) %*% diag(expit_d1(beta, cases.CC.inst2[, phase_II_index])) %*% cases.CC.inst2[, phase_II_index]) * ((no.of.cases.R.inst2/study.list[[3]])/no.of.cases.CC.inst2)
  Gamma.1.controls.inst1 = (t(controls.CC.inst1[, phase_I_index]) %*% diag(expit_d1(beta, controls.CC.inst1[, phase_II_index])) %*% controls.CC.inst1[, phase_II_index]) * ((no.of.controls.R.inst1/study.list[[3]])/no.of.controls.CC.inst1)
  Gamma.1.controls.inst2 = (t(controls.CC.inst2[, phase_I_index]) %*% diag(expit_d1(beta, controls.CC.inst2[, phase_II_index])) %*% controls.CC.inst2[, phase_II_index]) * ((no.of.controls.R.inst2/study.list[[3]])/no.of.controls.CC.inst2)
  Gamma.1 = Gamma.1.cases.inst1 + Gamma.1.cases.inst2 + Gamma.1.controls.inst1 + Gamma.1.controls.inst2
  
  #Gamma.2 = -(t(nested.case.control.study[, 3:6]) %*% diag(expit_d1(theta.CC, nested.case.control.study[, 3:6])) %*% nested.case.control.study[, 3:6])/study.list[[4]]
  
  
  
  Gamma.2.cases.inst1 = (t(cases.CC.inst1[, phase_II_index ]) %*% diag(myexpit_d1(rowSums(theta.cases.CC.inst1 * cases.CC.inst1[,phase_II_index])) * cases.CC.inst1[,7]) %*% cases.CC.inst1[, phase_II_index]) * ((no.of.cases.R.inst1/study.list[[3]])/no.of.cases.CC.inst1)
  Gamma.2.cases.inst2 = (t(cases.CC.inst2[, phase_II_index ]) %*% diag(myexpit_d1(rowSums(theta.cases.CC.inst2 * cases.CC.inst2[,phase_II_index])) * cases.CC.inst2[,7]) %*% cases.CC.inst2[, phase_II_index]) * ((no.of.cases.R.inst2/study.list[[3]])/no.of.cases.CC.inst2)
  Gamma.2.controls.inst1 = (t(controls.CC.inst1[, phase_II_index ]) %*% diag(myexpit_d1(rowSums(theta.controls.CC.inst1 * controls.CC.inst1[,phase_II_index])) * controls.CC.inst1[,7]) %*% controls.CC.inst1[, phase_II_index]) * ((no.of.controls.R.inst1/study.list[[3]])/no.of.controls.CC.inst1)
  Gamma.2.controls.inst2 = (t(controls.CC.inst2[, phase_II_index]) %*% diag(myexpit_d1(rowSums(theta.controls.CC.inst2 * controls.CC.inst2[,phase_II_index])) * controls.CC.inst2[,7]) %*% controls.CC.inst2[, phase_II_index]) * ((no.of.controls.R.inst2/study.list[[3]])/no.of.controls.CC.inst2)
  Gamma.2 = -(1/lambda.in.Delta)*(Gamma.2.cases.inst1 + Gamma.2.cases.inst2 + Gamma.2.controls.inst1 + Gamma.2.controls.inst2)
  
  
  Gamma.est = rbind(Gamma.1, Gamma.2)
  
  Gamma.est
  
}


#----------Calculating C_opt-----#
C_opt <- function(beta, study.list, lambda.in.Delta)
{
  # cohort.study = study.list[[5]]
  # cohort.study = study.list[[1]]
  # nested.case.control.study = study.list[[6]]
  # nested.case.control.study = study.list[[2]]
  # 
  # 
  # #fit.R = glm(cohort.study$Y ~ cohort.study$V3 + cohort.study$stage2 + cohort.study$stage3 + cohort.study$stage4 + cohort.study$age, family = binomial())
  # fit.R = glm(Y~V3, data = cohort.study, family = binomial())
  # theta.R = as.numeric(fit.R$coefficients)
  # #status.converged = fit.R$converged
  # 
  # no.of.cases.CC = sum(nested.case.control.study$Y)
  # no.of.controls.CC = study.list[[4]] - no.of.cases.CC
  # 
  # no.of.cases.R = sum(cohort.study$Y)
  # no.of.controls.R = study.list[[3]] - no.of.cases.R
  # 
  # cases.CC = nested.case.control.study[which(nested.case.control.study$Y == 1),]
  # controls.CC = nested.case.control.study[which(nested.case.control.study$Y == 0),]
  # 
  # Y.nested.case.control.study = nested.case.control.study$Y
  # Y.cases.CC = cases.CC$Y
  # Y.controls.CC = controls.CC$Y
  # 
  # cases.CC = as.matrix(cases.CC)
  # controls.CC = as.matrix(controls.CC)
  # nested.case.control.study = as.matrix(nested.case.control.study)
  # 
  # C_11.cases = t(cases.CC[, c(2,3)]) %*% diag((expit(beta, cases.CC[, c(2,6,4,5)]) - expit(theta.R, cases.CC[, c(2,3)]))^2) %*% cases.CC[, c(2,3)] * ((no.of.cases.R/study.list[[3]])/no.of.cases.CC)
  # C_11.controls = t(controls.CC[, c(2,3)]) %*% diag((expit(beta, controls.CC[, c(2,6,4,5)]) - expit(theta.R, controls.CC[, c(2,3)]))^2) %*% controls.CC[, c(2,3)] * ((no.of.controls.R/study.list[[3]])/no.of.controls.CC)
  # C_11 = C_11.cases + C_11.controls
  # 
  # p1 = study.list[[2]][study.list[[2]]$Y == 1, ]$sampling.fraction[1]
  # p2 = study.list[[2]][study.list[[2]]$Y == 0, ]$sampling.fraction[1]
  # correction.factor = c(log(p1/p2), rep(0, (length(beta) - 1)))
  # theta.CC = theta.cases = theta.controls = beta + correction.factor
  # 
  # C_21.cases = (t(cases.CC[, c(2,6,4,5)]) %*% diag((Y.cases.CC - expit(theta.CC, cases.CC[, c(2,6,4,5)]))*(expit(beta, cases.CC[, c(2,6,4,5)]) - expit(theta.R, cases.CC[, c(2,3)]))) %*% cases.CC[, c(2,3)]) * ((no.of.cases.R/study.list[[3]])/no.of.cases.CC)
  # C_21.controls = (t(controls.CC[, c(2,6,4,5)]) %*% diag((Y.controls.CC - expit(theta.CC, controls.CC[, c(2,6,4,5)]))*(expit(beta, controls.CC[, c(2,6,4,5)]) - expit(theta.R, controls.CC[, c(2,3)]))) %*% controls.CC[, c(2,3)]) * ((no.of.controls.R/study.list[[3]])/no.of.controls.CC)
  # C_21 = C_21.cases + C_21.controls
  # 
  # C_12 = t(C_21)
  # 
  # 
  # C_22.cases = (t(cases.CC[, c(2,6,4,5)]) %*% diag((Y.cases.CC - expit(theta.CC, cases.CC[, c(2,6,4,5)]))^2) %*% cases.CC[, c(2,6,4,5)]) * ((no.of.cases.R/study.list[[3]])/no.of.cases.CC)
  # C_22.controls = (t(controls.CC[, c(2,6,4,5)]) %*% diag((Y.controls.CC - expit(theta.CC, controls.CC[, c(2,6,4,5)]))^2) %*% controls.CC[, c(2,6,4,5)]) * ((no.of.controls.R/study.list[[3]])/no.of.controls.CC) 
  # C_22 = C_22.cases + C_22.controls
  # 
  # C_optim = rbind(cbind(C_11, C_12), cbind(C_21, C_22))
  # 
  # return(solve(C_optim, tol=1e-60))
  # 
  
  return(solve(Delta(beta, study.list, lambda.in.Delta) %*% Omega(beta, study.list, phase_I_var, phase_II_var) %*% t(Delta(beta, study.list, lambda.in.Delta)), tol=1e-60))
}



divison = function(a,b){b[a==1]/b[a==0]}


beta.true = c(-3, 0, log(1.3), log(1.3),  log(1.3))

probability.cases = matrix(NA,100,2)
sim.results = matrix(NA, 1000, 9)
for(i in 1:1000)
{
  set.seed(i)
  
  X = rmvnorm(10000, c(0,0,0), matrix(c(1,0.3,0.6,0.3,1,0.1,0.6,0.1,1),3,3))
  X1.cont = X[,1]
  X[which(X[,1] <= 0), 1] = 0
  X[which(X[,1] > 0), 1] = 1
  X.design = cbind(1, X, X1.cont)
  p = expit(beta.true, X.design)
  set.seed(i)
  Y = rbinom(10000, 1, p)
  
  cohort.study = data.frame(cbind(Y, X.design))
  
  fit.full.cohort = glm(Y ~ X1.cont + V4 + V5, family = binomial(), data = cohort.study)
  
  
  sampling.fraction.cases.strata1 = 1
  sampling.fraction.cases.strata2 = 1
  sampling.fraction.controls.strata1 = sum(cohort.study$Y == 1 & cohort.study$V3 == 1)/sum(cohort.study$Y == 0 & cohort.study$V3 == 0)
  sampling.fraction.controls.strata2 = sum(cohort.study$Y == 1 & cohort.study$V3 == 0)/sum(cohort.study$Y == 0 & cohort.study$V3 == 1)
  
  cohort.study$sampling.fraction = NA
  cohort.study$sampling.fraction[which(cohort.study$Y == 1 & cohort.study$V3 == 0)] = sampling.fraction.cases.strata1
  cohort.study$sampling.fraction[which(cohort.study$Y == 1 & cohort.study$V3 == 1)] = sampling.fraction.cases.strata2
  cohort.study$sampling.fraction[which(cohort.study$Y == 0 & cohort.study$V3 == 0)] = sampling.fraction.controls.strata1
  cohort.study$sampling.fraction[which(cohort.study$Y == 0 & cohort.study$V3 == 1)] = sampling.fraction.controls.strata2
  
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
  
  CC.size <- dim(CC.study)[1]
  cohort.size <- dim(cohort.study)[1]
  
  disease_var = "Y"
  samp_fraction_var = "sampling.fraction"
  strata_var = "V3"
  
  table.disease.strata.samp.frac = aggregate(as.formula(paste(samp_fraction_var, "~", paste(c(disease_var, strata_var), collapse = "+"))), data=CC.study, FUN="unique")
  strata_by_samp_frac = table.disease.strata.samp.frac %>% group_by(V3) %>% summarise(correction.factor = divison(Y,sampling.fraction))
  
  CC.study = full_join(CC.study, strata_by_samp_frac, by = strata_var)
  
  testthat::expect_equal(CC.size,dim(CC.study)[1])
  
  CC.study$NDS = NA
  CC.study$nDS = NA
  CC.study$NDS[which(CC.study$Y == 1 & CC.study$V3 == 0)] = length(R.cases.strata1)
  CC.study$NDS[which(CC.study$Y == 1 & CC.study$V3 == 1)] = length(R.cases.strata2)
  CC.study$NDS[which(CC.study$Y == 0 & CC.study$V3 == 0)] = length(R.controls.strata1)
  CC.study$NDS[which(CC.study$Y == 0 & CC.study$V3 == 1)] = length(R.controls.strata2)
  CC.study$nDS[which(CC.study$Y == 1 & CC.study$V3 == 0)] = sum(R.cases.strata1)
  CC.study$nDS[which(CC.study$Y == 1 & CC.study$V3 == 1)] = sum(R.cases.strata2)
  CC.study$nDS[which(CC.study$Y == 0 & CC.study$V3 == 0)] = sum(R.controls.strata1)
  CC.study$nDS[which(CC.study$Y == 0 & CC.study$V3 == 1)] = sum(R.controls.strata2)
  CC.study$sampling.fraction.est = CC.study$nDS/CC.study$NDS
  
  
  Q.study.list <- list(cohort.study, CC.study, cohort.size, CC.size)
  
  lambda.in.Delta = ((sampling.fraction.cases.strata1 * length(R.cases.strata1)) + (sampling.fraction.cases.strata2 * length(R.cases.strata2)) + (sampling.fraction.controls.strata1 * length(R.controls.strata1))  + (sampling.fraction.controls.strata2 * length(R.controls.strata2)))/Q.study.list[[3]] 
  
  
  correction_factor_var = "correction.factor"
  
  Q.beta.initial.cohort = as.numeric(fit.full.cohort$coefficients)
  #Q.beta.initial.cohort = append(Q.beta.initial.cohort, 0, after = 1)
  #Q.beta.initial.CC = as.numeric(fit.CC$coefficients)
  #Q.beta.initial.CC = append(Q.beta.initial.CC, 0, after = 1)
  fit.R = glm(Y~ X1.cont, data = cohort.study, family = binomial())
  theta.R = fit.R$coefficients
  
  phase_I_var = c("V2","X1.cont")
  phase_II_var = c("V2","X1.cont", "V4", "V5")
  sampling_fraction_est_var = c("sampling.fraction.est")
  no_of_iter_outer = 0
  eps = 1e-06
  C = diag(length(phase_I_var) + length(phase_II_var))
  res.sim = myoptim(Q.study.list, C, Q.beta.initial.cohort, eps, no_of_iter_outer, phase_I_var, phase_II_var, theta.R, correction_factor_var, sampling_fraction_est_var)
  if(res.sim$Status == 0)
  {
    sim.results[i, ] = rep(NA, 9)
  }else{
    beta.initial = as.numeric(res.sim$beta_optim)
    C = C_opt(as.numeric(res.sim$beta_optim), Q.study.list, lambda.in.Delta)
    repeat{
      res.sim = myoptim(Q.study.list, C, beta.initial, eps, no_of_iter_outer, phase_I_var, phase_II_var, theta.R, correction_factor_var, sampling_fraction_est_var)
      #res.sim = optimx(Q.beta.initial.cohort, Q, study.list = Q.study.list, C = C, theta_R = theta.R, phase_I_var=phase_I_var, phase_II_var=phase_II_var, correction_factor_var = correction_factor_var, sampling_fraction_est_var = sampling_fraction_est_var, method = "BFGS")
      #beta.next <- as.numeric(res.sim)[1:4]
      beta.next = as.numeric(res.sim$beta_optim)
      C = C_opt(beta.next, Q.study.list, lambda.in.Delta)
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
      #C = C_opt(beta.next, Q.study.list, lambda.in.Delta)
      #Var.est.formula.1 = solve((t(Gamma(beta.next, Q.study.list)) %*% C %*% Gamma(beta.next, Q.study.list)), tol=1e-30) %*% t(Gamma(beta.next, Q.study.list)) %*% C %*% Delta(Q.study.list)
      #Var.est.formula = (Var.est.formula.1 %*% Omega(beta.next, Q.study.list) %*% t(Var.est.formula.1))/Q.study.list[[4]]
      #sim.results[i,] = c(as.numeric(beta.next), as.numeric(diag(Var.est.formula)), as.numeric(Q.beta.initial.CC), as.numeric(diag(vcov(fit.CC))), status.convergence)
      #Var.est.formula = solve(t(Gamma(beta.next, Q.study.list)) %*% C %*% Gamma(beta.next, Q.study.list), tol = 1e-60 )/Q.study.list[[3]]
      Var.est.formula = solve(t(Gamma(beta.next, Q.study.list, phase_I_var, phase_II_var, lambda.in.Delta)) %*% C %*% Gamma(beta.next, Q.study.list, phase_I_var, phase_II_var, lambda.in.Delta), tol = 1e-60 )/Q.study.list[[3]]
      sim.results[i,] = c(as.numeric(beta.next), as.numeric(diag(Var.est.formula)), status.convergence)
    } else {
      sim.results[i, ] = rep(NA, 9)
    }
    #probability.cases[i,] = c(mean(cohort.study$Y), mean(CC.study$Y))
    no_of_iter_outer = res.sim$no.of.iter.outer
    
  }
  print(paste0("sim.number.",i))
}
saveRDS(sim.results, file="GENMETA_balanced_m_II.rds")

CI.L = sim.results[,1:4] - 1.96*sqrt(sim.results[,5:8])
CI.U = sim.results[,1:4] + 1.96*sqrt(sim.results[,5:8])
length(which(CI.L[,2] <= beta.true[3] & CI.U[,2] >= beta.true[3]))
length(which(CI.L[,3] <= beta.true[3] & CI.U[,3] >= beta.true[3]))
length(which(CI.L[,4] <= beta.true[3] & CI.U[,4] >= beta.true[3]))
apply(CI.U-CI.L, 2, mean)


