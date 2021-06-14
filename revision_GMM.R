library(survival)
library(mvtnorm)
library(magic)
library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)
library(expm)
library(dplyr)
sourceCpp("matrix_multiply.cpp")
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


# ---- Calculating qudratic form ------ #
Q <- function(beta_old, theta_R, study.list, phase_I_var, phase_II_var, correction_factor_var, sampling_fraction_est_var, C)
{
  phase_I = study.list[[2]][,match(phase_I_var, colnames(study.list[[2]]))]
  phase_II = study.list[[2]][,match(phase_II_var, colnames(study.list[[2]]))]
  sampling.fraction.est = study.list[[2]][,match(sampling_fraction_est_var, colnames(study.list[[2]]))]
  correction_factor_var_index = match(correction_factor_var, colnames(study.list[[2]]))
  ref_dat = as.matrix(phase_II)
  lambda = study.list[[4]]/study.list[[3]]
  beta_old_matrix = replicate(dim(ref_dat)[1], beta_old)
  beta_old_matrix[1, ] = beta_old_matrix[1,] + log(study.list[[2]][, correction_factor_var_index])
  #theta_CC = beta_old + correction_factor_var
  theta_CC = t(beta_old_matrix)
  U1 = t(as.matrix(phase_I)) %*% (((as.vector(1/(1 + exp(-ref_dat %*% beta_old))) - as.vector(1/(1 + exp(-as.matrix(phase_I) %*% theta_R))))/sampling.fraction.est) * lambda)
  U2 =  t(as.matrix(phase_II)) %*% as.vector(study.list[[2]]$rel - as.vector(1/(1 + exp(-rowSums(ref_dat * theta_CC)))))
  
  U = c(U1,U2)
  #C = diag(length(U))
  
  as.numeric((t(U) %*% C %*% U)/(study.list[[4]]^2))
  
}

mygradient = function(beta_old, theta_R, study.list, phase_I_var, phase_II_var, correction_factor_var, sampling_fraction_est_var, C)
{
  phase_I = as.matrix(study.list[[2]][,match(phase_I_var, colnames(study.list[[2]]))])
  phase_II = as.matrix(study.list[[2]][,match(phase_II_var, colnames(study.list[[2]]))])
  sampling.fraction.est = study.list[[2]][,match(sampling_fraction_est_var, colnames(study.list[[2]]))]
  ref_dat = phase_II
  
  #X_abdiag = Matrix(magic::adiag(as.matrix(phase_I), as.matrix(phase_II)), sparse = TRUE)
  
  #X_rbind = Matrix(rbind(ref_dat, ref_dat), sparse = TRUE)
  correction_factor_var_index = match(correction_factor_var, colnames(study.list[[2]]))
  
  beta_old_matrix = replicate(dim(ref_dat)[1], beta_old)
  beta_old_matrix[1, ] = beta_old_matrix[1,] + log(study.list[[2]][, correction_factor_var_index])
  #theta_CC = beta_old + correction_factor_var
  theta_CC = t(beta_old_matrix)
  
  lambda = study.list[[4]]/study.list[[3]]
  
  w_1 <- (as.vector((1/(1 + exp(-ref_dat %*% beta_old)))*(1/(1 + exp(ref_dat %*% beta_old))))/sampling.fraction.est) * lambda
  w_2 <- -as.vector((1/(1 + exp(-rowSums(ref_dat * theta_CC))))*(1/(1 + exp(rowSums(ref_dat * theta_CC)))))
  #W = Matrix(magic::adiag(diag(w_1), diag(w_2)), sparse = TRUE)
  
  #W <- diag(rep(1/W_temp_1, no_of_studies))
  #print(is.nan(W))
  p_1 = length(phase_I_var)
  p_2 = length(phase_II_var)
  C_11 = C[1:p_1, 1:p_1]
  C_12 = C[1:p_1, (p_1+1):(p_1+p_2)]
  C_21 = C[(p_1+1):(p_1+p_2), 1:p_1]
  C_22 = C[(p_1+1):(p_1+p_2), (p_1+1):(p_1+p_2)]
  
  b1 = ((as.vector(1/(1 + exp(-ref_dat %*% beta_old))) - as.vector(1/(1 + exp(-as.matrix(phase_I) %*% theta_R))))/sampling.fraction.est) * lambda
  b2 =  study.list[[2]]$rel - as.vector(1/(1 + exp(-rowSums(ref_dat * theta_CC))))
  b = c(b1, b2)
  
  e1 <- as.vector(1/(1 + exp(-ref_dat %*% beta_old)))
  e2 <- as.vector((exp(-ref_dat %*% beta_old) - 1)/(exp(-ref_dat %*% beta_old) + 1))
  e3 <- as.vector(1/(1 + exp(ref_dat %*% beta_old)))
  e4 = 1/sampling.fraction.est
  l <- e1*e2*e3*e4 *lambda
  #print(class(l))
  nan_indices <- which(l %in% NaN == TRUE)
  l[nan_indices] <- 0
  
  v_1 = (b1 %*% as.matrix(phase_I) %*% C_11 %*% t(as.matrix(phase_I)) %*% diag(l)) + (b2 %*% as.matrix(phase_II) %*% C_21 %*% t(as.matrix(phase_I)) %*% diag(l))
  #L1 <- diag(l)/study.list[[3]]
  
  
  e1 <- as.vector(1/(1 + exp(-rowSums(ref_dat * theta_CC))))
  e2 <- as.vector((exp(-rowSums(ref_dat * theta_CC)) - 1)/(exp(-rowSums(ref_dat * theta_CC)) + 1))
  e3 <- as.vector(1/(1 + exp(rowSums(ref_dat * theta_CC))))
  l <- e1*e2*e3
  #print(class(l))
  nan_indices <- which(l %in% NaN == TRUE)
  l[nan_indices] <- 0
  #print(class(b2))
  v_2 = (b1 %*% as.matrix(phase_I) %*% C_12 %*% t(as.matrix(phase_II)) %*% (-diag(l))) + (b2 %*% as.matrix(phase_II) %*% C_22 %*% t(as.matrix(phase_II)) %*% (-diag(l)))
  
  rm(l)
  
  
  Dn = (t(ref_dat) %*% Matrix(diag(w_1), sparse = T) %*% phase_I %*% C_11 %*% t(phase_I) %*% b1) +
    (t(ref_dat) %*% Matrix(diag(w_1), sparse = T) %*% phase_I %*% C_12 %*% t(phase_II) %*% b2) +
    (t(ref_dat) %*% Matrix(diag(w_2), sparse = T) %*% phase_II %*% C_21 %*% t(phase_I) %*% b1) +
    (t(ref_dat) %*% Matrix(diag(w_2), sparse = T) %*% phase_II %*% C_22 %*% t(phase_II) %*% b2)
  
  as.vector(2*Dn/study.list[[4]]^2)
  
}

myhessian = function(beta_old, theta_R, study.list, phase_I_var, phase_II_var, correction_factor_var, sampling_fraction_est_var, C)
{
  phase_I = as.matrix(study.list[[2]][,match(phase_I_var, colnames(study.list[[2]]))])
  phase_II = as.matrix(study.list[[2]][,match(phase_II_var, colnames(study.list[[2]]))])
  sampling.fraction.est = study.list[[2]][,match(sampling_fraction_est_var, colnames(study.list[[2]]))]
  ref_dat = phase_II
  
  correction_factor_var_index = match(correction_factor_var, colnames(study.list[[2]]))
  
  beta_old_matrix = replicate(dim(ref_dat)[1], beta_old)
  beta_old_matrix[1, ] = beta_old_matrix[1,] + log(study.list[[2]][, correction_factor_var_index])
  #theta_CC = beta_old + correction_factor_var
  theta_CC = t(beta_old_matrix)
  
  lambda = study.list[[4]]/study.list[[3]]
  
  w_1 <- (as.vector((1/(1 + exp(-ref_dat %*% beta_old)))*(1/(1 + exp(ref_dat %*% beta_old))))/sampling.fraction.est) * lambda
  w_2 <- -as.vector((1/(1 + exp(-rowSums(ref_dat * theta_CC))))*(1/(1 + exp(rowSums(ref_dat * theta_CC)))))
  
  
  #W <- diag(rep(1/W_temp_1, no_of_studies))
  #print(is.nan(W))
  p_1 = length(phase_I_var)
  p_2 = length(phase_II_var)
  C_11 = C[1:p_1, 1:p_1]
  C_12 = C[1:p_1, (p_1+1):(p_1+p_2)]
  C_21 = C[(p_1+1):(p_1+p_2), 1:p_1]
  C_22 = C[(p_1+1):(p_1+p_2), (p_1+1):(p_1+p_2)]
  
  b1 = ((as.vector(1/(1 + exp(-ref_dat %*% beta_old))) - as.vector(1/(1 + exp(-as.matrix(phase_I) %*% theta_R))))/sampling.fraction.est) * lambda
  b2 =  study.list[[2]]$rel - as.vector(1/(1 + exp(-rowSums(ref_dat * theta_CC))))
  b = c(b1, b2)
  
  e1 <- as.vector(1/(1 + exp(-ref_dat %*% beta_old)))
  e2 <- as.vector((exp(-ref_dat %*% beta_old) - 1)/(exp(-ref_dat %*% beta_old) + 1))
  e3 <- as.vector(1/(1 + exp(ref_dat %*% beta_old)))
  e4 = 1/sampling.fraction.est
  l <- e1*e2*e3*e4 *lambda
  #print(class(l))
  nan_indices <- which(l %in% NaN == TRUE)
  l[nan_indices] <- 0
  
  v_1 = (b1 %*% as.matrix(phase_I) %*% C_11 %*% t(as.matrix(phase_I)) %*% diag(l)) + (b2 %*% as.matrix(phase_II) %*% C_21 %*% t(as.matrix(phase_I)) %*% diag(l))
  #L1 <- diag(l)/study.list[[3]]
  
  
  e1 <- as.vector(1/(1 + exp(-rowSums(ref_dat * theta_CC))))
  e2 <- as.vector((exp(-rowSums(ref_dat * theta_CC)) - 1)/(exp(-rowSums(ref_dat * theta_CC)) + 1))
  e3 <- as.vector(1/(1 + exp(rowSums(ref_dat * theta_CC))))
  l <- e1*e2*e3
  #print(class(l))
  nan_indices <- which(l %in% NaN == TRUE)
  l[nan_indices] <- 0
  #print(class(b2))
  v_2 = (b1 %*% as.matrix(phase_I) %*% C_12 %*% t(as.matrix(phase_II)) %*% (-diag(l))) + (b2 %*% as.matrix(phase_II) %*% C_22 %*% t(as.matrix(phase_II)) %*% (-diag(l)))
  
  rm(l)
  
  J_n = t(ref_dat) %*% diag(w_1) %*% as.matrix(phase_I) %*% C_11 %*% t(t(ref_dat) %*% diag(w_1) %*% as.matrix(phase_I)) + 
    t(ref_dat) %*% diag(w_1) %*% as.matrix(phase_I) %*% C_12 %*% t(as.matrix(phase_II)) %*% diag(w_2) %*% ref_dat +
    t(t(ref_dat) %*% diag(w_1) %*% as.matrix(phase_I) %*% C_12 %*% t(as.matrix(phase_II)) %*% diag(w_2) %*% ref_dat)  +
    t(ref_dat) %*% diag(w_2) %*% as.matrix(phase_II) %*% C_22 %*% t(t(ref_dat) %*% diag(w_2) %*% as.matrix(phase_II)) +
    (t(ref_dat) %*% (diag(as.numeric(v_1)) + diag(as.numeric(v_2))) %*% ref_dat)
  
  (2*J_n)/(study.list[[4]]^2)
  
  
  
}

Omega = function(beta, study.list, phase_I_var, phase_II_var)
{
  cohort.study = study.list[[1]]
  nested.case.control.study = study.list[[2]]
  strata = cohort.study$strata
  phase_I.indices = match(phase_I_var, colnames(study.list[[2]]))
  phase_II.indices = match(phase_II_var, colnames(study.list[[2]]))
  
  sampling_fraction_index = match("sampling.fraction.est", colnames(study.list[[2]]))
  fit.R = glm(as.formula(paste("rel", "~", paste(phase_I_var[-1], collapse = "+"))), data = cohort.study, family = binomial())
  theta.R = fit.R$coefficients
  
  Y_index_phase_I = match("rel", colnames(cohort.study))
  Y_index_phase_II = match("rel", colnames(nested.case.control.study))
  correction_factor_var_index = match(correction_factor_var, colnames(study.list[[2]]))
  
  beta_old_matrix = replicate(dim(study.list[[2]])[1], beta)
  beta_old_matrix[1, ] = beta_old_matrix[1,] + log(study.list[[2]][, correction_factor_var_index])
  theta_CC = t(beta_old_matrix)
  
  Omega_11 = list()
  Omega_22 = list()
  Omega_33 = list()
  Omega_12 = list()
  Omega_13 = list()
  Omega_23 = list()
  for(i in 1:length(unique(strata)))
  {
    data_I = cohort.study[which(cohort.study$strata == i), ]
    data_II = nested.case.control.study[which(nested.case.control.study$strata == i), ]
    theta.cases.CC.inst1 = theta_CC[which(nested.case.control.study$strata == i), ]
    Y.cases.CC.inst1 = as.numeric(data_II[, Y_index_phase_II])
    
    data_II = as.matrix(data_II)
    Omega_11[[i]] = (t(data_II[, phase_I.indices ]) %*% diag((expit(beta, data_II[, phase_II.indices]) - expit(theta.R, data_II[, phase_I.indices ]))^2*((nrow(data_I)/nrow(data_II))/data_II[,sampling_fraction_index])) %*% data_II[, phase_I.indices ])/study.list[[3]]
    Omega_22[[i]] = (t(data_II[,  phase_II.indices]) %*% diag((Y.cases.CC.inst1 - myexpit(rowSums(theta.cases.CC.inst1 * data_II[, phase_II.indices])))^2 * data_II[,sampling_fraction_index]) %*% data_II[,  phase_II.indices])  * ((nrow(data_I)/study.list[[3]])/nrow(data_II))
    Omega_33[[i]] = t(data_II[, phase_I.indices ]) %*% diag((Y.cases.CC.inst1 - expit(theta.R, data_II[, phase_I.indices ]))^2) %*% data_II[, phase_I.indices ] * ((nrow(data_I)/study.list[[3]])/nrow(data_II))
    Omega_12[[i]] = t(data_II[, phase_I.indices]) %*% diag((expit(beta, data_II[,  phase_II.indices]) - expit(theta.R, data_II[, phase_I.indices ]))*(Y.cases.CC.inst1 - myexpit(rowSums(theta.cases.CC.inst1 * data_II[, phase_II.indices])))) %*% data_II[,  phase_II.indices] * ((nrow(data_I)/study.list[[3]])/nrow(data_II))
    Omega_23[[i]] = (t(data_II[,  phase_II.indices]) %*% diag((Y.cases.CC.inst1 - myexpit(rowSums(theta.cases.CC.inst1 * data_II[, phase_II.indices])))*(Y.cases.CC.inst1 - expit(theta.R, data_II[, phase_I.indices ])) * data_II[,sampling_fraction_index]) %*% data_II[, phase_I.indices ]) * ((nrow(data_I)/study.list[[3]])/nrow(data_II))
    Omega_13[[i]] = t(data_II[, phase_I.indices ]) %*% diag((expit(beta, data_II[,  phase_II.indices]) - expit(theta.R, data_II[, phase_I.indices ]))*(Y.cases.CC.inst1 - expit(theta.R, data_II[, phase_I.indices ]))) %*% data_II[, phase_I.indices ] * ((nrow(data_I)/study.list[[3]])/nrow(data_II))
  }
  
  Omega_11 = Reduce('+', Omega_11)
  Omega_22 = Reduce('+', Omega_22)
  Omega_33 = Reduce('+', Omega_33)
  Omega_13 = Reduce('+', Omega_13)
  Omega_12 = Reduce('+', Omega_12)
  Omega_23 = Reduce('+', Omega_23)
  
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


Gamma = function(beta, study.list, phase_I_var, phase_II_var)
{
  cohort.study = study.list[[1]]
  nested.case.control.study = study.list[[2]]
  strata = cohort.study$strata
  phase_I.indices = match(phase_I_var, colnames(study.list[[2]]))
  phase_II.indices = match(phase_II_var, colnames(study.list[[2]]))
  
  sampling_fraction_index = match("sampling.fraction.est", colnames(study.list[[2]]))
  fit.R = glm(as.formula(paste("rel", "~", paste(phase_I_var[-1], collapse = "+"))), data = cohort.study, family = binomial())
  theta.R = fit.R$coefficients
  
  correction_factor_var_index = match(correction_factor_var, colnames(study.list[[2]]))
  
  beta_old_matrix = replicate(dim(study.list[[2]])[1], beta)
  beta_old_matrix[1, ] = beta_old_matrix[1,] + log(study.list[[2]][, correction_factor_var_index])
  theta_CC = t(beta_old_matrix)
  
  Gamma.1 = list()
  Gamma.2 = list()
  for(i in 1:length(unique(strata)))
  {
    data_I = cohort.study[which(cohort.study$strata == i), ]
    data_II = nested.case.control.study[which(nested.case.control.study$strata == i), ]
    theta.cases.CC.inst1 = theta_CC[which(nested.case.control.study$strata == i), ]
    
    data_II = as.matrix(data_II)
    Gamma.1[[i]] = (t(data_II[, phase_I.indices]) %*% diag(expit_d1(beta, data_II[,  phase_II.indices])) %*% data_II[,  phase_II.indices]) * ((nrow(data_I)/study.list[[3]])/nrow(data_II))
    Gamma.2[[i]] = (t(data_II[,  phase_II.indices]) %*% diag(myexpit_d1(rowSums(theta.cases.CC.inst1 * data_II[, phase_II.indices])) * data_II[,sampling_fraction_index]) %*% data_II[,  phase_II.indices]) * ((nrow(data_I)/study.list[[3]])/nrow(data_II))
  }
  Gamma.1 = Reduce("+", Gamma.1)
  Gamma.2 = Reduce("+", Gamma.2)
  Gamma.2 = -(1/lambda.in.Delta)*(Gamma.2)
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






sim.results = matrix(NA, 1000, 15)
for(i in 1:1000)
{
  cohort.study = nwtco
  cohort.study$age = as.numeric(cohort.study$age/12)
  cohort.study$age.cat = if_else(cohort.study$age>=1, 1,0)
  
  cohort.study$stage.collapse = if_else(cohort.study$stage==4 | cohort.study$stage==3, 1, 0)
  cohort.study$stage.collapse = factor(cohort.study$stage.collapse)
  
  #cohort.study$age.cat = cut(cohort.study$age, breaks = c(min(cohort.study$age),median(cohort.study$age), max(cohort.study$age)), include.lowest = T)
  cohort.study$strata = 1
  cohort.study$strata[which(cohort.study$rel==0 & cohort.study$age.cat==0 & cohort.study$instit == 1 & (cohort.study$stage==1|cohort.study$stage==2))] = 2
  cohort.study$strata[which(cohort.study$rel==0 & cohort.study$age.cat==1 & cohort.study$instit == 1 & (cohort.study$stage==1|cohort.study$stage==2))] = 3
  cohort.study$strata[which(cohort.study$rel==0 & cohort.study$age.cat==1 & cohort.study$instit == 1 & (cohort.study$stage==3|cohort.study$stage==4))] = 4
  # cohort.study$strata[which(cohort.study$rel==0 & cohort.study$age.cat==0 & cohort.study$instit == 1 & (cohort.study$stage==1|cohort.study$stage==2))] = 2
  # cohort.study$strata[which(cohort.study$rel==0 & cohort.study$age.cat==1 & cohort.study$instit == 1 & (cohort.study$stage==1|cohort.study$stage==2))] = 3
  # cohort.study$strata[which(cohort.study$rel==0 & cohort.study$age.cat==1 & cohort.study$instit == 1 & (cohort.study$stage==3|cohort.study$stage==4))] = 4
  # cohort.study$strata[which(cohort.study$rel==0 & cohort.study$age.cat==0 & cohort.study$instit == 1 & (cohort.study$stage==1|cohort.study$stage==2))] = 2
  # cohort.study$strata[which(cohort.study$rel==0 & cohort.study$age.cat==1 & cohort.study$instit == 1 & (cohort.study$stage==1|cohort.study$stage==2))] = 3
  # cohort.study$strata[which(cohort.study$rel==0 & cohort.study$age.cat==1 & cohort.study$instit == 1 & (cohort.study$stage==3|cohort.study$stage==4))] = 4
  # cohort.study$strata[which(cohort.study$rel==0 & cohort.study$age.cat==0 & cohort.study$instit == 1 & (cohort.study$stage==1|cohort.study$stage==2))] = 2
  # cohort.study$strata[which(cohort.study$rel==0 & cohort.study$age.cat==1 & cohort.study$instit == 1 & (cohort.study$stage==1|cohort.study$stage==2))] = 3
  # cohort.study$strata[which(cohort.study$rel==0 & cohort.study$age.cat==1 & cohort.study$instit == 1 & (cohort.study$stage==3|cohort.study$stage==4))] = 4
  # 
  cohort.study$R = 0
  cohort.study$R[which(cohort.study$strata == 1)] = 1
  set.seed(i)
  cohort.study$R[which(cohort.study$strata == 2)] = rbinom(length(which(cohort.study$strata == 2)), 1, 90/length(which(cohort.study$strata == 2)))
  set.seed(i)
  cohort.study$R[which(cohort.study$strata == 3)] = rbinom(length(which(cohort.study$strata == 3)), 1, 120/length(which(cohort.study$strata == 3)))
  set.seed(i)
  cohort.study$R[which(cohort.study$strata == 4)] = rbinom(length(which(cohort.study$strata == 4)), 1, 90/length(which(cohort.study$strata == 4)))
  
  
  cohort.study$strata = as.numeric(as.factor(paste0(cohort.study$rel, ":", cohort.study$instit, ":", cohort.study$age.cat, ":", cohort.study$stage.collapse)))
  
  #cohort.study$histol <- factor(cohort.study$histol)
  #cohort.study$instit <- factor(cohort.study$instit)
  cohort.study$stage <- factor(cohort.study$stage)
  cohort.study$age = as.numeric(cohort.study$age)
  cohort.study$study = as.factor(cohort.study$study)
  #est.wts.data = estWeights(cohort.study, formula=~instit*stage.impute+age.impute + study,
  #              subset=~I(!is.na(histol)))
  #est.wts.fit = svyglm(histol~instit*stage.impute + age.impute + study, design = est.wts.data, family=quasibinomial())
  cohort.study$in.subsample = TRUE
  cohort.study$in.subsample[which(cohort.study$R==0)] = FALSE
  
  cohort.study$one = 1
  table_Nds = aggregate(one~rel+stage.collapse+instit+age.cat, FUN=sum, data = cohort.study)
  colnames(table_Nds)[5] = "NDS"
  
  CC.study.new = cohort.study[which(cohort.study$in.subsample == T), ]
  table_nds = aggregate(one~rel+stage.collapse+instit+age.cat, FUN=sum, data = CC.study.new)
  colnames(table_nds)[5] = "nDS"
  
  CC.study.new = inner_join(inner_join(CC.study.new, table_Nds, by = c("rel", "stage.collapse", "instit", "age.cat")), table_nds, by = c("rel", "stage.collapse", "instit", "age.cat"))
  CC.study.new$sampling.fraction.est = CC.study.new$nDS/CC.study.new$NDS
  divison = function(a,b){unique(b[a==1])/unique(b[a==0])}
  strata_by_samp_frac = CC.study.new %>% group_by(instit, age.cat, stage.collapse) %>% summarise(correction.factor = divison(rel,sampling.fraction.est))
  CC.study.new = inner_join(CC.study.new, strata_by_samp_frac, by = c("stage.collapse", "instit", "age.cat"))
  CC.study.new$instit = factor(CC.study.new$instit)
  CC.study.new$histol = factor(CC.study.new$histol)
  CC.study.new$stage.collapse = factor(CC.study.new$stage.collapse)
  
  CC.study.design.matrix <- data.frame(model.matrix(~ rel + instit + histol + stage.collapse + age + instit*stage.collapse*age + histol*stage.collapse*age + correction.factor + sampling.fraction.est + strata, data = CC.study.new))
  
  lambda.in.Delta = nrow(CC.study.new)/nrow(cohort.study)
  
  sampling_fraction_var = c("sampling.fraction.est")
  correction_factor_var = "correction.factor"
  
  cohort.study$instit = factor(cohort.study$instit)
  cohort.study$histol = factor(cohort.study$histol)
  cohort.study$stage.collapse = factor(cohort.study$stage.collapse)
  cohort.study.design.matrix <- data.frame(model.matrix(~ rel + instit + histol + stage.collapse + age + instit*stage.collapse*age + histol*stage.collapse*age + strata, data = cohort.study))
  phase_I_var = c("X.Intercept.","instit2", "stage.collapse1", "age", "instit2.age", "instit2.stage.collapse1", "stage.collapse1.age")
  phase_II_var = c("X.Intercept.","histol2","stage.collapse1","age", "histol2.age","histol2.stage.collapse1", "stage.collapse1.age")
  fit.full.cohort = glm(as.formula(paste("rel", "~", paste(phase_II_var[-1], collapse = "+"))), family = binomial(), data = cohort.study.design.matrix)
  Q.beta.initial.cohort = as.numeric(fit.full.cohort$coefficients)
  
  CC.size <- dim(CC.study.new)[1]
  cohort.size <- dim(cohort.study)[1]
  Q.study.list <- list(cohort.study.design.matrix, CC.study.design.matrix, cohort.size, CC.size,  cohort.study, CC.study.new)
  
  fit.R = glm(as.formula(paste("rel", "~", paste(phase_I_var[-1], collapse = "+"))), data = cohort.study.design.matrix, family = binomial())
  theta.R = fit.R$coefficients
  
  
  no_of_iter_outer = 0
  eps = 1e-06
  outer_iter = 0
  C = C_opt(Q.beta.initial.cohort, Q.study.list, phase_I_var, phase_II_var,lambda.in.Delta)
  #C = diag(length(phase_I_var) + length(phase_II_var))
  res.sim = myoptim(Q.study.list, C, Q.beta.initial.cohort, eps, no_of_iter_outer, phase_I_var, phase_II_var, theta.R, correction_factor_var, sampling_fraction_var)
  #res.sim = optimx(Q.beta.initial.cohort, Q, study.list = Q.study.list, C = C, theta_R = theta.R, phase_I_var=phase_I_var, phase_II_var=phase_II_var, correction_factor_var = correction_factor_var, sampling_fraction_est_var = sampling_fraction_var, method = "BFGS")
  
  if(res.sim$Status == 0)
  {
    sim.results[i, ] = rep(NA, 15)
  }else{
    beta.initial = as.numeric(res.sim$beta_optim)
    #beta.initial = as.numeric(res.sim[1:7])
    C = C_opt(beta.initial, Q.study.list, phase_I_var, phase_II_var,lambda.in.Delta)
    repeat{
      res.sim = myoptim(Q.study.list, C, beta.initial, eps, no_of_iter_outer, phase_I_var, phase_II_var, theta.R, correction_factor_var, sampling_fraction_var)
      #res.sim = optimx(Q.beta.initial.cohort, Q, study.list = Q.study.list, C = C, theta_R = theta.R, phase_I_var=phase_I_var, phase_II_var=phase_II_var, correction_factor_var = correction_factor_var, sampling_fraction_est_var = sampling_fraction_var, method = "BFGS")
      #beta.next <- as.numeric(res.sim)[1:7]
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
      sim.results[i, ] = rep(NA, 15)
    }
    #probability.cases[i,] = c(mean(cohort.study$Y), mean(CC.study$Y))
    no_of_iter_outer = res.sim$no.of.iter.outer
    
  }
  print(paste0("sim.number.",i))
  
}

MSE_GMM = apply((sim.results[,1:7] - t(matrix(rep(Q.beta.initial.cohort, 1000), nrow = 7)))^2, 2, mean)
names(MSE_GMM) = names(fit.full.cohort$coefficients)
