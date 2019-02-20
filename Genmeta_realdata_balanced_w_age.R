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










cohort.study = nwtco
cohort.study$histol <- factor(cohort.study$histol)
cohort.study$instit <- factor(cohort.study$instit)
cohort.study$stage <- factor(cohort.study$stage)
cohort.study$age  = as.numeric(scale(cohort.study$age/12))
cohort.strata.table <- table(cohort.study$rel, cohort.study$instit)

sampling.fraction.cases.strata1 = 1
sampling.fraction.cases.strata2 = 1
sampling.fraction.controls.strata1 = (sum(which(cohort.study$rel == 1)) - sum(which(cohort.study$rel == 0 & cohort.study$instit == 2)))/sum(which(cohort.study$rel == 0)) 
sampling.fraction.controls.strata2 = 1

cohort.study$sampling.fraction = NA
cohort.study$sampling.fraction[which(cohort.study$rel == 1 & cohort.study$instit == 1)] = sampling.fraction.cases.strata1
cohort.study$sampling.fraction[which(cohort.study$rel == 1 & cohort.study$instit == 2)] = sampling.fraction.cases.strata2
cohort.study$sampling.fraction[which(cohort.study$rel == 0 & cohort.study$instit == 1)] = sampling.fraction.controls.strata1
cohort.study$sampling.fraction[which(cohort.study$rel == 0 & cohort.study$instit == 2)] = sampling.fraction.controls.strata2

divison = function(a,b){b[a==1]/b[a==0]}
start.time = Sys.time()
sim.results = matrix(NA, 1000, 19)
cohort.study.design.matrix <- data.frame(model.matrix(~ rel + instit + histol + stage + age + instit*stage*age + histol*stage*age, data = cohort.study))
phase_I_var = c("X.Intercept.","instit2", "stage2", "stage3", "stage4", "age", "instit2.stage2","instit2.stage3", "instit2.stage4")
phase_II_var = c("X.Intercept.","histol2","stage2","stage3","stage4","age","histol2.stage2","histol2.stage3", "histol2.stage4")
fit.full.cohort = glm(as.formula(paste("rel", "~", paste(phase_II_var[-1], collapse = "+"))), family = binomial(), data = cohort.study.design.matrix)
Q.beta.initial.cohort = as.numeric(fit.full.cohort$coefficients)

for(i in 1:1000)
{
  
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
  
  CC.study <- CC.study[order(CC.study$seqno), ]
  CC.size <- dim(CC.study)[1]
  cohort.size <- dim(cohort.study)[1]
  
  disease_var = "rel"
  #samp_fraction_var = "sampling.fraction"
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
    sim.results[i, ] = rep(NA, 19)
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
      sim.results[i, ] = rep(NA, 19)
    }
    #probability.cases[i,] = c(mean(cohort.study$Y), mean(CC.study$Y))
    no_of_iter_outer = res.sim$no.of.iter.outer
    
  }
  print(paste0("sim.number.",i))
}
end.time = Sys.time()
saveRDS(sim.results, file="real_data_GENMETA_w_age_balanced.rds")

