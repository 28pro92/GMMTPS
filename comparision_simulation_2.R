#----- Where there is measurement error---#
library(mvtnorm)
library(magic)
library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)
library(expm)
sourceCpp("matrix_multiply.cpp")


expit <- function(beta, X)
{
  return(as.numeric(1/(1 + exp(-as.matrix(X) %*% beta))))
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
  U2 =  t(as.matrix(phase_II)) %*% as.vector(study.list[[2]]$Y - as.vector(1/(1 + exp(-rowSums(ref_dat * theta_CC)))))
  
  U = c(U1,U2)
  #C = diag(length(U))
  
  as.numeric(t(U) %*% C %*% U)
  
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
  b2 =  study.list[[2]]$Y - as.vector(1/(1 + exp(-rowSums(ref_dat * theta_CC))))
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
  
  as.vector(2*Dn)
  
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
  b2 =  study.list[[2]]$Y - as.vector(1/(1 + exp(-rowSums(ref_dat * theta_CC))))
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
  
  2*J_n
  
  
  
}






#----- Calculating SE ------#

#---- Calculating Omega -----#

Omega = function(beta, study.list, phase_I_var, phase_II_var)
{
  cohort.study = study.list[[1]]
  nested.case.control.study = study.list[[2]]
  
  phase_I_index = match(phase_I_var, colnames(study.list[[2]]))
  phase_II_index = match(phase_II_var, colnames(study.list[[2]]))
  
  #fit.R = glm(cohort.study$Y ~ cohort.study$V3 + cohort.study$stage2 + cohort.study$stage3 + cohort.study$stage4 + cohort.study$age, family = binomial())
  fit.R = glm(Y~X1.cont, data = cohort.study, family = binomial())
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
  nested.case.control.study = as.matrix(nested.case.control.study)
  
  #correction.factor = c(log((no.of.cases.CC * no.of.controls.R)/(no.of.cases.R * no.of.controls.CC)), rep(0, (length(beta) - 1)))
  p1 = study.list[[2]][study.list[[2]]$Y == 1, ]$sampling.fraction[1]
  p2 = study.list[[2]][study.list[[2]]$Y == 0, ]$sampling.fraction[1]
  correction.factor = c(log(p1/p2), rep(0, (length(beta) - 1)))
  theta.CC = theta.cases = theta.controls = beta + correction.factor
  
  Omega_11.cases.inst1 = crossprod((cases.CC.inst1[, phase_I_index]*(expit(beta, cases.CC.inst1[,phase_II_index]) - expit(theta.R, cases.CC.inst1[, phase_I_index]))^2*((no.of.cases.R.inst1/no.of.cases.CC.inst1)/p1)), cases.CC.inst1[,phase_I_index])/study.list[[3]]
  Omega_11.cases.inst2 = crossprod((cases.CC.inst2[, phase_I_index]*(expit(beta, cases.CC.inst2[,phase_II_index]) - expit(theta.R, cases.CC.inst2[, phase_I_index]))^2*((no.of.cases.R.inst2/no.of.cases.CC.inst2)/p1)), cases.CC.inst2[,phase_I_index])/study.list[[3]]
  #Omega_11.cases.inst1 = (t(cases.CC.inst1[, phase_I_index]) %*% diag((expit(beta, cases.CC.inst1[,phase_II_index]) - expit(theta.R, cases.CC.inst1[, phase_I_index]))^2*((no.of.cases.R.inst1/no.of.cases.CC.inst1)/p1)) %*% cases.CC.inst1[,phase_I_index])/study.list[[3]]
  #Omega_11.cases.inst2 = (t(cases.CC.inst2[, phase_I_index]) %*% diag((expit(beta, cases.CC.inst2[,phase_II_index]) - expit(theta.R, cases.CC.inst2[, phase_I_index]))^2*((no.of.cases.R.inst2/no.of.cases.CC.inst2)/p1)) %*% cases.CC.inst2[, phase_I_index])/study.list[[3]]
  Omega_11.controls.inst1  = crossprod((controls.CC.inst1[, phase_I_index]*(expit(beta, controls.CC.inst1[,phase_II_index]) - expit(theta.R, controls.CC.inst1[, phase_I_index]))^2*((no.of.controls.R.inst1/no.of.controls.CC.inst1)/p2)), controls.CC.inst1[,phase_I_index])/study.list[[3]]
  Omega_11.controls.inst2  = crossprod((controls.CC.inst2[, phase_I_index]*(expit(beta, controls.CC.inst2[,phase_II_index]) - expit(theta.R, controls.CC.inst2[, phase_I_index]))^2*((no.of.controls.R.inst2/no.of.controls.CC.inst2)/p2)), controls.CC.inst2[,phase_I_index])/study.list[[3]]
  #Omega_11.controls.inst1 = (t(controls.CC.inst1[,phase_I_index]) %*% diag((expit(beta, controls.CC.inst1[,phase_II_index]) - expit(theta.R, controls.CC.inst1[, phase_I_index]))^2*((no.of.controls.R.inst1/no.of.controls.CC.inst1)/p2)) %*% controls.CC.inst1[, phase_I_index])/study.list[[3]]
  #Omega_11.controls.inst2 = (t(controls.CC.inst2[,phase_I_index]) %*% diag((expit(beta, controls.CC.inst2[,phase_II_index]) - expit(theta.R, controls.CC.inst2[, phase_I_index]))^2*((no.of.controls.R.inst2/no.of.controls.CC.inst2)/p2)) %*% controls.CC.inst2[, phase_I_index])/study.list[[3]]
  Omega_11 = Omega_11.cases.inst1 + Omega_11.cases.inst2 + Omega_11.controls.inst1 + Omega_11.controls.inst2
  
  
  #Omega_22 = (t(nested.case.control.study[, 3:6]) %*% diag((Y.nested.case.control.study - expit(theta.CC, nested.case.control.study[, 3:6]))^2) %*% nested.case.control.study[, 3:6])/study.list[[3]]
  Omega_22.cases.inst1 = crossprod(cases.CC.inst1[, phase_II_index] * (Y.cases.CC.inst1 - expit(theta.CC, cases.CC.inst1[, phase_II_index]))^2,cases.CC.inst1[, phase_II_index]) * (((no.of.cases.R.inst1/study.list[[3]])/no.of.cases.CC.inst1)*p1)
  Omega_22.cases.inst2 = crossprod(cases.CC.inst2[, phase_II_index] * (Y.cases.CC.inst2 - expit(theta.CC, cases.CC.inst2[, phase_II_index]))^2,cases.CC.inst2[, phase_II_index]) * (((no.of.cases.R.inst2/study.list[[3]])/no.of.cases.CC.inst2)*p1)
  #Omega_22.cases.inst1 = (t(cases.CC.inst1[, phase_II_index]) %*% diag((Y.cases.CC.inst1 - expit(theta.CC, cases.CC.inst1[, phase_II_index]))^2) %*% cases.CC.inst1[, phase_II_index]) * (((no.of.cases.R.inst1/study.list[[3]])/no.of.cases.CC.inst1)*p1)
  #Omega_22.cases.inst2 = (t(cases.CC.inst2[,phase_II_index]) %*% diag((Y.cases.CC.inst2 - expit(theta.CC, cases.CC.inst2[, phase_II_index]))^2) %*% cases.CC.inst2[, phase_II_index]) * (((no.of.cases.R.inst2/study.list[[3]])/no.of.cases.CC.inst2)*p1)
  Omega_22.controls.inst1 = crossprod(controls.CC.inst1[, phase_II_index] * (Y.controls.CC.inst1 - expit(theta.CC, controls.CC.inst1[, phase_II_index]))^2,controls.CC.inst1[, phase_II_index]) * (((no.of.controls.R.inst1/study.list[[3]])/no.of.controls.CC.inst1)*p2)
  Omega_22.controls.inst2 = crossprod(controls.CC.inst2[, phase_II_index] * (Y.controls.CC.inst2 - expit(theta.CC, controls.CC.inst2[, phase_II_index]))^2,controls.CC.inst2[, phase_II_index]) * (((no.of.controls.R.inst2/study.list[[3]])/no.of.controls.CC.inst2)*p2)
  #Omega_22.controls.inst1 = (t(controls.CC.inst1[, phase_II_index]) %*% diag((Y.controls.CC.inst1 - expit(theta.CC, controls.CC.inst1[, phase_II_index]))^2) %*% controls.CC.inst1[, phase_II_index]) * (((no.of.controls.R.inst1/study.list[[3]])/no.of.controls.CC.inst1)*p2)
  #Omega_22.controls.inst2 = (t(controls.CC.inst2[,phase_II_index]) %*% diag((Y.controls.CC.inst2 - expit(theta.CC, controls.CC.inst2[, phase_II_index]))^2) %*% controls.CC.inst2[, phase_II_index]) * (((no.of.controls.R.inst2/study.list[[3]])/no.of.controls.CC.inst2)*p2)
  Omega_22 = Omega_22.cases.inst1 + Omega_22.cases.inst2 + Omega_22.controls.inst1 + Omega_22.controls.inst2
  
  
  Omega_33.cases.inst1 = crossprod(cases.CC.inst1[, phase_I_index] * (Y.cases.CC.inst1 - expit(theta.R, cases.CC.inst1[, phase_I_index]))^2, cases.CC.inst1[, phase_I_index]) * ((no.of.cases.R.inst1/study.list[[3]])/no.of.cases.CC.inst1)
  Omega_33.cases.inst2 = crossprod(cases.CC.inst2[, phase_I_index] * (Y.cases.CC.inst2 - expit(theta.R, cases.CC.inst2[, phase_I_index]))^2, cases.CC.inst2[, phase_I_index]) * ((no.of.cases.R.inst2/study.list[[3]])/no.of.cases.CC.inst2)
  #Omega_33.cases.inst1 = t(cases.CC.inst1[, phase_I_index]) %*% diag((Y.cases.CC.inst1 - expit(theta.R, cases.CC.inst1[, phase_I_index]))^2) %*% cases.CC.inst1[, phase_I_index] * ((no.of.cases.R.inst1/study.list[[3]])/no.of.cases.CC.inst1)
  #Omega_33.cases.inst2 = t(cases.CC.inst2[, phase_I_index]) %*% diag((Y.cases.CC.inst2 - expit(theta.R, cases.CC.inst2[,phase_I_index]))^2) %*% cases.CC.inst2[, phase_I_index] * ((no.of.cases.R.inst2/study.list[[3]])/no.of.cases.CC.inst2)
  Omega_33.controls.inst1 = crossprod(controls.CC.inst1[, phase_I_index] * (Y.controls.CC.inst1- expit(theta.R, controls.CC.inst1[, phase_I_index]))^2, controls.CC.inst1[, phase_I_index]) * ((no.of.controls.R.inst1/study.list[[3]])/no.of.controls.CC.inst1)
  Omega_33.controls.inst2 = crossprod(controls.CC.inst2[, phase_I_index] * (Y.controls.CC.inst2- expit(theta.R, controls.CC.inst2[, phase_I_index]))^2, controls.CC.inst2[, phase_I_index]) * ((no.of.controls.R.inst2/study.list[[3]])/no.of.controls.CC.inst2)
  #Omega_33.controls.inst1 = t(controls.CC.inst1[, phase_I_index]) %*% diag((Y.controls.CC.inst1 - expit(theta.R, controls.CC.inst1[, phase_I_index]))^2) %*% controls.CC.inst1[, phase_I_index] * ((no.of.controls.R.inst1/study.list[[3]])/no.of.controls.CC.inst1)
  #Omega_33.controls.inst2 = t(controls.CC.inst2[, phase_I_index]) %*% diag((Y.controls.CC.inst2 - expit(theta.R, controls.CC.inst2[, phase_I_index]))^2) %*% controls.CC.inst2[, phase_I_index] * ((no.of.controls.R.inst2/study.list[[3]])/no.of.controls.CC.inst2)
  Omega_33 = Omega_33.cases.inst1 + Omega_33.cases.inst2 + Omega_33.controls.inst1 + Omega_33.controls.inst2
  #Omega_33 = as.matrix(solve(vcov(fit.R)* study.list[[3]]))
  
  Omega_12.cases.inst1 = t(cases.CC.inst1[, phase_I_index]) %*% diag((expit(beta, cases.CC.inst1[, phase_II_index]) - expit(theta.R, cases.CC.inst1[, phase_I_index]))*(Y.cases.CC.inst1 - expit(theta.CC, cases.CC.inst1[, phase_II_index]))) %*% cases.CC.inst1[, phase_II_index] * ((no.of.cases.R.inst1/study.list[[3]])/no.of.cases.CC.inst1)
  Omega_12.cases.inst2 = t(cases.CC.inst2[, phase_I_index]) %*% diag((expit(beta, cases.CC.inst2[, phase_II_index]) - expit(theta.R, cases.CC.inst2[, phase_I_index]))*(Y.cases.CC.inst2 - expit(theta.CC, cases.CC.inst2[, phase_II_index]))) %*% cases.CC.inst2[, phase_II_index] * ((no.of.cases.R.inst2/study.list[[3]])/no.of.cases.CC.inst2)
  Omega_12.controls.inst1 = t(controls.CC.inst1[, phase_I_index]) %*% diag((expit(beta, controls.CC.inst1[, phase_II_index]) - expit(theta.R, controls.CC.inst1[,phase_I_index]))*(Y.controls.CC.inst1 - expit(theta.CC, controls.CC.inst1[, phase_II_index]))) %*% controls.CC.inst1[, phase_II_index] * ((no.of.controls.R.inst1/study.list[[3]])/no.of.controls.CC.inst1)
  Omega_12.controls.inst2= t(controls.CC.inst2[,phase_I_index]) %*% diag((expit(beta, controls.CC.inst2[, phase_II_index]) - expit(theta.R, controls.CC.inst2[,phase_I_index]))*(Y.controls.CC.inst2 - expit(theta.CC, controls.CC.inst2[, phase_II_index]))) %*% controls.CC.inst2[,phase_II_index] * ((no.of.controls.R.inst2/study.list[[3]])/no.of.controls.CC.inst2)
  Omega_12 = Omega_12.cases.inst1 + Omega_12.cases.inst2 + Omega_12.controls.inst1 + Omega_12.controls.inst2
  
  #Omega_23 = (t(nested.case.control.study[, 3:6]) %*% diag((Y.nested.case.control.study - expit(theta.CC, nested.case.control.study[, 3:6]))*(Y.nested.case.control.study - expit(theta.R, nested.case.control.study[, 3:4]))) %*% nested.case.control.study[, 3:4])/study.list[[3]]
  Omega_23.cases.inst1 = (t(cases.CC.inst1[, phase_II_index]) %*% diag((Y.cases.CC.inst1 - expit(theta.CC, cases.CC.inst1[, phase_II_index]))*(Y.cases.CC.inst1 - expit(theta.R, cases.CC.inst1[, phase_I_index]))) %*% cases.CC.inst1[, phase_I_index]) * (((no.of.cases.R.inst1/study.list[[3]])/no.of.cases.CC.inst1)*p1)
  Omega_23.cases.inst2 = (t(cases.CC.inst2[, phase_II_index]) %*% diag((Y.cases.CC.inst2 - expit(theta.CC, cases.CC.inst2[, phase_II_index]))*(Y.cases.CC.inst2 - expit(theta.R, cases.CC.inst2[, phase_I_index]))) %*% cases.CC.inst2[,phase_I_index]) * (((no.of.cases.R.inst2/study.list[[3]])/no.of.cases.CC.inst2)*p1)
  Omega_23.controls.inst1 = (t(controls.CC.inst1[, phase_II_index]) %*% diag((Y.controls.CC.inst1 - expit(theta.CC, controls.CC.inst1[, phase_II_index]))*(Y.controls.CC.inst1 - expit(theta.R, controls.CC.inst1[, phase_I_index]))) %*% controls.CC.inst1[, phase_I_index]) * (((no.of.controls.R.inst1/study.list[[3]])/no.of.controls.CC.inst1)*p2)
  Omega_23.controls.inst2 = (t(controls.CC.inst2[, phase_II_index]) %*% diag((Y.controls.CC.inst2 - expit(theta.CC, controls.CC.inst2[, phase_II_index]))*(Y.controls.CC.inst2 - expit(theta.R, controls.CC.inst2[, phase_I_index]))) %*% controls.CC.inst2[, phase_I_index]) * (((no.of.controls.R.inst2/study.list[[3]])/no.of.controls.CC.inst2)*p2)
  Omega_23 = Omega_23.cases.inst1 + Omega_23.cases.inst2 + Omega_23.controls.inst1 + Omega_23.controls.inst2
  
  Omega_13.cases.inst1 = t(cases.CC.inst1[, phase_I_index]) %*% diag((expit(beta, cases.CC.inst1[, phase_II_index]) - expit(theta.R, cases.CC.inst1[, phase_I_index]))*(Y.cases.CC.inst1 - expit(theta.R, cases.CC.inst1[, phase_I_index]))) %*% cases.CC.inst1[, phase_I_index] * ((no.of.cases.R.inst1/study.list[[3]])/no.of.cases.CC.inst1)
  Omega_13.cases.inst2 = t(cases.CC.inst2[, phase_I_index]) %*% diag((expit(beta, cases.CC.inst2[, phase_II_index]) - expit(theta.R, cases.CC.inst2[, phase_I_index]))*(Y.cases.CC.inst2 - expit(theta.R, cases.CC.inst2[, phase_I_index]))) %*% cases.CC.inst2[, phase_I_index] * ((no.of.cases.R.inst2/study.list[[3]])/no.of.cases.CC.inst2)
  Omega_13.controls.inst1 = t(controls.CC.inst1[, phase_I_index]) %*% diag((expit(beta, controls.CC.inst1[, phase_II_index]) - expit(theta.R, controls.CC.inst1[, phase_I_index]))*(Y.controls.CC.inst1 - expit(theta.R, controls.CC.inst1[, phase_I_index]))) %*% controls.CC.inst1[, phase_I_index] * ((no.of.controls.R.inst1/study.list[[3]])/no.of.controls.CC.inst1)
  Omega_13.controls.inst2 = t(controls.CC.inst2[, phase_I_index]) %*% diag((expit(beta, controls.CC.inst2[, phase_II_index]) - expit(theta.R, controls.CC.inst2[, phase_I_index]))*(Y.controls.CC.inst2 - expit(theta.R, controls.CC.inst2[, phase_I_index]))) %*% controls.CC.inst2[,phase_I_index] * ((no.of.controls.R.inst2/study.list[[3]])/no.of.controls.CC.inst2)
  Omega_13 = Omega_13.cases.inst1 + Omega_13.cases.inst2 + Omega_13.controls.inst1 + Omega_13.controls.inst2
  
  Omega.est = rbind(cbind(Omega_11, Omega_12, Omega_13), cbind(t(Omega_12), Omega_22, Omega_23), cbind(t(Omega_13), t(Omega_23), Omega_33))
  
  Omega.est
  
  
}




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
  p1 = study.list[[2]][study.list[[2]]$Y == 1, ]$sampling.fraction[1]
  p2 = study.list[[2]][study.list[[2]]$Y == 0, ]$sampling.fraction[1]
  correction.factor = c(log(p1/p2), rep(0, (length(beta) - 1)))
  theta.CC = theta.cases = theta.controls = beta + correction.factor
  
  #print(head(cases.CC[, 3:6]))
  #print(class(cases.CC[, 3:6]))
  #print(expit_d1(beta, cases.CC[, 3:6]))
  #print(((no.of.cases.R/study.list[[3]])/no.of.cases.CC))
  
  
  Gamma.1.cases.inst1 = (t(cases.CC.inst1[, phase_I_index]) %*% diag(expit_d1(beta, cases.CC.inst1[, phase_II_index])) %*% cases.CC.inst1[, phase_II_index]) * ((no.of.cases.R.inst1/study.list[[3]])/no.of.cases.CC.inst1)
  Gamma.1.cases.inst2 = (t(cases.CC.inst2[, phase_I_index]) %*% diag(expit_d1(beta, cases.CC.inst2[, phase_II_index])) %*% cases.CC.inst2[, phase_II_index]) * ((no.of.cases.R.inst2/study.list[[3]])/no.of.cases.CC.inst2)
  Gamma.1.controls.inst1 = (t(controls.CC.inst1[,phase_I_index]) %*% diag(expit_d1(beta, controls.CC.inst1[, phase_II_index])) %*% controls.CC.inst1[, phase_II_index]) * ((no.of.controls.R.inst1/study.list[[3]])/no.of.controls.CC.inst1)
  Gamma.1.controls.inst2 = (t(controls.CC.inst2[,phase_I_index]) %*% diag(expit_d1(beta, controls.CC.inst2[, phase_II_index])) %*% controls.CC.inst2[, phase_II_index]) * ((no.of.controls.R.inst2/study.list[[3]])/no.of.controls.CC.inst2)
  Gamma.1 = Gamma.1.cases.inst1 + Gamma.1.cases.inst2 + Gamma.1.controls.inst1 + Gamma.1.controls.inst2
  
  #Gamma.2 = -(t(nested.case.control.study[, 3:6]) %*% diag(expit_d1(theta.CC, nested.case.control.study[, 3:6])) %*% nested.case.control.study[, 3:6])/study.list[[4]]
  
  Gamma.2.cases.inst1 = (t(cases.CC.inst1[,phase_II_index]) %*% diag(expit_d1(theta.CC, cases.CC.inst1[, phase_II_index])) %*% cases.CC.inst1[, phase_II_index]) * (((no.of.cases.R.inst1/study.list[[3]])/no.of.cases.CC.inst1)*p1)
  Gamma.2.cases.inst2 = (t(cases.CC.inst2[, phase_II_index]) %*% diag(expit_d1(theta.CC, cases.CC.inst2[, phase_II_index])) %*% cases.CC.inst2[, phase_II_index]) * (((no.of.cases.R.inst2/study.list[[3]])/no.of.cases.CC.inst2)*p1)
  Gamma.2.controls.inst1 = (t(controls.CC.inst1[, phase_II_index]) %*% diag(expit_d1(theta.CC, controls.CC.inst1[, phase_II_index])) %*% controls.CC.inst1[, phase_II_index]) * (((no.of.controls.R.inst1/study.list[[3]])/no.of.controls.CC.inst1)*p2)
  Gamma.2.controls.inst2 = (t(controls.CC.inst2[, phase_II_index]) %*% diag(expit_d1(theta.CC, controls.CC.inst2[, phase_II_index])) %*% controls.CC.inst2[, phase_II_index]) * (((no.of.controls.R.inst2/study.list[[3]])/no.of.controls.CC.inst2)*p2)
  Gamma.2 = -(1/lambda.in.Delta)*(Gamma.2.cases.inst1 + Gamma.2.cases.inst2 + Gamma.2.controls.inst1 + Gamma.2.controls.inst2)
  
  
  Gamma.est = rbind(Gamma.1, Gamma.2)
  
  Gamma.est
  
}


#----------Calculating C_opt-----#
C_opt <- function(beta, study.list, lambda.in.Delta)
{
  return(solve(Delta(beta, study.list, lambda.in.Delta) %*% Omega(beta, study.list, phase_I_var, phase_II_var) %*% t(Delta(beta, study.list, lambda.in.Delta)), tol=1e-60))
}









beta.true = c(-3, 0, log(1.3), log(1.3),log(1.3))

probability.cases = rep(NA,100)
sim.results = matrix(NA, 1000, 9)
start.time.loop = Sys.time()
for(i in 1:1000)
{
  set.seed(i)
  
  X = rmvnorm(10000, c(0,0,0), matrix(c(1,0.3,0.6,0.3,1,0.1,0.6,0.1,1),3,3))
  X1.cont = X[,1]
  X[,1] = ifelse(X[,1]>0, 1, 0)
  #X[which(X[,1] <= 0), 1] = 0
  #X[which(X[,1] > 0), 1] = 1
  X.design = cbind(1, X, X1.cont)
  p = expit(beta.true, X.design)
  set.seed(i)
  Y = rbinom(10000, 1, p)
  probability.cases[i] = mean(Y)
  cohort.study = data.frame(cbind(Y, X.design))
  
  fit.full.cohort = glm(Y ~ X1.cont + V4 + V5, family = binomial(), data = cohort.study)
  
  
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
  
  CC.study = CC.study %>% mutate(correction.factor = CC.study[CC.study$Y == 1, ]$sampling.fraction[1]/CC.study[CC.study$Y == 0, ]$sampling.fraction[1])
 

  
  
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
  
  
  
  CC.size <- dim(CC.study)[1]
  cohort.size <- dim(cohort.study)[1]
  
  Q.study.list <- list(cohort.study, CC.study, cohort.size, CC.size)
  lambda.in.Delta = ((sampling.fraction.cases.strata1 * length(R.cases.strata1)) + (sampling.fraction.cases.strata2 * length(R.cases.strata2)) + (sampling.fraction.controls.strata1 * length(R.controls.strata1))  + (sampling.fraction.controls.strata2 * length(R.controls.strata2)))/Q.study.list[[3]] 
  
  Q.beta.initial.cohort = as.numeric(fit.full.cohort$coefficients)
  #Q.beta.initial.cohort = append(Q.beta.initial.cohort, 0, after = 1)
  #Q.beta.initial.CC = as.numeric(fit.CC$coefficients)
  #Q.beta.initial.CC = append(Q.beta.initial.CC, 0, after = 1)
  fit.R = glm(Y~ X1.cont, data = cohort.study, family = binomial())
  theta.R = fit.R$coefficients
  
  phase_I_var = c("V2","X1.cont")
  phase_II_var = c("V2","X1.cont", "V4", "V5")
  sampling_fraction_est_var = c("sampling.fraction.est")
  correction_factor_var = "correction.factor"
  no_of_iter_outer = 0
  eps = 1e-06
  C = diag(length(phase_I_var) + length(phase_II_var))
  start.time = Sys.time()
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
  end.time = Sys.time()
  print(paste0("sim.number.",i))
}
end.time.loop = Sys.time()
saveRDS(sim.results, file="GENMETA_case_control_m_II.rds")

CI.L = sim.results[,1:4] - 1.96*sqrt(sim.results[,5:8])
CI.U = sim.results[,1:4] + 1.96*sqrt(sim.results[,5:8])
length(which(CI.L[,2] <= beta.true[3] & CI.U[,2] >= beta.true[3]))
length(which(CI.L[,3] <= beta.true[3] & CI.U[,3] >= beta.true[3]))
length(which(CI.L[,4] <= beta.true[3] & CI.U[,4] >= beta.true[3]))
apply(CI.U-CI.L, 2, mean)






