myoptim <- function(study.list, C, initial_val, threshold, no_of_iter_outer, phase_I_var, phase_II_var, theta_R, correction_factor_var, sampling_fraction_est_var)
{
  beta_old <- as.vector(initial_val)
  eps_inner <- 0
  phase_I = as.matrix(study.list[[2]][,match(phase_I_var, colnames(study.list[[2]]))])
  phase_II = as.matrix(study.list[[2]][,match(phase_II_var, colnames(study.list[[2]]))])
  ref_dat = phase_II
  threshold_optim = threshold
  status = 1
  #X_abdiag = Matrix(magic::adiag(as.matrix(phase_I), as.matrix(phase_II)), sparse = TRUE)
  #sampling_fraction = study.list[[2]][,match(sampling_fraction_var, colnames(study.list[[2]]))]
  #X_rbind = Matrix(rbind(ref_dat, ref_dat), sparse = TRUE)
  correction_factor_var_index = match(correction_factor_var, colnames(study.list[[2]]))
  sampling.fraction.est = study.list[[2]][,match(sampling_fraction_est_var, colnames(study.list[[2]]))]
  p_1 = length(phase_I_var)
  p_2 = length(phase_II_var)
  C_11 = C[1:p_1, 1:p_1]
  C_12 = C[1:p_1, (p_1+1):(p_1+p_2)]
  C_21 = t(C_12)
  C_22 = C[(p_1+1):(p_1+p_2), (p_1+1):(p_1+p_2)]
  lambda = study.list[[4]]/study.list[[3]]
  #p1 = study.list[[2]][study.list[[2]]$Y == 1, ]$sampling.fraction[1]
  #p2 = study.list[[2]][study.list[[2]]$Y == 0, ]$sampling.fraction[1]
  #correction.factor = c(log(p1/p2), rep(0, (length(beta_old) - 1)))
    iter = 0
    continue = TRUE
    while(continue)
    {
      beta_old_matrix = replicate(dim(ref_dat)[1], beta_old)
      beta_old_matrix[1, ] = beta_old_matrix[1,] + log(study.list[[2]][, correction_factor_var_index])
      #theta_CC = beta_old + correction_factor_var
      theta_CC = t(beta_old_matrix)

      
      
      w_1 <- (as.vector((1/(1 + exp(-ref_dat %*% beta_old)))*(1/(1 + exp(ref_dat %*% beta_old))))/sampling.fraction.est) * lambda
      w_2 <- -as.vector((1/(1 + exp(-rowSums(ref_dat * theta_CC))))*(1/(1 + exp(rowSums(ref_dat * theta_CC)))))
      #W = Matrix(magic::adiag(diag(w_1), diag(w_2)), sparse = TRUE)
      
      #W <- diag(rep(1/W_temp_1, no_of_studies))
      #print(is.nan(W))
      
      
      b1 = ((as.vector(1/(1 + exp(-ref_dat %*% beta_old))) - as.vector(1/(1 + exp(-as.matrix(phase_I) %*% theta_R))))/sampling.fraction.est) * lambda
      b2 =  study.list[[2]]$Y - as.vector(1/(1 + exp(-rowSums(ref_dat * theta_CC))))
      #b = c(b1, b2)
      
      e1 <- as.vector(1/(1 + exp(-ref_dat %*% beta_old)))
      e2 <- as.vector((exp(-ref_dat %*% beta_old) - 1)/(exp(-ref_dat %*% beta_old) + 1))
      e3 <- as.vector(1/(1 + exp(ref_dat %*% beta_old)))
      #e4 = 1/sampling.fraction.est
      l <- (e1*e2*e3*lambda)/sampling.fraction.est
      #print(class(l))
      #nan_indices <- which(l %in% NaN == TRUE)
      l[which(l %in% NaN == TRUE)] = 0
      
      #v_1 = (b1 %*% as.matrix(phase_I) %*% C_11 %*% t(as.matrix(phase_I)) %*% diag(l)) + (b2 %*% as.matrix(phase_II) %*% C_21 %*% t(as.matrix(phase_I)) %*% diag(l))
      v_1 = (b1 %*% phase_I %*% tcrossprod(C_11, phase_I * (l))) + (b2 %*% phase_II %*% tcrossprod(C_21, phase_I * (l)))
      #L1 <- diag(l)/study.list[[3]]
      
      
      e1 <- as.vector(1/(1 + exp(-rowSums(ref_dat * theta_CC))))
      e2 <- as.vector((exp(-rowSums(ref_dat * theta_CC)) - 1)/(exp(-rowSums(ref_dat * theta_CC)) + 1))
      e3 <- as.vector(1/(1 + exp(rowSums(ref_dat * theta_CC))))
      l <- e1*e2*e3
      #print(class(l))
      nan_indices <- which(l %in% NaN == TRUE)
      l[nan_indices] <- 0
      #print(class(b2))
      #v_2 = (b1 %*% phase_I %*% C_12 %*% t(phase_II) %*% (-diag(l))) + (b2 %*% phase_II %*% C_22 %*% t(phase_II) %*% (-diag(l)))
      v_2 = (b1 %*% phase_I %*% tcrossprod(C_12, phase_II * (-l))) + (b2 %*% phase_II %*% tcrossprod(C_22, phase_II * (-l)))
      #L2 <- -diag(l)/study.list[[4]]
      
      #v = c(v_1, v_2)
      #rm(v_1)
      #rm(v_2)
      #L = Matrix(magic::adiag(L1, L2), sparse=TRUE)
      
      #rm(L1)
      #rm(L2)
      rm(l)
      
      #b1 = ((as.vector(1/(1 + exp(-ref_dat %*% beta_old))) - as.vector(1/(1 + exp(-as.matrix(phase_I) %*% theta_R))))/sampling_fraction)/study.list[[3]]
      #b2 =  (study.list[[2]]$Y - as.vector(1/(1 + exp(-rowSums(ref_dat * theta_CC)))))/dim(phase_II)[1]
      #b = c(b1, b2)
      #start.time = Sys.time()
      #V2 = Matrix(diag(as.vector(t(b) %*% X_abdiag %*% C %*% t(X_abdiag) %*% L)), sparse = TRUE)
      #Sys.time() - start.time
      #start.time = Sys.time()
      #V = eigenMapMatMult(as.matrix(t(X_abdiag)), as.matrix(L))
      #V = eigenMapMatMult(as.matrix(C), V)
      #V = eigenMapMatMult(as.matrix(X_abdiag), V)
      #V = eigenMapMatMult(as.matrix(t(b)), V)
      #V = diag(as.vector(V))
      #Sys.time() - start.time
      #print(dim(Dn_1))
      #print(length(Dn_2))
      #Dn <- t(X_rbind) %*% W %*% X_abdiag %*% C %*% t(X_abdiag) %*% b
      
   
      
      Dn = crossprod((ref_dat * w_1),phase_I) %*% C_11 %*% crossprod(phase_I, b1) +
           crossprod((ref_dat * w_1),phase_I) %*% C_12 %*% crossprod(phase_II, b2) + 
           crossprod((ref_dat * w_2),phase_II) %*% C_21 %*% crossprod(phase_I, b1) +
           crossprod((ref_dat * w_2),phase_II) %*% C_22 %*% crossprod(phase_II, b2)
      
      # Dn = (t(ref_dat) %*% diag(w_1) %*% phase_I %*% C_11 %*% t(phase_I) %*% b1) +
      #      (t(ref_dat) %*% diag(w_1) %*% phase_I %*% C_12 %*% t(phase_II) %*% b2) +
      #      (t(ref_dat) %*% diag(w_2) %*% phase_II %*% C_21 %*% t(phase_I) %*% b1) +
      #      (t(ref_dat) %*% diag(w_2) %*% phase_II %*% C_22 %*% t(phase_II) %*% b2)
      #      
     
      #sqrt_C = sqrtm(C)
      #W_star_first <- W %*% X_abdiag %*% sqrtm(C) %*% t(X_abdiag) %*% W
      #W_star_first = Matrix(tcrossprod(W %*% X_abdiag %*% sqrtm(C)), sparse = TRUE)
      #print(class(W_star_first))
      #W_star <- W_star_first + V

      #rm(V)
      # print(W_star_first)
      # 
      # if(is.symmetric.matrix(W_star_first) != TRUE)
      #   print("W_star_first Not symmetric")
      # 
      # if(is.symmetric.matrix(W_star_first) != TRUE)
      #   print("W_star_second Not symmetric")
      # 
      # if(is.symmetric.matrix(W_star) != TRUE)
      #   print("W_star Not symmetric")
      # 
      # if(is.negative.definite(W_star, tol=1e-8) == TRUE)
      #   print("W_star Negative definite")
      #Define Jacobian here
      #J_n <- t(X_rbind) %*% W_star %*% X_rbind
      
      #J_n = eigenMapMatMult(eigenMapMatMult(as.matrix(t(X_rbind)), as.matrix(W_star)), as.matrix(X_rbind))
      J_n = crossprod((ref_dat * w_1), phase_I) %*% tcrossprod(C_11, crossprod((ref_dat*w_1), phase_I)) + 
      crossprod((ref_dat * w_1), phase_I) %*% C_12 %*% crossprod((phase_II * w_2), ref_dat) +
      crossprod((ref_dat * w_2), phase_II) %*% C_21 %*% crossprod((phase_I * w_1), ref_dat) + 
      crossprod((ref_dat * w_2), phase_II) %*% C_22 %*% crossprod((phase_II * w_2), ref_dat) +
      crossprod((ref_dat * (as.numeric(v_1) + as.numeric(v_2))), ref_dat)
      

        
        
        
      # J_n = t(ref_dat) %*% diag(w_1) %*% as.matrix(phase_I) %*% C_11 %*% t(t(ref_dat) %*% diag(w_1) %*% as.matrix(phase_I)) +
      #   t(ref_dat) %*% diag(w_1) %*% as.matrix(phase_I) %*% C_12 %*% t(as.matrix(phase_II)) %*% diag(w_2) %*% ref_dat +
      #   t(t(ref_dat) %*% diag(w_1) %*% as.matrix(phase_I) %*% C_12 %*% t(as.matrix(phase_II)) %*% diag(w_2) %*% ref_dat)  +
      #   t(ref_dat) %*% diag(w_2) %*% as.matrix(phase_II) %*% C_22 %*% t(t(ref_dat) %*% diag(w_2) %*% as.matrix(phase_II)) +
      #   (t(ref_dat) %*% (diag(as.numeric(v_1)) + diag(as.numeric(v_2))) %*% ref_dat)
      
      # if(is.symmetric.matrix(J_n) != TRUE)
      # {
      #   print("J_n Not symmetric")
      #   print(class(W_star))
      #   print(class(X_rbind))
      #   if(is.symmetric.matrix(W_star) != TRUE)
      #     print("W_star Not symmetric")
      # }
        
      #print(class(J_n))
      #print(J_n)
      #print(is.nan(J_n_beta))
      #print(det(J_n))
      if(det(J_n) == 0)
      { 
        beta_old <- rep(NA, ncol(ref_dat))
        status = 0
        #print(det(J_n))
        #print("The Jacobian is singular")
        break;
      }
      beta_new <- beta_old - qr.solve(J_n, Dn)
      #print(D_n_beta_t)
      #print(beta_old)
      #print(beta_new)
      eps_inner <- sqrt(sum((beta_new - beta_old)^2))
      #print(eps_inner)
      beta_old <- as.numeric(beta_new)
      iter = iter + 1
      #print("Number of iterations \n")
      #print(iter)
      if(eps_inner < threshold_optim || iter > 2000)
      {
         continue <- FALSE
         if(iter >= 2000 && eps_inner >= threshold_optim)
         {
           status = 0
         }
      }

    }
    no_of_iter_outer <- no_of_iter_outer + 1

    return(list("beta_optim" = beta_old,"iter_IRWLS" = iter - 1, "Status" = status, "no.of.iter.outer" = no_of_iter_outer))
  }

