t2_1 = readRDS("/Users/prosenjitkundu/Desktop/Code_biometrics_github/Simulation/Table1_results/GENMETA_case_control_table1.rds")
t2_2 = readRDS("/Users/prosenjitkundu/Desktop/Code_biometrics_github/Simulation/Table1_results/missreg_case_control_table1.rds")

beta.true = c(-3, 0, log(1.3), log(1.3),  log(1.3))
beta.true.matrix = t(matrix(rep(beta.true[-2], 1000), 4, 1000))
apply((t2_1[, 1:4]-beta.true.matrix)^2, 2, mean)
CI.L = t2_1[,1:4] - 1.96*sqrt(t2_1[,5:8])
CI.U = t2_1[,1:4] + 1.96*sqrt(t2_1[,5:8])
length(which(CI.L[,2] <= beta.true[3] & CI.U[,2] >= beta.true[3]))
length(which(CI.L[,3] <= beta.true[3] & CI.U[,3] >= beta.true[3]))
length(which(CI.L[,4] <= beta.true[3] & CI.U[,4] >= beta.true[3]))
apply(CI.U-CI.L, 2, mean)

apply(t2_2[, 1:4], 2, var)/apply(t2_1[, 1:4], 2, var)

