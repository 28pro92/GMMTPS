# GMMTPS
It contains all the codes for simulation and real data analysis for calculating the GMM estimator in two-phase studies which we call as GMMTPS estimator.Some of the R (version 4.1.1) packages the user needs to install include *survival*, *survey*, *OrdNor*, *magic*, *dplyr*, *ggplot2* and *ggpubr*.
## Simulation
This folder contains all the codes associated with the simulation studies in the paper.
### Table 1 Results in the paper:
The Table1_results folder contains all the codes evaluating the performance of GMMTPS estimator and its comparison with SPML estimator.  
#### GMMTPS estimator
- *simulation_NEW_X_2_S_inPhase_I_case_control.R* file is the main file evaluating the GMMTPS estimator. This file uses the *myoptim.R* file for the GMM (two-step) optimization. The user needs to run *simulation_NEW_X_2_S_inPhase_I_case_control.R* for reproducing the results in Table 1 for the GMMTPS estimator.
-  Make sure to change the path of the file *myoptim.R* accordingly that is sourced in *simulation_NEW_X_2_S_inPhase_I_case_control.R* file.
- The results after the user runs the *simulation_NEW_X_2_S_inPhase_I_case_control.R* file are stored in *GENMETA_case_control_table1.rds*. It contains 1000 rows (indicating 1000 simulations) and contains 27 columns (first 13 are the GMMTPS estimates for the parameters associated with the full model , 14 to 26 contain the variance of the estimates computed using asymptotic formula and 27th column indicates if the algorithm converged) 
#### SPML estimator
- *missreg_comparision_simu_1.R* is the main file evaluating the SPML estimator using the missreg package by Chris Wild.
- *missreg_comparision_simu_1.R* uses five files which are sourced in the file itself. The five files are *bin2stg.R*, *binmods.R*, *utilities.R*, *MLEfn.R* and *MLInf.R*. These files are provided by Prof. Chris Wild.
- *missreg_case_control_table1.rds* contains the results for 1000 simulations evaluating the SPML estimator.
### Table 2 Results in the paper:
## Real Data Analysis
