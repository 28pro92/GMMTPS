# GMMTPS
It contains all the codes for simulation and real data analysis for calculating the GMM estimator in two-phase studies which we call as GMMTPS estimator.
## Simulation
This folder contains all the codes associated with the simulation studies in the paper.
### Table 1 Results in the paper:
The Table1_results folder contains all the codes evaluating the performance of GMMTPS estimator and its comparison with SPML estimator.  
#### GMMTPS estimator
- *simulation_NEW_X_2_S_inPhase_I_case_control.R* file is the main file evaluating the GMMTPS estimator. This file uses the *myoptim.R* file for the GMM (two-step) optimization. The user needs to run *simulation_NEW_X_2_S_inPhase_I_case_control.R* for reproducing the results in Table 1 for the GMMTPS estimator. Make sure to change the path of the file *myoptim.R* accordingly that is sourced in *simulation_NEW_X_2_S_inPhase_I_case_control.R* file.
