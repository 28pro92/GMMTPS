# GMMTPS
It contains all the codes for simulation and real data analysis for calculating the GMM estimator in two-phase studies which we call as GMMTPS estimator.Some of the R (version 4.1.1, Platform: x86_64-apple-darwin17.0 (64-bit), Running under: macOS Big Sur 10.16) packages the user needs to install include *survival*, *survey*, *OrdNor*, *magic*, *dplyr*, *ggplot2* and *ggpubr*.
## Simulation
*Simulation* folder contains all the codes associated with the simulation studies in the paper.
### Table 1 Results in the paper:
The Table1_results folder contains all the codes evaluating the performance of GMMTPS estimator and its comparison with SPML estimator.  
#### GMMTPS
- *simulation_NEW_X_2_S_inPhase_I_case_control.R* file is the main file evaluating the GMMTPS estimator. This file uses the *myoptim.R* file for the GMM (two-step) optimization. The user needs to run *simulation_NEW_X_2_S_inPhase_I_case_control.R* for reproducing the results in Table 1 for the GMMTPS estimator.
-  Make sure to change the path of the file *myoptim.R* accordingly that is sourced in *simulation_NEW_X_2_S_inPhase_I_case_control.R* file.
- The results after the user runs the *simulation_NEW_X_2_S_inPhase_I_case_control.R* file are stored in *GENMETA_case_control_table1.rds*. It contains 1000 rows (indicating 1000 simulations) and contains 27 columns (first 13 are the GMMTPS estimates for the parameters associated with the full model , 14 to 26 contain the variance of the estimates computed using asymptotic formula and 27th column indicates if the algorithm converged) 
#### SPML
- *missreg_comparision_simu_1.R* is the main file evaluating the SPML estimator using the missreg package by Chris Wild.
- *missreg_comparision_simu_1.R* uses five files which are sourced in the file itself. The five files are *bin2stg.R*, *binmods.R*, *utilities.R*, *MLEfn.R* and *MLInf.R*. These files are provided by Prof. Chris Wild.
- *missreg_case_control_table1.rds* contains the results for 1000 simulations evaluating the SPML estimator.
### Table 2 Results in the paper:
The *Table2_results* folder contains all the codes evaluating the performance of GMMTPS estimator and its comparison with SPML estimator for reproducing the results in Table 2 of the paper. The simulations are performed imitating the structure of *nwtco* data from National Wilm's Tumor Study. The *nwtco* data is available in the *survival* package in R.
#### Table2_results
- *genmeta_simulation_case_control.R* file is the main file evaluating the GMMTPS estimator in the case-control design setup. The phase-I model contains the main terms as indicated in the paper. The results are stored in *genmeta_realdata_simulation_CC_w_main_effects_phase_I.rds*.
- *genmeta_simulation_case_control_with_added_interaction.R* file is the main file evaluating the GMMTPS estimator in the case-control design setup. The phase-I model is a saturated model including interaction terms in addition to the main terms as indicated in the paper. The results are stored in *genmeta_realdata_simulation_CC_w_additional_interaction.rds*.
- *genmeta_simulation_balanced.R* file is the main file evaluating the GMMTPS estimator in the case-control design setup. The phase-I model contains the main terms as indicated in the paper. The results are stored in *genmeta_realdata_simulation_balanced_w_main_effects_phase_I.rds*.
- *genmeta_simulation_balanced_with_added_interaction.R* file is the main file evaluating the GMMTPS estimator in the case-control design setup. The phase-I model is a saturated model including interaction terms in addition to the main terms as indicated in the paper. The results are stored in *genmeta_realdata_simulation_balanced_w_additional_interaction.rds*. 
- *results_in_table_2.R* file calculates the relative effciency of GMMTPS estimator compared to the SPML estimator from the *.rds* files, ones from our method and the others from SPML method which are stored in *missreg3* folder described below.
##### missreg3
This sub-folder contains the codes for evaluating the SPML estimator according the same simulaton scheme.
- *real_data_simulation_genord_case_control_added_interaction.R* is the main file evaluating the SPML estimator in the case-control design setup. The results are stored in *missreg_real_data_simulation_genord_cc_added_interaction.rds*
- *real_data_simulation_genord_balanced_added_interaction.R* is the main file evaluating the SPML estimator in the balanced design setup. The results are stored in *missreg_real_data_simulation_genord_balanced_added_interaction.rds*
## Real Data Analysis
The *Real_data_analysis* folder contains the code for comparing the performance of three methods namely GMMTPS, SPML and calibration in the *nwtco* data from the National Wilms Tumor Study. *figures1and2.R* contains the code for reproducing figures 1 and 2 in the paper.
### GMMTPS
- *Genmeta_realdata_casecontrol_w_age.R* contains the code for evaluating the GMMTPS estimator in the *nwtco* data from the *survival* package under the case-control design discussed in Table 3 of the paper. The results are based on 1000 simualtions of phase-II data. They are stored in *real_data_GENMETA_w_age_case_control.rds*. The corresponding code file and result file for the balanced design are *Genmeta_realdata_balanced_w_age.R* and *real_data_GENMETA_w_age_balanced.rds*, respectively.
### SPML
*missreg* sub-folder contains the code files for estimating using SPML method under the same set-up as in GMMTPS.  
- *missreg_real_data_age_case_control.R* contains the code for evaluating the GMMTPS estimator in the *nwtco* data from the *survival* package under the case-control design discussed in Table 3 of the paper. The results are based on 1000 simualtions of phase-II data. They are stored in *missreg_realdata_case_control_w_age.rds*. The corresponding code file and result file for the balanced design are *missreg_real_data_age_balanced.R* and *missreg_realdata_balanced_w_age.rds*, respectively.
### Calibration
The following *.R* files contain the code for estimating using the calibration approache described in *Kulich et.al.* under the same set-up as in GMMTPS. *Raking* was used as the calibration function from the *survey* package.
- *calibrate_cc.R* contains the code for evaluating the calibration estimator in the *nwtco* data from the *survival* package under the case-control design discussed in Table 3 of the paper. The results are based on 1000 simualtions of phase-II data. They are stored in *calibrate_cc.rds*. The corresponding code file and result file for the balanced design are *calibrate_balanced.R* and *calibrate_balanced.rds*, respectively.

# References
