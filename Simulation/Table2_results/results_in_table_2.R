#---Table 2 Results for the relative efficiency---#
# Case-control with main effects in phase-I
t1_1 = readRDS("/Users/prosenjitkundu/Desktop/Code_biometrics_github/Simulation/Table2_results/genmeta_realdata_simulation_CC_w_main_effects_phase_I.rds")
t1_2 = readRDS("/Users/prosenjitkundu/Desktop/Code_biometrics_github/Simulation/Table2_results/missreg3/missreg_real_data_simulation_genord_cc_added_interaction.rds")

apply(t1_2[, 1:13], 2, var)/apply(t1_1[, 1:13], 2, var)

# Case-control with added interaction in phase-I
t1_1 = readRDS("/Users/prosenjitkundu/Desktop/Code_biometrics_github/Simulation/Table2_results/genmeta_realdata_simulation_CC_w_additional_interaction.rds")
apply(t1_2[, 1:13], 2, var)/apply(t1_1[, 1:13], 2, var)



# Balance with main effects in phase-I
t1_1 = readRDS("/Users/prosenjitkundu/Desktop/Code_biometrics_github/Simulation/Table2_results/genmeta_realdata_simulation_balanced_w_main_effects_phase_I.rds")
t1_2 = readRDS("/Users/prosenjitkundu/Desktop/Code_biometrics_github/Simulation/Table2_results/missreg3/missreg_real_data_simulation_genord_balanced_added_interaction.rds")
apply(t1_2[, 1:13], 2, var)/apply(t1_1[, 1:13], 2, var)

# Balance with added interaction in phase-I
t1_1 = readRDS("/Users/prosenjitkundu/Desktop/Code_biometrics_github/Simulation/Table2_results/genmeta_realdata_simulation_balanced_w_additional_interaction.rds")
apply(t1_2[, 1:13], 2, var)/apply(t1_1[, 1:13], 2, var)
