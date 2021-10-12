library(survival)
library(survey)
cohort.study = nwtco
cohort.study$age = as.numeric(scale(cohort.study$age/12))

cohort.study$histol = factor(cohort.study$histol)
cohort.study$stage <- factor(cohort.study$stage)
beta_true = glm(rel~histol*stage + histol*age + stage*age , data = cohort.study, family = quasibinomial())$coefficients
results.npmle.survey = readRDS("/Users/prosenjitkundu/Desktop/Code_biometrics_github/Real_data_analysis/calibrate_cc.rds")
results.npmle.survey = results.npmle.survey[, c(1:8, 12, 9:11)]
sim.results.GENMETA.cc = readRDS("/Users/prosenjitkundu/Desktop/Code_biometrics_github/Real_data_analysis/real_data_GENMETA_w_age_case_contol.rds")
sim.results.missreg.cc = readRDS("/Users/prosenjitkundu/Desktop/Code_biometrics_github/Real_data_analysis/missreg/missreg_realdata_case_control_w_age.rds")
na.index = c(which(results.npmle.survey[, 8]> 3), which(is.na(results.npmle.survey[,1]) == T), which(is.na(sim.results.GENMETA.cc[,1])==T), which(is.na(sim.results.missreg.cc[,1])==T))
effect.simulations = 1000 - length(unique(na.index))
beta_true_matrix = matrix(rep(beta_true[-1], each = effect.simulations), nrow=effect.simulations)
MSE.cc.calibrate = apply((results.npmle.survey[-na.index, 1:12]-beta_true_matrix)^2, 2, mean)
genmeta_CC = apply((sim.results.GENMETA.cc[-na.index, 2:13]-beta_true_matrix)^2, 2, mean)
missreg_CC = apply((sim.results.missreg.cc[-na.index, 2:13]-beta_true_matrix)^2, 2, mean)

Method = rep(c("GMM", "SPML", "Calibration"), 12)
var.names = rep(c("UH","Stage-II","Stage-III","Stage-IV","Age","UHxStage-II","UHxStage-III", "UHxStage-IV", "UHxAge", "Stage-IIxAge", "Stage-IIIxAge", "Stage-IVxAge"), each=3)

data.plot.mse.cc = data.frame(cbind(c(rbind(genmeta_CC, missreg_CC, MSE.cc.calibrate)), var.names, Method))
colnames(data.plot.mse.cc) = c("MSE", "Variable.Name", "Method")
data.plot.mse.cc$MSE = as.numeric(data.plot.mse.cc$MSE)

p1 = ggplot(data.plot.mse.cc[1:3,c(1,3)], aes(x = Method, MSE), ylab = "") + geom_bar(stat = "identity", width = .4) + scale_y_continuous(limits=c(0,max(data.plot.mse.cc[1:3,c(1,3)][,1]))) + xlab("UH") 
p2 = ggplot(data.plot.mse.cc[4:6,c(1,3)], aes(Method, MSE), ylab = "") + geom_bar(stat = "identity", width = .4) + scale_y_continuous(limits=c(0,max(data.plot.mse.cc[4:6,c(1,3)][,1]))) + xlab("Stage-II") 
p3 = ggplot(data.plot.mse.cc[7:9,c(1,3)], aes(Method, MSE), ylab = "") + geom_bar(stat = "identity", width = .4) + scale_y_continuous(limits=c(0,max(data.plot.mse.cc[7:9,c(1,3)][,1]))) + xlab("Stage-III") 
p4 = ggplot(data.plot.mse.cc[10:12,c(1,3)], aes(Method, MSE), ylab = "") + geom_bar(stat = "identity", width = .4) + scale_y_continuous(limits=c(0,max(data.plot.mse.cc[10:12,c(1,3)][,1]))) + xlab("Stage-IV") 
p5 = ggplot(data.plot.mse.cc[13:15,c(1,3)], aes(Method, MSE), ylab = "") + geom_bar(stat = "identity", width = .4) + scale_y_continuous(limits=c(0,max(data.plot.mse.cc[13:15,c(1,3)][,1]))) + xlab("Age") 
p6 = ggplot(data.plot.mse.cc[16:18,c(1,3)], aes(Method, MSE), ylab = "") + geom_bar(stat = "identity", width = .4) + scale_y_continuous(limits=c(0,max(data.plot.mse.cc[16:18,c(1,3)][,1]))) + xlab("UHxStage-II") 
p7 = ggplot(data.plot.mse.cc[19:21,c(1,3)], aes(Method, MSE), ylab = "") + geom_bar(stat = "identity", width = .4) + scale_y_continuous(limits=c(0,max(data.plot.mse.cc[19:21,c(1,3)][,1]))) + xlab("UHxStage-III") 
p8 = ggplot(data.plot.mse.cc[22:24,c(1,3)], aes(Method, MSE), ylab = "") + geom_bar(stat = "identity", width = .4) + scale_y_continuous(limits=c(0,max(data.plot.mse.cc[22:24,c(1,3)][,1]))) + xlab("UHxStage-IV")  
p9 = ggplot(data.plot.mse.cc[25:27,c(1,3)], aes(Method, MSE), ylab = "") + geom_bar(stat = "identity", width = .4) + scale_y_continuous(limits=c(0,max(data.plot.mse.cc[25:27,c(1,3)][,1]))) + xlab("UHxAge") 
p10 = ggplot(data.plot.mse.cc[28:30,c(1,3)], aes(Method, MSE), ylab = "") + geom_bar(stat = "identity", width = .4) + scale_y_continuous(limits=c(0,max(data.plot.mse.cc[28:30,c(1,3)][,1]))) + xlab("Stage-IIxAge")  
p11 = ggplot(data.plot.mse.cc[31:33,c(1,3)], aes(Method, MSE), ylab = "") + geom_bar(stat = "identity", width = .4) + scale_y_continuous(limits=c(0,max(data.plot.mse.cc[31:33,c(1,3)][,1]))) + xlab("Stage-IIIxAge") 
p12 = ggplot(data.plot.mse.cc[34:36,c(1,3)], aes(Method, MSE), ylab = "") + geom_bar(stat = "identity", width = .4) + scale_y_continuous(limits=c(0,max(data.plot.mse.cc[34:36,c(1,3)][,1]))) + xlab("Stage-IVxAge") 
mse.plot.cc = ggarrange(p1 + rremove("y.title"), p2 + rremove("y.title"), p3 + rremove("y.title"), p4 + rremove("y.title"), p5 + rremove("y.title"), p6 + rremove("y.title"), p7 + rremove("y.title"), p8 + rremove("y.title"), p9 + rremove("y.title"), p10 + rremove("y.title"), p11 + rremove("y.title"), p12 + rremove("y.title"), ncol = 4, nrow = 3)
#png("/Users/prosenjitkundu/Dropbox/Biometrics_two_phase/Resubmission1/CC_stage_instit_interaction_imp_model_revision.png", res = 1200)
annotate_figure(mse.plot.cc, left = text_grob("MSE", rot = 90, size=12))



#Balanced
cohort.study = nwtco
cohort.study$age = as.numeric(scale(cohort.study$age/12))

cohort.study$histol = factor(cohort.study$histol)
cohort.study$stage <- factor(cohort.study$stage)
beta_true = glm(rel~histol*stage + histol*age + stage*age , data = cohort.study, family = quasibinomial())$coefficients
results.npmle.survey = readRDS("/Users/prosenjitkundu/Desktop/Code_biometrics_github/Real_data_analysis/calibrate_balanced.rds")
results.npmle.survey = results.npmle.survey[, c(1:8, 12, 9:11)]
sim.results.GENMETA.balanced = readRDS("/Users/prosenjitkundu/Desktop/Code_biometrics_github/Real_data_analysis/real_data_GENMETA_w_age_balanced.rds")
sim.results.missreg.balanced = readRDS("/Users/prosenjitkundu/Desktop/Code_biometrics_github/Real_data_analysis/missreg/missreg_realdata_balanced_w_age.rds")
na.index = c(which(results.npmle.survey[, 8]> 3), which(is.na(results.npmle.survey[,1]) == T), which(is.na(sim.results.GENMETA.balanced[,1])==T), which(is.na(sim.results.missreg.balanced[,1])==T))
effect.simulations = 1000 - length(unique(na.index))
if(length(na.index) > 0)
{
  beta_true_matrix = matrix(rep(beta_true[-1], each = effect.simulations), nrow=effect.simulations)
  MSE.balanced.calibrate = apply((results.npmle.survey[-na.index, 1:12]-beta_true_matrix)^2, 2, mean)
  genmeta_balanced = apply((sim.results.GENMETA.balanced[-na.index, 2:13]-beta_true_matrix)^2, 2, mean)
  missreg_balanced = apply((sim.results.missreg.balanced[-na.index, 2:13]-beta_true_matrix)^2, 2, mean)
}else{
  beta_true_matrix = matrix(rep(beta_true[-1], each = effect.simulations), nrow=effect.simulations)
  MSE.balanced.calibrate = apply((results.npmle.survey[, 1:12]-beta_true_matrix)^2, 2, mean)
  genmeta_balanced = apply((sim.results.GENMETA.balanced[, 2:13]-beta_true_matrix)^2, 2, mean)
  missreg_balanced = apply((sim.results.missreg.balanced[, 2:13]-beta_true_matrix)^2, 2, mean)
}

Method = rep(c("GMM", "SPML", "Calibration"), 12)
var.names = rep(c("UH","Stage-II","Stage-III","Stage-IV","Age","UHxStage-II","UHxStage-III", "UHxStage-IV", "UHxAge", "Stage-IIxAge", "Stage-IIIxAge", "Stage-IVxAge"), each=3)

data.plot.mse.balanced = data.frame(cbind(c(rbind(genmeta_balanced, missreg_balanced, MSE.balanced.calibrate)), var.names, Method))
colnames(data.plot.mse.balanced) = c("MSE", "Variable.Name", "Method")
data.plot.mse.balanced$MSE = as.numeric(data.plot.mse.balanced$MSE)

p1 = ggplot(data.plot.mse.balanced[1:3,c(1,3)], aes(x = Method, MSE), ylab = "") + geom_bar(stat = "identity", width = .4) + scale_y_continuous(limits=c(0,max(data.plot.mse.balanced[1:3,c(1,3)][,1]))) + xlab("UH") 
p2 = ggplot(data.plot.mse.balanced[4:6,c(1,3)], aes(Method, MSE), ylab = "") + geom_bar(stat = "identity", width = .4) + scale_y_continuous(limits=c(0,max(data.plot.mse.balanced[4:6,c(1,3)][,1]))) + xlab("Stage-II") 
p3 = ggplot(data.plot.mse.balanced[7:9,c(1,3)], aes(Method, MSE), ylab = "") + geom_bar(stat = "identity", width = .4) + scale_y_continuous(limits=c(0,max(data.plot.mse.balanced[7:9,c(1,3)][,1]))) + xlab("Stage-III") 
p4 = ggplot(data.plot.mse.balanced[10:12,c(1,3)], aes(Method, MSE), ylab = "") + geom_bar(stat = "identity", width = .4) + scale_y_continuous(limits=c(0,max(data.plot.mse.balanced[10:12,c(1,3)][,1]))) + xlab("Stage-IV") 
p5 = ggplot(data.plot.mse.balanced[13:15,c(1,3)], aes(Method, MSE), ylab = "") + geom_bar(stat = "identity", width = .4) + scale_y_continuous(limits=c(0,max(data.plot.mse.balanced[13:15,c(1,3)][,1]))) + xlab("Age") 
p6 = ggplot(data.plot.mse.balanced[16:18,c(1,3)], aes(Method, MSE), ylab = "") + geom_bar(stat = "identity", width = .4) + scale_y_continuous(limits=c(0,max(data.plot.mse.balanced[16:18,c(1,3)][,1]))) + xlab("UHxStage-II") 
p7 = ggplot(data.plot.mse.balanced[19:21,c(1,3)], aes(Method, MSE), ylab = "") + geom_bar(stat = "identity", width = .4) + scale_y_continuous(limits=c(0,max(data.plot.mse.balanced[19:21,c(1,3)][,1]))) + xlab("UHxStage-III") 
p8 = ggplot(data.plot.mse.balanced[22:24,c(1,3)], aes(Method, MSE), ylab = "") + geom_bar(stat = "identity", width = .4) + scale_y_continuous(limits=c(0,max(data.plot.mse.balanced[22:24,c(1,3)][,1]))) + xlab("UHxStage-IV")  
p9 = ggplot(data.plot.mse.balanced[25:27,c(1,3)], aes(Method, MSE), ylab = "") + geom_bar(stat = "identity", width = .4) + scale_y_continuous(limits=c(0,max(data.plot.mse.balanced[25:27,c(1,3)][,1]))) + xlab("UHxAge") 
p10 = ggplot(data.plot.mse.balanced[28:30,c(1,3)], aes(Method, MSE), ylab = "") + geom_bar(stat = "identity", width = .4) + scale_y_continuous(limits=c(0,max(data.plot.mse.balanced[28:30,c(1,3)][,1]))) + xlab("Stage-IIxAge")  
p11 = ggplot(data.plot.mse.balanced[31:33,c(1,3)], aes(Method, MSE), ylab = "") + geom_bar(stat = "identity", width = .4) + scale_y_continuous(limits=c(0,max(data.plot.mse.balanced[31:33,c(1,3)][,1]))) + xlab("Stage-IIIxAge") 
p12 = ggplot(data.plot.mse.balanced[34:36,c(1,3)], aes(Method, MSE), ylab = "") + geom_bar(stat = "identity", width = .4) + scale_y_continuous(limits=c(0,max(data.plot.mse.balanced[34:36,c(1,3)][,1]))) + xlab("Stage-IVxAge") 
mse.plot.balanced = ggarrange(p1 + rremove("y.title"), p2 + rremove("y.title"), p3 + rremove("y.title"), p4 + rremove("y.title"), p5 + rremove("y.title"), p6 + rremove("y.title"), p7 + rremove("y.title"), p8 + rremove("y.title"), p9 + rremove("y.title"), p10 + rremove("y.title"), p11 + rremove("y.title"), p12 + rremove("y.title"), ncol = 4, nrow = 3)
#png("/Users/prosenjitkundu/Dropbox/Biometrics_two_phase/Resubmission1/CC_stage_instit_interaction_imp_model_revision.png", res = 1200)
annotate_figure(mse.plot.balanced, left = text_grob("MSE", rot = 90, size=12))
