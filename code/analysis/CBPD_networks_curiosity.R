library(dplyr)
library(stats)
library(parallel)
library(lm.beta)
library(packrat)
library(summarytools)
#library(psych)
#library(polycor)
#library(psy)
#library(nFactors)
#library(MASS)
library(stringi)
library(PerformanceAnalytics)
# Load Data ---------------------------------------------------------------
parcellation="schaefer400_"
#pipelines=c('nogsr_spkreg_fd0.5dvars1.75_drpvls','gsr_spkreg_fd0.5dvars1.75_drpvls', "gsr_censor_5contig_fd0.5dvars1.75_drpvls", "gsr_censor_5contig_fd1.25dvars2_drpvls", "nogsr_spkreg_fd1.25dvars2_drpvls")
pipelines=c('nogsr_spkreg_fd0.5dvars1.75_drpvls', 'gsr_spkreg_fd0.5dvars1.75_drpvls', "nogsr_spkreg_fd1.25dvars2_drpvls")
#runs=c("run1", "run2")
pipeline= 'gsr_spkreg_fd0.5dvars1.75_drpvls'
run="both"
#for (run in runs){
for (pipeline in pipelines){
subjdata_dir="~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/subjectData/"
outdir=paste0("~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/",pipeline)
cluster_outdir=paste0("/data/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/imageData/",pipeline)
netdata_dir=paste0("~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/",pipeline)
cluster_mounted_data="/data/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/imageData/"

#ages <- read.csv(paste0(subjdata_dir,"CBPD_data_DMD_2020.02.24_revised.csv")) #find most recent CBPD data wherever it is
main <- read.csv(paste0(subjdata_dir,"CBPD_data_DMD_2020.02.24_revised.csv")) #find most recent CBPD data wherever it is
withavgweight <- read.csv(paste0(netdata_dir,"/n125_long_inc_within_between_Yeo7_avgruns_",parcellation,pipeline,"_withmodulpartcoef_with_QA.csv"))

#filter out some unneeded variables
main <- main %>% dplyr::select(., -c(work_life_balance_questionnaire_timestamp:home_complete, colorado_child_temperament_index_timestamp:ccti_sum, cbcl_18mo_admin:cbcl_total_t,Date_coder1:Total_ddpi_parent_bookreading_average))

#filter out extra variables in average weight
withavgweight$record_id <- as.character(withavgweight$ID)
withavgweight$ID <- as.character(withavgweight$ID)
#merge in main CBPD data
main$ID <- paste0("sub-",as.character(main$record_id))
main$ID <- sub("_","",main$ID)
main <- merge(withavgweight,main, by="ID")
# #keep only ages, merge those in
# ages$ID <- paste0("sub-",ages$record_id)
# ages <- ages %>% dplyr::select(., ID,dob_entered:age_ques)
# main <- merge(main,ages, by="ID")

#take out the participant with the glitter in her hair and artifact in rest?
main <- main %>% filter(.,ID!="sub-CBPD0020")

# Data Cleaning -----------------------------------------------------------
#make factor variables as factors
main$male <- factor(main$male, levels=c(0,1), labels=c("Female", "Male"))
main$part_coef <- (main$part_coef_neg+main$part_coef_pos)/2

#make race simpler
main$race2 <- ifelse(main$race_americanindian==1, 3, ifelse(main$race_asian == 1, 3, ifelse(main$race_hawaiian==1, 3, ifelse(main$race_other==1, 3,(ifelse(main$race_black==1, 2, ifelse(main$race_white==1, 1, ifelse(main$race_americanindian+main$race_asian+main$race_black+main$race_hawaiian+main$race_white > 1, 3, NA))))))))
#Multiracial is other
main$race2 <- factor(main$race2, labels=c("White", "Black", "Other"))
main$ethnicity <- factor(main$ethnicity, labels=c("Not Hispanic or Latino","Hispanic or Latino"))
#main <- main %>% select(.,-c(cbcl_18mo_admin:cbcl_6yr_complete))
#main <- main %>% select(.,-c(has_diagnoses:letterword_identification_comple))

main$ses_composite <- as.numeric(scale(main$parent1_edu)+scale(main$income_median))

#make categorical child aces
main$aces3category <- ifelse(main$childaces_sum_ignorenan == 0, 0, ifelse(main$childaces_sum_ignorenan == 1, 1, ifelse(main$childaces_sum_ignorenan==2, 2, ifelse(main$childaces_sum_ignorenan>=3, 3, NA))))
main$aces3category <- factor(main$aces3category, labels=c("None"," or One", "Two", "Three+"))

#Take out people with mean motion over 1, or more than 50% of frames censored (pctSpikesFD), or max motion > 10 mm
main_filt <- main %>% filter(., pctVolsCensored< 0.5 & relMaxRMSMotion < 10 & fd_mean_avg < 1 & fd_perc_2mm_avg < 10)
#check that people still have enough data left
main_filt$nremainingvols<- main_filt$totalSizet-main_filt$totalnVolCensored
hist(main_filt$nremainingvols)
main_filt <- main_filt %>% filter(., nremainingvols > 120)

#main_filt <- main %>% filter(., pctSpikesFD< 0.5 & relMaxRMSMotion < 10 & relMeanRMSMotion< 1)
#filter out only kids under 8
main_under9 <- filter(main_filt, age_scan <=9)

#take only the first data we have from a subject
main_filt$base_ID <- stri_split(main_filt$record_id.y,fixed="_", simplify=T)[,1]
main_filt$longitudinal_visit_num[is.na(main_filt$longitudinal_visit_num)] <- 1

#either take run 1 or run 2 depending on the value of the variable
main_unique <- main_filt  %>% group_by(base_ID) %>% arrange(longitudinal) %>% filter(row_number() == 1) %>% ungroup

#then melt it into wide format and average across participants with more than one run?
#view(dfSummary(main_unique))

# I-D Type Curiosity ------------------------------------------------------
#add the two together for a sum
main_unique$epcuriosity_sum_both <- main_unique$epcuriosity_dtype_sum+main_unique$epcuriosity_itype_sum
# Look at each column and correct for multiple comparisons
measures=select(main_unique,one_of("age_scan","fd_mean_avg","avgweight","epcuriosity_dtype_sum","epcuriosity_itype_sum")) %>% ungroup()
#measures=select(main_unique,contains("epcuriosity"))

# Distributions
hist(main_unique$epcuriosity_dtype_sum)
hist(main_unique$epcuriosity_itype_sum)
chart.Correlation(measures)

#D-type
lm_between_sys_dtype<- lm(mean_between_sys~age_scan+ses_composite+male+fd_mean+avgweight+pctVolsCensored+totalSizet+epcuriosity_dtype_sum, data=main_unique)
summary(lm_between_sys_dtype)
lm_within_sys_dtype<- lm(mean_within_sys~age_scan+male+fd_mean+ses_composite+avgweight+pctVolsCensored+totalSizet+epcuriosity_dtype_sum, data=main_unique)
summary(lm_within_sys_dtype)
visreg(lm_within_sys_dtype)

#I-type
lm_between_sys_itype<- lm(mean_between_sys~age_scan+male+fd_mean+ses_composite+avgweight+pctVolsCensored+totalSizet+epcuriosity_itype_sum, data=main_unique)
summary(lm_between_sys_itype)
lm_within_sys_itype<- lm(mean_within_sys~age_scan+male+fd_mean+avgweight+ses_composite+pctVolsCensored+totalSizet+epcuriosity_itype_sum, data=main_unique)
summary(lm_within_sys_itype)

#Sum
lm_between_sys_itype<- lm(mean_between_sys~age_scan+male+fd_mean+avgweight+ses_composite+pctVolsCensored+totalSizet+epcuriosity_sum_both, data=main_unique)
summary(lm_between_sys_itype)
lm_within_sys_itype<- lm(mean_within_sys~age_scan+male+fd_mean+avgweight+ses_composite+pctVolsCensored+totalSizet+epcuriosity_sum_both, data=main_unique)
summary(lm_within_sys_itype)

#do a hacky cross-validation--modify this code
final_dat<-data.frame(matrix(ncol = 2, nrow = length(play_brain$record_id)))
colnames(final_dat) <- c("beta","sig")
for (i in 1:length(play_brain$record_id)){
  df <- play_brain[-i,]
  model<-lm(litnum_freq_play_avg~sys4to7_clean+age_behav+male+fd_mean_avg+avgweight+totalSizet+SES+race_white+race_asian+race_other_haw_ai+ethnicity, data=df)
  final_dat$beta[i]<-model$coefficients[2]
  final_dat$sig[i]<-summary(model)$coefficients[2,4]
}


#Across all nets I-type
covariates="~ epcuriosity_itype_sum+age_scan+male+fd_mean_avg+avgweight+pctSpikesFD_avg+totalSizet"

#make a dataframe with no repeats of net comparisons
main_unique <- dplyr::select(main_unique, -c(sys2to1,sys3to1,sys3to2,sys4to1,sys4to2,sys4to3,sys5to1,sys5to2,sys5to3,sys5to4,sys6to1,sys6to2,sys6to3,sys6to4,sys6to5,sys7to1,sys7to2,sys7to3,sys7to4,sys7to5,sys7to6))
#run and compare for multiple comparisons again
m <- mclapply(names((dplyr::select(ungroup(main_unique),sys1to1:sys7to7))), function(sys) {as.formula(paste(sys, covariates, sep=""))},mc.cores=2)
networks_Age_pvals <- mclapply(m, function(sys) { summary(lm(formula = sys,data=main_unique))$coefficients[2,4]},mc.cores=1)
networks_Age_pvals <- as.data.frame(networks_Age_pvals)
networks_Age_pvals <- t(networks_Age_pvals)
networks_Age_pvals <- as.data.frame(networks_Age_pvals)
#bonferroni correct
networks_Age_pvals_fdr <- p.adjust(networks_Age_pvals$V1, method="fdr")
networks_age_pvals_fdr <- data.frame(networks_Age_pvals_fdr,names((dplyr::select(ungroup(main_unique),sys1to1:sys7to7))))
colnames(networks_age_pvals_fdr) <- c("pvalue", "network")
#FDR correction shows DMN to attentional networks and visual to dorsal attention is marginal.
print(networks_age_pvals_fdr)

#Across all nets D-type
covariates="~ epcuriosity_dtype_sum+age_scan+male+fd_mean_avg+avgweight+pctSpikesFD_avg+totalSizet"

#make a dataframe with no repeats of net comparisons
main_unique <- dplyr::select(main_unique, -c(sys2to1,sys3to1,sys3to2,sys4to1,sys4to2,sys4to3,sys5to1,sys5to2,sys5to3,sys5to4,sys6to1,sys6to2,sys6to3,sys6to4,sys6to5,sys7to1,sys7to2,sys7to3,sys7to4,sys7to5,sys7to6))
#run and compare for multiple comparisons again
m <- mclapply(names((dplyr::select(ungroup(main_unique),sys1to1:sys7to7))), function(sys) {as.formula(paste(sys, covariates, sep=""))},mc.cores=2)
networks_Age_pvals <- mclapply(m, function(sys) { summary(lm(formula = sys,data=main_unique))$coefficients[2,4]},mc.cores=1)
networks_Age_pvals <- as.data.frame(networks_Age_pvals)
networks_Age_pvals <- t(networks_Age_pvals)
networks_Age_pvals <- as.data.frame(networks_Age_pvals)
#bonferroni correct
networks_Age_pvals_fdr <- p.adjust(networks_Age_pvals$V1, method="fdr")
networks_age_pvals_fdr <- data.frame(networks_Age_pvals_fdr,names((dplyr::select(ungroup(main_unique),sys1to1:sys7to7))))
colnames(networks_age_pvals_fdr) <- c("pvalue", "network")
#FDR correction shows DMN to attentional networks and visual to dorsal attention is marginal.
print(networks_age_pvals_fdr)

# Question-asking and networks --------------------------------------------
#merge in data
q_asking <- read.csv("~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/subjectData/question_asking_both_coders_parent_qs_and_longitudinal_042320.csv", header = T, skip=1)
summary(q_asking)

q_asking$Total_q_kid_allpcit_merge <- ifelse(is.na(q_asking$Total_q_kid_allpcit_average),q_asking$Total_q_kid_allpcit_coder2, q_asking$Total_q_kid_allpcit_average)
q_asking$Total_pq_parent_bookreading_merge <- ifelse(is.na(q_asking$Total_pq_parent_bookreading_average),q_asking$Total_pq_parent_bookreading_coder2, q_asking$Total_pq_parent_bookreading_average)
q_asking$Total_q_parent_bookreading_merge <- ifelse(is.na(q_asking$Total_q_parent_bookreading_average),q_asking$Total_q_parent_bookreading_coder2, q_asking$Total_q_parent_bookreading_average)

q_asking$ID <- sub("_","",q_asking$record_id)
q_asking$ID <- paste0("sub-",as.character(q_asking$ID))
main_unique <- left_join(main_unique,q_asking, by="ID")
main <- left_join(main,q_asking, by="ID")
x <- main %>% select(ID,avgweight,Total_q_kid_allpcit_merge:Total_q_parent_bookreading_merge) #did it work?
Total_pq_parent_bookreading_average

# Look at each column and correct for multiple comparisons
measures=select(main_unique,one_of("age_scan","ses_composite","fd_mean_avg","avgweight","Total_q_kid_allpcit_merge","Total_pq_parent_bookreading_merge", "Total_q_parent_bookreading_merge"))
#measures=select(main_unique,contains("epcuriosity"))

# Distributions of q-asking
hist(main_unique$Total_q_kid_allpcit_merge)
hist(main_unique$Total_pq_parent_bookreading_merge)
chart.Correlation(measures)
#how does it relate to age and parent report in the full sample
main$epcuriosity_sum_both <- main$epcuriosity_dtype_sum+main$epcuriosity_itype_sum
summary(lm(epcuriosity_itype_sum~age_behav.x+ses_composite+race2+Total_q_kid_allpcit_merge,data=main))
summary(lm(epcuriosity_dtype_sum~age_behav.x+ses_composite+race2+Total_q_kid_allpcit_merge,data=main))
summary(lm(epcuriosity_sum_both~age_behav.x+ses_composite+race2+Total_q_kid_allpcit_merge,data=main))
summary(lm(epcuriosity_sum_both~age_behav.x+ses_composite+race2+Total_q_kid_allpcit_merge,data=main))
summary(lm(epcuriosity_q9~age_behav.x+ses_composite+race2+Total_q_kid_allpcit_merge,data=main))

summary(lm(Total_q_kid_allpcit_merge~age_behav.x+ses_composite+race2,data=main))
summary(lm(Total_pq_parent_bookreading_merge~age_behav.x+ses_composite,data=main))
summary(lm(Total_q_parent_bookreading_merge~age_behav.x+ses_composite,data=main))

#Kid q's
lm_between_sys_kidq<- lm(Total_q_kid_allpcit_merge~age_scan+male+fd_mean+avgweight+pctVolsCensored+totalSizet+mean_between_sys, data=main_unique)
summary(lm_between_sys_kidq)
lm_within_sys_kidq<- lm(Total_q_kid_allpcit_merge~age_scan+male+fd_mean+avgweight+pctVolsCensored+totalSizet+mean_within_sys , data=main_unique)
summary(lm_within_sys_kidq)

visreg(lm_within_sys_kidq)
visreg(lm_between_sys_kidq)

#Across all nets kid q's
covariates="~ Total_q_kid_allpcit_merge+age_scan+male+fd_mean_avg+avgweight+pctSpikesFD_avg+totalSizet"

#make a dataframe with no repeats of net comparisons
main_unique <- dplyr::select(main_unique, -c(sys2to1,sys3to1,sys3to2,sys4to1,sys4to2,sys4to3,sys5to1,sys5to2,sys5to3,sys5to4,sys6to1,sys6to2,sys6to3,sys6to4,sys6to5,sys7to1,sys7to2,sys7to3,sys7to4,sys7to5,sys7to6))
#run and compare for multiple comparisons again
m <- mclapply(names((dplyr::select(ungroup(main_unique),sys1to1:sys7to7))), function(sys) {as.formula(paste(sys, covariates, sep=""))},mc.cores=2)
networks_Age_pvals <- mclapply(m, function(sys) { summary(lm(formula = sys,data=main_unique))$coefficients[2,4]},mc.cores=1)
networks_Age_pvals <- as.data.frame(networks_Age_pvals)
networks_Age_pvals <- t(networks_Age_pvals)
networks_Age_pvals <- as.data.frame(networks_Age_pvals)
#bonferroni correct
networks_Age_pvals_fdr <- p.adjust(networks_Age_pvals$V1, method="fdr")
networks_age_pvals_fdr <- data.frame(networks_Age_pvals_fdr,names((dplyr::select(ungroup(main_unique),sys1to1:sys7to7))))
colnames(networks_age_pvals_fdr) <- c("pvalue", "network")
#FDR correction shows DMN to attentional networks and visual to dorsal attention is marginal.
print(networks_age_pvals_fdr)
