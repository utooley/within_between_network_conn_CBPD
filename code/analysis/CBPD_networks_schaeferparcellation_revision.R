library(dplyr)
library(stats)
library(parallel)
library(lm.beta)
library(packrat)
library(summarytools)
library(psych)
library(polycor)
library(psy)
library(stringi)
#library(nFactors)
#library(MASS)
# Load Data ---------------------------------------------------------------
#pipelines=c('nogsr_spkreg_fd0.5dvars1.75_drpvls','gsr_spkreg_fd0.5dvars1.75_drpvls', "gsr_censor_5contig_fd0.5dvars1.75_drpvls", "gsr_censor_5contig_fd1.25dvars2_drpvls", "nogsr_spkreg_fd1.25dvars2_drpvls")
pipelines=c('nogsr_spkreg_fd0.5dvars1.75_drpvls', 'gsr_spkreg_fd0.5dvars1.75_drpvls', "nogsr_spkreg_fd1.25dvars2_drpvls")
#runs=c("run1", "run2")
run="both"
#for (run in runs){
for (pipeline in pipelines){
parcellation="schaefer400_"
pipeline= 'gsr_spkreg_fd0.25dvars1.75_drpvls'
subjdata_dir="~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/subjectData/"
outdir=paste0("~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/",pipeline)
cluster_outdir=paste0("~/cbica/home/tooleyu/projects/in_progress/within_between_network_conn_CBPD/data/imageData/",pipeline)
netdata_dir=paste0("~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/",pipeline)
cluster_mounted_data="~/cbica/home/tooleyu/projects/in_progress/within_between_network_conn_CBPD/data/imageData/"

#ages <- read.csv(paste0(subjdata_dir,"CBPD_data_DMD_2020.02.24_revised.csv")) #find most recent CBPD data wherever it is
main <- read.csv(paste0(subjdata_dir,"CBPD_data_DMD_2020.05.22_UNOFFICIAL_mergedPITMRIdropout.csv")) #find most recent CBPD data wherever it is
withavgweight <- read.csv(paste0(netdata_dir,"/n125_long_inc_within_between_Yeo7_avgruns_",parcellation,pipeline,"_withmodulpartcoef_with_QA.csv"))

#filter out some unneeded variables
#main <- main %>% dplyr::select(., -c(work_life_balance_questionnaire_timestamp:les_c_other_explain, colorado_child_temperament_index_timestamp:ccti_sum, cbcl_18mo_admin:cbcl_total_t))

#filter out extra variables in average weight
withavgweight$record_id <- as.character(withavgweight$ID)
withavgweight$ID <- as.character(withavgweight$ID)
#merge in main CBPD data
main$ID <- paste0("sub-",as.character(main$record_id))
main$ID <- sub("_","",main$ID)
main <- merge(withavgweight,main, by="ID")

#take out the participant with the glitter in her hair and artifact in rest?
main <- main %>% filter(.,ID!="sub-CBPD0020")

#Take out anyone with a 1 in the exclude column
main <- main %>% filter(.,exclude!=1)

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

main$avg_parentedu <- main %>% select(c("parent1_edu","parent2_edu")) %>% rowMeans(., na.rm=T) #average parent education first
main$ses_composite <- main %>% select(income_median, avg_parentedu) %>% scale(.) %>% rowMeans(., na.rm = T)

#make categorical child aces
main$aces3category <- ifelse(main$childaces_sum_ignorenan == 0, 0, ifelse(main$childaces_sum_ignorenan == 1, 1, ifelse(main$childaces_sum_ignorenan==2, 2, ifelse(main$childaces_sum_ignorenan>=3, 3, NA))))
main$aces3category <- factor(main$aces3category, labels=c("None"," or One", "Two", "Three+"))

#make a num uncensored vols variable
main$totalUncensoredVols <- main$totalSizet-main$totalnVolCensored

#Take out people with mean motion over 1, or more than 50% of frames censored (pctSpikesFD), or max motion > 10 mm
main_filt <- main %>% filter(., pctVolsCensored< 0.3 & relMaxRMSMotion < 10 & fd_mean_avg < 1 & totalSizet > 130)
length(unique(main_filt$ID))
main_filt <- main %>% filter(., relMaxRMSMotion < 10 & fd_mean_avg < 1 & fd_perc_2mm_avg < 10 & totalSizet > 130)
main_filt <- main %>% filter(., pctVolsCensored< 0.5 & relMaxRMSMotion < 10 & fd_mean_avg < 1 & fd_perc_2mm_avg < 10 & totalSizet > 130 & totalUncensoredVols > 100)

#criteria allyson suggests, not using different FD thresholds for censoring and for the % of volumes over a threshold that warrants subject exclusion
main_filt <- main %>% filter(., fd_mean_avg < 1  & pctSpikesFD_avg < 0.3)
length(unique(main_filt$ID)) #n=64 unique with new pipeline

#for new preprocessing, don't filter out based on spikes?
main_filt <- main %>% filter(., fd_mean_avg < 0.5)
#n=94 unique

#or filter only when > .50 have spikes?
main_filt <- main %>% filter(., fd_mean_avg < 1  & pctSpikesFD_avg < 0.5)
#n=85 unique

#or filter only if < 110 volumes remain?
main_filt <- main %>% filter(., fd_mean_avg < 1  &  totalUncensoredVols > 100)
#n=82 unique 

#criteria reviewer suggests, using 0.5mm FD threshold and for the % of volumes over a threshold that warrants subject exclusion
#main_filt <- main %>% filter(., fd_mean_avg < 0.5  & pctSpikesFD_avg < 0.3)
length(unique(main_filt$ID))

#filter out only kids under 8
main_under9 <- filter(main_filt, age_scan <=9)

#take only the first data we have from a subject
main_filt$base_ID <- stri_split(main_filt$record_id.y,fixed="_", simplify=T)[,1]
main_filt$longitudinal_visit_num[is.na(main_filt$longitudinal_visit_num)] <- 1

#only take one set of brain data from each participant, no matter the timepoint.
main_unique <- main_filt  %>% group_by(base_ID) %>% arrange(longitudinal) %>% filter(row_number() == 1)
dim(main_unique)

#Check differences between sample included and excluded 
#remember to keep this person in the exclude sample, the participant with the glitter in her hair and artifact in rest
main <- main %>% filter(.,ID!="sub-CBPD0020") 
#criteria allyson suggests, not using different FD thresholds for censoring and for the % of volumes over a threshold that warrants subject exclusion
main_excluded <- main %>% filter(., fd_mean_avg > 1  | pctSpikesFD_avg > 0.3 | ID == "sub-CBPD0020")

#take only the first data we have from a subject
main_excluded$base_ID <- stri_split(main_excluded$record_id.y,fixed="_", simplify=T)[,1]
main_excluded$longitudinal_visit_num[is.na(main_excluded$longitudinal_visit_num)] <- 1

#only take one set of brain data from each participant, no matter the timepoint.
main_exclude <- main_excluded  %>% group_by(base_ID) %>% arrange(longitudinal) %>% filter(row_number() == 1)
dim(main_exclude)
t.test(main_exclude$age_scan, main_unique$age_scanh); wilcox.test(main_exclude$age_scan, main_unique$age_scan)

main_unique$matrix_reasoning_both <- ifelse(is.na(main_unique$wppsi_matrix_valid),main_unique$wisc_matrix_scaled,ifelse(main_unique$wppsi_matrix_valid==0,NA,main_unique$wppsi_matrix_scaled))
main_unique$matrix_reasoning_both_raw <- ifelse(is.na(main_unique$wppsi_matrix_valid),main_unique$wisc_matrix_raw,ifelse(main_unique$wppsi_matrix_valid==0,NA,main_unique$wppsi_matrix_raw))
main_unique$matrix_reasoning_both;main_unique$matrix_reasoning_both_raw
main_exclude$matrix_reasoning_both <- ifelse(is.na(main_exclude$wppsi_matrix_valid),main_exclude$wisc_matrix_scaled,ifelse(main_exclude$wppsi_matrix_valid==0,NA,main_exclude$wppsi_matrix_scaled))
main_exclude$matrix_reasoning_both_raw <- ifelse(is.na(main_exclude$wppsi_matrix_valid),main_exclude$wisc_matrix_raw,ifelse(main_exclude$wppsi_matrix_valid==0,NA,main_exclude$wppsi_matrix_raw))
main_exclude$matrix_reasoning_both;main_exclude$matrix_reasoning_both_raw

t.test(main_exclude$matrix_reasoning_both, main_unique$matrix_reasoning_both); wilcox.test(main_exclude$matrix_reasoning_both, main_unique$matrix_reasoning_both)
# Load Replicate in Schaefer200 Data ---------------------------------------------------------------
parcellation="schaefer200_"
subjdata_dir="~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/subjectData/"
outdir=paste0("~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/",pipeline)
netdata_dir=paste0("~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/",pipeline)
cluster_mounted_data="~/cbica/home/tooleyu/projects/in_progress/within_between_network_conn_CBPD/data/imageData/"

main_replicate <- read.csv(paste0(subjdata_dir,"CBPD_data_DMD_2020.05.22_UNOFFICIAL_mergedPITMRIdropout.csv")) #find most recent CBPD data wherever it is
withavgweight <- read.csv(paste0(netdata_dir,"/n125_long_inc_within_between_Yeo7_avgruns_",parcellation,pipeline,"_withmodulpartcoef_with_QA.csv"))

#filter out extra variables in average weight
withavgweight$record_id <- withavgweight$ID
#merge in main_replicate CBPD data
main_replicate$ID <- paste0("sub-",main_replicate$record_id)
main_replicate$ID <- sub("_","",main_replicate$ID)
main_replicate <- merge(withavgweight,main_replicate, by="ID")

#take out the participant with the glitter in her hair and artifact in rest?
main_replicate <- main_replicate %>% filter(.,ID!="sub-CBPD0020")

#Take out anyone with a 1 in the exclude column
main_replicate <- main_replicate %>% filter(.,exclude!=1)

# Data Cleaning -----------------------------------------------------------
#make factor variables as factors
main_replicate$male <- factor(main_replicate$male, levels=c(0,1), labels=c("Female", "Male"))
main_replicate$part_coef <- (main_replicate$part_coef_neg+main_replicate$part_coef_pos)/2

#make race simpler
main_replicate$race2 <- ifelse(main_replicate$race_americanindian==1, 3, ifelse(main_replicate$race_asian == 1, 3, ifelse(main_replicate$race_hawaiian==1, 3, ifelse(main_replicate$race_other==1, 3,(ifelse(main_replicate$race_black==1, 2, ifelse(main_replicate$race_white==1, 1, ifelse(main_replicate$race_americanindian+main_replicate$race_asian+main_replicate$race_black+main_replicate$race_hawaiian+main_replicate$race_white > 1, 3, NA))))))))
main_replicate$race2 <- factor(main_replicate$race2, labels=c("White", "Black", "Other"))
main_replicate$ethnicity <- factor(main_replicate$ethnicity, labels=c("Not Hispanic or Latino","Hispanic or Latino"))
#filter out some unneeded variables
main_replicate <- main_replicate %>% dplyr::select(., -c(colorado_child_temperament_index_timestamp:ccti_sum))
#main_replicate <- main_replicate %>% select(.,-c(cbcl_18mo_admin:cbcl_6yr_complete))
#main_replicate <- main_replicate %>% select(.,-c(has_diagnoses:letterword_identification_comple))

main_replicate$avg_parentedu <- main_replicate %>% select(c("parent1_edu","parent2_edu")) %>% rowMeans(., na.rm=T) #average parent education first
main_replicate$ses_composite <- main_replicate %>% select(income_median, avg_parentedu) %>% scale(.) %>% rowMeans(., na.rm = T)

#make categorical child aces
main_replicate$aces3category <- ifelse(main_replicate$childaces_sum_ignorenan <= 1, 0, ifelse(main_replicate$childaces_sum_ignorenan == 2, 1, ifelse(main_replicate$childaces_sum_ignorenan>=3, 2, NA)))
main_replicate$aces3category <- factor(main_replicate$aces3category, labels=c("None or One", "Two", "Three+"))

#make a num uncensored vols variable
main_replicate$totalUncensoredVols <- main_replicate$totalSizet-main_replicate$totalnVolCensored

#Take out people with mean motion over 1, or more than 50% of frames censored (pctSpikesFD), or max motion > 10 mm
main_replicate_filt <- main_replicate %>% filter(., pctVolsCensored< 0.5 & relMaxRMSMotion < 10 & fd_mean_avg < 1 & fd_perc_2mm_avg < 10) #old threshold
main_replicate_filt <- main_replicate %>% filter(., pctVolsCensored< 0.5 & relMaxRMSMotion < 10 & fd_mean_avg < 1 & fd_perc_2mm_avg < 10 & totalSizet > 130)
main_replicate_filt <- main_replicate %>% filter(., pctVolsCensored< 0.5 & relMaxRMSMotion < 10 & fd_mean_avg < 1 & fd_perc_2mm_avg < 10 & totalSizet > 130 & totalUncensoredVols > 100)

#criteria allyson suggests, not using different FD thresholds for censoring and for the % of volumes over a threshold that warrants subject exclusion
main_replicate_filt <- main_replicate %>% filter(., fd_mean_avg < 1 & pctSpikesFD_avg < 0.3)

# #criteria reviewer suggests, using 0.5mm FD threshold and for the % of volumes over a threshold that warrants subject exclusion
# main_replicate_filt <- main_replicate %>% filter(., fd_mean_avg < 0.5  & pctSpikesFD_avg < 0.3)
# length(unique(main_replicate_filt$ID))
# 
# #for new preprocessing, don't filter out based on spikes?
main_replicate_filt <- main_replicate %>% filter(., fd_mean_avg < 0.5 )
# #n=94 unique
# 
# #or filter only when > .50 have spikes?
main_replicate_filt <- main_replicate %>% filter(., fd_mean_avg < 1  & pctSpikesFD_avg < 0.5)
# #n=85 unique
# 
# #or filter only if < 100 volumes remain?
main_replicate_filt <- main_replicate %>% filter(., fd_mean_avg < 1  &  totalUncensoredVols > 100)
# #n=82 unique 

#filter out only kids under 8
main_replicate_under9 <- filter(main_replicate_filt, age_scan <=9)

#take only the first data we have from a subject
main_replicate_filt$base_ID <- stri_split(main_replicate_filt$record_id.y,fixed="_", simplify=T)[,1]
main_replicate_filt$longitudinal_visit_num[is.na(main_replicate_filt$longitudinal_visit_num)] <- 1

#only take one set of brain data from each participant, no matter the timepoint.
main_replicate_unique <- main_replicate_filt  %>% group_by(base_ID) %>% arrange(longitudinal) %>% filter(row_number() == 1)
dim(main_replicate_unique)

# Age and networks --------------------------------------------------------

# Look at each column and correct for multiple comparisons
covariates="~ age_scan+male+fd_mean_avg+avgweight+totalSizet"

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

#make a dataframe with no repeats of net comparisons
main_replicate_unique <- dplyr::select(main_replicate_unique, -c(sys2to1,sys3to1,sys3to2,sys4to1,sys4to2,sys4to3,sys5to1,sys5to2,sys5to3,sys5to4,sys6to1,sys6to2,sys6to3,sys6to4,sys6to5,sys7to1,sys7to2,sys7to3,sys7to4,sys7to5,sys7to6))
#run and compare for multiple comparisons again
m <- mclapply(names((dplyr::select(ungroup(main_unique),sys1to1:sys7to7))), function(sys) {as.formula(paste(sys, covariates, sep=""))},mc.cores=2)
networks_Age_pvals <- mclapply(m, function(sys) { summary(lm(formula = sys,data=main_replicate_unique))$coefficients[2,4]},mc.cores=1)
networks_Age_pvals <- as.data.frame(networks_Age_pvals)
networks_Age_pvals <- t(networks_Age_pvals)
networks_Age_pvals <- as.data.frame(networks_Age_pvals)
#bonferroni correct
networks_Age_pvals_fdr_replicate <- p.adjust(networks_Age_pvals$V1, method="fdr")
networks_age_pvals_fdr_replicate <- data.frame(networks_Age_pvals_fdr_replicate,names((dplyr::select(ungroup(main_unique),sys1to1:sys7to7))))
#FDR correction shows DMN to attentional networks and visual to dorsal attention is marginal.
colnames(networks_age_pvals_fdr_replicate) <- c("pvalue", "network")
print(networks_age_pvals_fdr_replicate)

# Separate Networks and Effects of SES ------------------------------------
covariates="~ age_scan+male+fd_mean_avg+avgweight+totalSizet+ses_composite+race2+ethnicity"
#make a dataframe with no repeats of net comparisons
#main_unique <- dplyr::select(main_unique, -c(sys2to1,sys3to1,sys3to2,sys4to1,sys4to2,sys4to3,sys5to1,sys5to2,sys5to3,sys5to4,sys6to1,sys6to2,sys6to3,sys6to4,sys6to5,sys7to1,sys7to2,sys7to3,sys7to4,sys7to5,sys7to6))
#run and compare for multiple comparisons again
m <- mclapply(names((dplyr::select(ungroup(main_unique),sys1to1:sys7to7))), function(sys) {as.formula(paste(sys, covariates, sep=""))},mc.cores=2)
networks_ses_pvals <- mclapply(m, function(sys) { summary(lm(formula = sys,data=main_unique))$coefficients[7,4]},mc.cores=1)
networks_ses_pvals <- as.data.frame(networks_ses_pvals)
networks_ses_pvals <- t(networks_ses_pvals)
networks_ses_pvals <- as.data.frame(networks_ses_pvals)
#bonferroni correct
networks_ses_pvals_fdr <- p.adjust(networks_ses_pvals$V1, method="fdr")
networks_ses_pvals_fdr <- data.frame(networks_ses_pvals_fdr,names((dplyr::select(ungroup(main_unique),sys1to1:sys7to7))))
colnames(networks_ses_pvals_fdr) <- c("pvalue", "network")
networks_ses_pvals_fdr

#make a dataframe with no repeats of net comparisons
#main_replicate_unique <- dplyr::select(main_replicate_unique, -c(sys2to1,sys3to1,sys3to2,sys4to1,sys4to2,sys4to3,sys5to1,sys5to2,sys5to3,sys5to4,sys6to1,sys6to2,sys6to3,sys6to4,sys6to5,sys7to1,sys7to2,sys7to3,sys7to4,sys7to5,sys7to6))
#run and compare for multiple comparisons again
m <- mclapply(names((dplyr::select(ungroup(main_unique),sys1to1:sys7to7))), function(sys) {as.formula(paste(sys, covariates, sep=""))},mc.cores=2)
networks_ses_pvals <- mclapply(m, function(sys) { summary(lm(formula = sys,data=main_replicate_unique))$coefficients[7,4]},mc.cores=1)
networks_ses_pvals <- as.data.frame(networks_ses_pvals)
networks_ses_pvals <- t(networks_ses_pvals)
networks_ses_pvals <- as.data.frame(networks_ses_pvals)
#bonferroni correct
networks_ses_pvals_fdr_replicate <- p.adjust(networks_ses_pvals$V1, method="fdr")
networks_ses_pvals_fdr_replicate <- data.frame(networks_ses_pvals_fdr_replicate,names((dplyr::select(ungroup(main_unique),sys1to1:sys7to7))))
colnames(networks_ses_pvals_fdr_replicate) <- c("pvalue", "network")
networks_ses_pvals_fdr_replicate

# Save models for use in markdown file ------------------------------------
if(run=='run1'){
  outfile=paste0(outdir,"/CBPD_n91_schaefer400_run1.Rdata")
} else if (run=='run2'){
  outfile=paste0(outdir,"/CBPD_n91_schaefer400_run2.Rdata")
} else if (run=='average'){
  outfile=paste0(outdir,"/CBPD_n64_schaefer400_avgrunsFD1outliers30.Rdata")
  outfile2=paste0(cluster_outdir,"/CBPD_n64_schaefer400_avgrunsFD1outliers30.Rdata")
}else {outfile=paste0(outdir,"/CBPD_n83_schaefer400_allruns_FD0.5mmcensored0.2mmnooutliers.Rdata")
outfile2=paste0(cluster_outdir,"/CBPD_n83_schaefer400_allrunsFD0.5mmcensored0.2mmnooutliers.Rdata")
}
save(main, main_filt, main_unique,file=outfile)
save(main, main_filt, main_unique,file=outfile2)

save(main, main_filt, main_unique, networks_age_pvals_fdr,networks_ses_pvals_fdr, networks_age_pvals_fdr_replicate,networks_ses_pvals_fdr_replicate, main_replicate, main_replicate_filt, main_replicate_unique, file=outfile)
save(main, main_filt, main_unique, networks_age_pvals_fdr,networks_ses_pvals_fdr, networks_age_pvals_fdr_replicate,networks_ses_pvals_fdr_replicate, main_replicate, main_replicate_filt, main_replicate_unique, file=outfile2)
}

# Age and measures of segregation -----------------------------------------
#Need to control for the amount of good data a participant has, size_t

#look at mean within and between with age
lm_within_sys_age <- lm(mean_within_sys~age_scan+male+fd_mean_avg+avgweight+totalSizet, data=main_unique)
summary(lm_within_sys_age)
lm.beta(lm_within_sys_age)
lm.beta(l)
lm_between_sys_age <- lm(mean_between_sys~age_scan+male+fd_mean_avg+avgweight+totalSizet, data=main_unique)
summary(lm_between_sys_age)
lm.beta(lm_between_sys_age)
lm_segreg_age<- lm(system_segreg~age_scan+male+fd_mean_avg+avgweight+pctVolsCensored+totalSizet, data=main_unique)
summary(lm_segreg_age)
lm.beta(lm_segreg_age)
lm_modul_age <- lm(modul_avg~age_scan+male+fd_mean_avg+avgweight+pctVolsCensored+totalSizet, data=main_unique)
summary(lm_modul_age)
lm.beta(lm_modul_age)
lm_part_coef_age <- lm(part_coef~age_scan+male+fd_mean_avg+avgweight+pctVolsCensored+totalSizet, data=main_unique)
summary(lm_part_coef_age)
lm.beta(lm_part_coef_age)

#Does number of communities detected with modul change with age?
lm_num_comms_age <- lm(num_comms_modul_avg~age_scan+male+fd_mean_avg+avgweight+pctVolsCensored+totalSizet, data=main_unique)
summary(lm_num_comms_age) #Decreases slightly with age
#look at average weight with age
l <- lm(avgweight~age_scan+male+fd_mean_avg+pctVolsCensored+totalSizet, data=main_unique)
summary(l)
lm.beta(l)

# Longitudinal random effects models --------------------------------------
library(nlme)
library(lme4)
#melt data
main_filt_long <- dplyr::select(main_filt, c(ID:record_id.y,age_scan,male,fd_mean_avg,avgweight,pctVolsCensored,totalSizet,race2,ethnicity, longitudinal_visit_num,base_ID))
#need to filter out only the first datapoint from each timepoint
main_filt_long <- main_filt_long  %>% group_by(ID) %>% filter(row_number() == 1)
main_long_only <- main_filt_long %>% group_by(base_ID) %>% filter(any(longitudinal_visit_num==2))

library(lattice)
xyplot(mean_within_sys ~ age_scan | base_ID, main_long_only, aspect = "xy",
            type = c("g", "p", "r"), index.cond = function(x,y) coef(lm(y ~ x))[1])

#look at mean within and between with age
lm_within_sys_age <- lmer(sys4to7~scale(age_scan)+male+scale(fd_mean_avg)+scale(avgweight)+scale(pctVolsCensored)+scale(totalSizet)+(1|base_ID), data=main_filt)
summary(lm_within_sys_age)


lm.beta(lm_within_sys_age)
~age_scan+male+fd_mean+avgweight+pctVolsCensored+totalSizet
lm.beta(l)
lm_between_sys_age <- lm(mean_between_sys~age_scan+male+fd_mean+avgweight+pctVolsCensored+totalSizet, data=main_unique)
summary(lm_between_sys_age)
lm.beta(l)
lm_segreg_age<- lm(system_segreg~age_scan+male+fd_mean+avgweight+pctVolsCensored+totalSizet, data=main_unique)
summary(lm_segreg_age)
lm_modul_age <- lm(modul_avg~age_scan+male+fd_mean+avgweight+pctVolsCensored+totalSizet, data=main_unique)
summary(lm_modul_age)
lm_part_coef_age <- lm(part_coef~age_scan+male+fd_mean+avgweight+pctVolsCensored+totalSizet, data=main_unique)
summary(lm_part_coef_age)

 # Non-linear effects of age? ----------------------------------------------
gam_part_coef_age <- gam(part_coef~s(age_scan)+male+fd_mean+avgweight+pctSpikesFD+size_t, data=main_unique)
l <- gam(modul~s(ageAtScan1cent)+sex+avgweight+envSES, data=master, method = "REML")

#This is in the RMarkdown document.

# Separate Networks and Effects of Age ------------------------------------
networks <- dplyr::select(main_unique, sys1to1:sys7to7) %>% ungroup()
nets <- colnames(networks);nets <- nets[-1]
#look at non-linear interaction between age and envSES 
for (net in nets){
  name<-paste0("lm_",net)
  formula<-formula(paste0(net, '~age_scan+male+fd_mean_avg+avgweight+totalSizet'))
  assign(name, lm(formula, data=main_unique))
  #p_val[net] <- summary(name)$coefficients[2,4]
}
summary(lm_sys1to1) #increase with age
lm.beta(lm_sys1to1)
summary(lm_sys1to2)
summary(lm_sys1to3) #increase with age
lm.beta(lm_sys1to3)
summary(lm_sys1to4)
summary(lm_sys1to5)
summary(lm_sys1to6)
summary(lm_sys1to7)
summary(lm_sys2to2) #increase
lm.beta(lm_sys2to2)
summary(lm_sys2to3)
summary(lm_sys2to4)
summary(lm_sys2to5)
summary(lm_sys2to6) #effect of SES composite
summary(lm_sys2to7)
summary(lm_sys3to3)
lm.beta(lm_sys3to3)
summary(lm_sys4to4)
lm.beta(lm_sys4to4)
summary(lm_sys5to5)
lm.beta(lm_sys5to5)
summary(lm_sys6to6)
lm.beta(lm_sys6to6)
summary(lm_sys3to6)
summary(lm_sys4to6)
summary(lm_sys5to6)
summary(lm_sys6to6)
summary(lm_sys6to7)
summary(lm_sys7to7)#increase with age
lm.beta(lm_sys7to7)
summary(lm_sys3to7) #marginal decrease with age
lm.beta(lm_sys3to7)
summary(lm_sys4to7)
lm.beta(lm_sys4to7)
visreg(lm_sys4to7) #strong decrease with age
summary(lm_sys7to5)
main_unique$within_system <- (main_unique$sys1to1+main_unique$sys2to2+main_unique$sys3to3
+main_unique$sys4to4+main_unique$sys5to5+main_unique$sys6to6+main_unique$sys7to7)/7
main


l <- lm(litnum_freq_art_avg~age_scan+ses_composite, data=main_unique)
networks <- dplyr::select(main_unique, litnum_freq_lit_avg:litnum_freq_finemotor_avg)
nets <- colnames(networks)
#look at non-linear interaction between age and envSES 
for (net in nets){
  name<-paste0("lm_",net)
  formula<-formula(paste0(net, '~age_behav*ses_composite'))
  assign(name, lm(formula, data=main_unique))
  #p_val[net] <- summary(name)$coefficients[2,4]
}
summary(lm_litnum_freq_art_avg)
summary(lm_litnum_freq_finemotor_avg)
visreg(lm_litnum_freq_finemotor_avg,"age_behav", by="ses_composite")
summary(lm_litnum_freq_pretend_avg)
summary(lm_litnum_freq_spatial_avg)
summary(lm_litnum_freq_lit_avg)
visreg(lm_litnum_freq_lit_avg,"age_behav", by="ses_composite")
summary(lm_litnum_freq_music_avg)

# Environmental effects on networks -------------------------------------------------

#look at segreg measures with SES
lm_within_sys_age <- lm(mean_within_sys~age_scan+male+fd_mean+avgweight+pctSpikesFD+size_t+ses_composite+race2+ethnicity, data=main_replicate_unique)
summary(lm_within_sys_age)
lm.beta(l)
lm_between_sys_age <- lm(mean_between_sys~age_scan+male+fd_mean+avgweight+pctSpikesFD+size_t+ses_composite+race2+ethnicity, data=main_replicate_unique)
summary(lm_between_sys_age)
lm.beta(l)
lm_segreg_age<- lm(system_segreg~age_scan+male+fd_mean+avgweight+pctSpikesFD+size_t+ses_composite+race2+ethnicity, data=main_replicate_unique)
summary(lm_segreg_age)
lm_modul_age <- lm(modul_avg~age_scan+male+fd_mean+avgweight+pctSpikesFD+size_t+ses_composite+race2+ethnicity, data=main_replicate_unique)
summary(lm_modul_age)
lm_part_coef_age <- lm(part_coef~age_scan+male+fd_mean+avgweight+pctSpikesFD+size_t+ses_composite+race2+ethnicity, data=main_replicate_unique)
summary(lm_part_coef_age)

visreg(lm_between_sys_age)
visreg(lm_within_sys_age)
visreg(lm_segreg_age)
visreg(lm_modul_age)
visreg(lm_part_coef_age)
visreg(lm_num_comms_age)

#look at whether there is an interaction
lm_within_sys_age_income <- lm(mean_within_sys~age_scan*ses_composite+race2+ethnicity+male+fd_mean+avgweight+pctSpikesFD+size_t, data=main_replicate_unique)
summary(lm_within_sys_age_income)
lm.beta(l)
lm_between_sys_age_income <- lm(mean_between_sys~age_scan*ses_composite+race2+ethnicity+male+fd_mean+avgweight+pctSpikesFD+size_t, data=main_replicate_unique)
summary(lm_between_sys_age_income)
lm.beta(l)
lm_segreg_age<- lm(system_segreg~age_scan*ses_composite+race2+ethnicity+male+fd_mean+avgweight+pctSpikesFD+size_t, data=main_replicate_unique)
summary(lm_segreg_age)
lm_modul_age <- lm(modul_avg~age_scan*ses_composite+race2+ethnicity+male+fd_mean+avgweight+pctSpikesFD+size_t, data=main_replicate_unique)
summary(lm_modul_age)
lm_part_coef_age <- lm(part_coef~age_scan*ses_composite+race2+ethnicity+male+fd_mean+avgweight+pctSpikesFD+size_t, data=main_replicate_unique)
summary(lm_part_coef_age)

visreg(lm_within_sys_age, "age_scan", by="income_median")
visreg(lm_between_sys_age, "age_scan", by="income_median")

# Specific Networks and Environment ---------------------------------------
nets=c("sys1to3", "sys3to7", "sys4to7")
for (net in nets){
  name<-paste0("lm_",net)
  formula<-formula(paste0(net, '~age_scan+male+fd_mean+avgweight+pctSpikesFD+size_t+ses_composite+race2+ethnicity'))
  assign(name, lm(formula, data=main_replicate_unique))
  #p_val[net] <- summary(name)$coefficients[2,4]
}
summary(lm_sys1to3)
summary(lm_sys3to7)
summary(lm_sys4to7)
visreg(lm_sys1to3)
visreg(lm_sys3to7)

for (net in nets){
  name<-paste0("lm_",net)
  formula<-formula(paste0(net, '~age_scan*ses_composite+male+fd_mean+avgweight+pctSpikesFD+size_t+race2+ethnicity'))
  assign(name, lm(formula, data=main_replicate_unique))
  #p_val[net] <- summary(name)$coefficients[2,4]
}
summary(lm_sys1to3)
summary(lm_sys3to7)
summary(lm_sys4to7)


# Parcel Specificity of Age and Age x SES Effects -------------------------

#read in file with edges by subject
edge_weights <- read.csv(paste0(cluster_mounted_data,"zedges_for_each_subj_081219.csv"))
edge_weights <- as.data.frame(edge_weights)
#join to master with covariates
edge_weights<-dplyr::rename(edge_weights, ID=Var1, run=Var2)
edge_weights<-right_join(edge_weights, main, by =c("ID", "run"))

#run it on all edges
covariates=" ~ age_scan+male+fd_mean+avgweight+pctSpikesFD+size_t"
#put the betas/pvals back into a matrix
m <- mclapply(names(edge_weights[,3:79803]), function(x) {as.formula(paste( x, covariates, sep=""))},mc.cores=4)
edge_pvals_age <- mclapply(m, function(x) { summary(lm(formula = x,data=edge_weights))$coefficients[2,4]},mc.cores=4)
edge_pvals_age <- as.data.frame(edge_pvals_age)
edge_pvals_age <- t(edge_pvals_age)  
edge_pvals_age <- as.data.frame(edge_pvals_age)
colnames(edge_pvals_age) <- "edge_pvals_age"
edge_pvals_age$Node_index <- 1:79800

#then, pull out only DMN to DAN edges to FDR correct 

#write out pvals into a matrix
my_matrix<-matrix(0,400,400)
my_matrix[upper.tri(my_matrix, diag=FALSE)]<-edge_pvals_age$edge_pvals_age
#write out a file to read in brainnetviewer
write.csv(my_matrix,paste0(cluster_mounted_data,"edge_pvals_mat_age_effect_081219.csv"))
write.csv(edge_pvals_age,paste0(cluster_mounted_data,"edge_pvals_raw_age_effect_081219.csv"))


#FDR CORRECT the pvals
m <- mclapply(names(edge_weights[,3:64263]), function(x) {as.formula(paste(x, covariates, sep=""))},mc.cores=2)
edge_pvals_agexses <- mclapply(m, function(x) { summary(lm(formula = x,data=edge_weights))$coefficients[8,4]},mc.cores=2)
edge_pvals_agexses <- as.data.frame(edge_pvals_agexses)
edge_pvals_agexses <- t(edge_pvals_agexses)
edge_pvals_agexses <- as.data.frame(edge_pvals_agexses)
colnames(edge_pvals_agexses) <- "edge_pvals_agexses"
edge_pvals_agexses$Node_index <- 1:64261
fdr_corrected<-p.adjust(edge_pvals_agexses$edge_pvals_agexses, method = "fdr")
sig_edges_AgexSES<-cbind(edge_pvals_agexses, fdr_corrected)
fdr_sig_edges_AgexSES<-sig_edges_AgexSES[sig_edges_AgexSES$fdr_corrected < 0.05,]
#write out pvals into a matrix
my_matrix<-matrix(0,359,359)
my_matrix[upper.tri(my_matrix, diag=FALSE)]<-edge_betas_agexses$edge_betas_agexses
#write out a file to read in brainnetviewer
write.csv(my_matrix,paste0(clustcodir,"edge_pvals_agexses_int.csv"))

#write out the betas by edge back into a matrix
my_matrix<-matrix(0,359,359)
my_matrix[upper.tri(my_matrix, diag=FALSE)]<-edge_betas_agexses$edge_betas_agexses
#write out a file to read in brainnetviewer
write.csv(my_matrix,paste0(clustcodir,"edge_betas_agexses_int_scaled.csv"))

#read in the matrix after it's made
my_matrix <- read.csv(paste0(clustcodir,"edge_betas_agexses_int.csv"))
my_matrix <- my_matrix[,-1]
#look at whether they're higher within significant nodes or to the rest of the brain?
#read in significant nodes for agexses
library(tibble)
fdr_sig_nodes_lm_AgexSES <- read.csv(paste0(clustcodir,"nodewise_pvals_for_GAMs_and_lms.csv"))
fdr_sig_nodes_only<-fdr_sig_nodes_lm_AgexSES[fdr_sig_nodes_lm_AgexSES$fdr_pvals_lm_agexses < 0.05,]
dim(fdr_sig_nodes_only)
colnames(my_matrix)=c(1:359)
my_matrix_edges <- as.data.frame(my_matrix)
my_matrix_edges <- rownames_to_column(my_matrix_edges, var="Node_Index")
sig_node_betas_only <- my_matrix_edges %>% select(one_of(as.character(fdr_sig_nodes_only$X))) %>% filter(my_matrix_edges$Node_Index %in% fdr_sig_nodes_only$X)
#get the mean of betas within only the significant nodes
mean(sig_node_betas_only[sig_node_betas_only!=0])
#get the mean of betas outside the significant nodes
fdr_non_sig_nodes<-fdr_sig_nodes_lm_AgexSES[fdr_sig_nodes_lm_AgexSES$fdr_pvals_lm_agexses > 0.05,]
non_sig_node_betas_only <- my_matrix_edges %>% select(one_of(as.character(fdr_non_sig_nodes$X))) %>% filter(my_matrix_edges$Node_Index %in% fdr_non_sig_nodes$X)
mean(non_sig_node_betas_only[non_sig_node_betas_only!=0])

#put the yeo edges as column and index names, sort by them and calculate the average beta within each network
yeomapping <- read.csv(paste0(clustcodir,"Glasser_to_Yeo_no52.csv"))
names(my_matrix_edges)=c(yeomapping$Yeo_Parcellation.7YeoPNC.nii.gz_0...)
my_matrix_edges$Yeo_mapping<-yeomapping$Yeo_Parcellation.7YeoPNC.nii.gz_0...

test<-my_matrix_edges[,1:359] 
test>0.10

#Yeo1
Yeo_1_edges <- my_matrix_edges[my_matrix_edges$Yeo_mapping== 1,]
Yeo_1_edges<- Yeo_1_edges[names(Yeo_1_edges)==1]
mean(Yeo_1_edges!=0)
#Yeo2
Yeo_2_edges <- my_matrix_edges[my_matrix_edges$Yeo_mapping== 2,]
Yeo_2_edges<- Yeo_2_edges[names(Yeo_2_edges)==2]
mean(Yeo_2_edges!=0)
#Yeo1
Yeo_3_edges <- my_matrix_edges[my_matrix_edges$Yeo_mapping== 3,]
Yeo_3_edges<- Yeo_3_edges[names(Yeo_3_edges)==3]
mean(Yeo_3_edges!=0)
#Yeo4
Yeo_4_edges <- my_matrix_edges[my_matrix_edges$Yeo_mapping== 4,]
Yeo_4_edges<- Yeo_4_edges[names(Yeo_4_edges)==4]
mean(Yeo_4_edges!=0)
#Yeo5
Yeo_5_edges <- my_matrix_edges[my_matrix_edges$Yeo_mapping== 5,]
Yeo_5_edges<- Yeo_5_edges[names(Yeo_5_edges)==5]
mean(Yeo_5_edges!=0)
#Yeo6
Yeo_6_edges <- my_matrix_edges[my_matrix_edges$Yeo_mapping== 6,]
Yeo_6_edges<- Yeo_6_edges[names(Yeo_6_edges)==6]
mean(Yeo_6_edges!=0)
#Yeo7
Yeo_7_edges <- my_matrix_edges[my_matrix_edges$Yeo_mapping== 7,]
Yeo_7_edges<- Yeo_7_edges[names(Yeo_7_edges)==7]
mean(Yeo_7_edges!=0)


# Regress on Each Edge-Distance Dependence --------------------------------
#read in the pvals
my_matrix <- read.csv(paste0(clustcodir,"edge_betas_agexses_int.csv"))
my_matrix <- my_matrix[,-1]
#read in the euclidean distance matrix of distances
distances <- read.table("~/Documents/bassett_lab/tooleyEnviNetworks/parcels/Glasser/GlasserEuclideanDistanceMatrix.txt")
#remove parcel 52
distances <- select(distances, -V103)
distances <- distances[-103,]
#select each pval by its corresponding distance, and then plot
distances
#vectorize my_matrix and vectorize distances, plot one by the other.
edge_pvals_agexses_vector <- my_matrix[upper.tri(my_matrix, diag=FALSE)]
distances_vector <- distances[upper.tri(distances, diag=FALSE)]
#plot the edge weight pvals by distance
plot(distances_vector, sig_edges_AgexSES$edge_pvals_agexses, col=alpha("blue",0.2), cex=0.2)

# Exploratory Factor Analysis of Stress -----------------------------------
#Child ACES, PSS, WLBQ items
#per Allyson on 8-9, used ACES sum score, PSS sum score, and WLBQ individual items for stress, along with neighborhood safety 
#sum the two neighborhood safety scores
stress_variables <- main_unique %>% dplyr::select(.,ID,pss_sum, childaces_sum_ignorenan, wlbq_q1:wlbq_q12,neighborhood_q1,neighborhood_q2) %>% data.frame()
#stress_variables <- main %>% dplyr::select(.,childaces_q1:childaces_q10b_notreversed, pss_q1:pss_q10,wlbq_q1:wlbq_q12)
#remove variables with no variance
#stress_variables <- stress_variables%>% dplyr::select(.,-c(childaces_q8,childaces_q9,childaces_q6b,childaces_q6a,childaces_q6d,  pss_q4_notreversed,pss_q5_notreversed,pss_q8_notreversed))
#stress_variables <- stress_variables%>% dplyr::select(.,-c(pss_q4_notreversed,pss_q5_notreversed,pss_q8_notreversed, pss_q7_notreversed))
stress_variables <- stress_variables[complete.cases(stress_variables),]
stress_variables_id <- stress_variables$ID
stress_variables <- dplyr::select(stress_variables, -ID)
view(dfSummary(stress_variables))
#since this is a mix of continuous and dichotomous variables, need to get a correlation matrix in some other way
#convert variables with < 4 categories to factor? Easier if they're continuous!
# col_names <- sapply(stress_variables, function(col) length(unique(col)) <= 4)
# stress_variables[ , col_names] <- lapply(stress_variables[ , col_names] , factor)

#one factor analysis with factanal()
fit <- factanal(stress_variables, factors = 2, rotation = "varimax", scores = "regression")
fit
scree.plot(fit$correlation)
#a different scree plot
ev <- eigen(cor(stress_variables)) # get eigenvalues
ap <- parallel(subject=nrow(stress_variables),var=ncol(stress_variables), rep=100, cent=.05)
nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)
plotnScree(nS) #looks like should be extracting 3 or 4 factors
firststressfactor <- fit$scores[,2]
temp <- data.frame(stress_variables_id,firststressfactor)
colnames(temp) <- c("ID","firststressfactor")
main_unique<- right_join(temp,main_unique,by="ID")
main_replicate_unique<- right_join(temp,main_replicate_unique,by="ID")

#another factor analysis with fa()
#main_unique <- main_unique %>% group_by(ID) %>% filter(row_number() == 1)
fa.2 <- fa(stress_variables, nfactors = 2, rotate = "varimax")
print(fa.2, cut = 0.1)
factorscores <- factor.scores(stress_variables, fa.2)$scores[,2]
temp <- data.frame(stress_variables_id,factorscores)
colnames(temp) <- c("ID","firststressfactor1")
main_unique<- right_join(temp,main_unique,by="ID")
main_replicate_unique<- right_join(temp,main_replicate_unique,by="ID")

#trying to do factor analysis with dichotomous variables, ignore for now.
#mixedCor(data=stress_variables,d=1:9,p=10:17, c=18:29) #for some reason it doesn't like pss_q5_notreversed or pss_q8_notreversed, take those out-not working
cormat <- polycor::hetcor(stress_variables)
fit <- factanal(covmat = cormat$correlations, factors = 3, rotation = "varimax")
fa.2 <- fa(r = cormat$correlations, nfactors = 3, n.obs = nrow(stress_variables))

# Exploratory Factor Analysis of Cognitive Stimulation -----------------------------------
#HOME, Literacy and Numeracy Questionnaire for cog stimulation, take only those from litnum that are applicable
#recode HOME into subscales
## RECODE HOME INTO SUBSCALES per that found online on "Early Childhood HOME Record Form"
main_unique$home_learningmats <- main_unique$home_q1+ main_unique$home_q2 + main_unique$home_q3 +main_unique$home_q4 + main_unique$home_q6+main_unique$home_q7 + main_unique$home_q8 +main_unique$home_q9 + main_unique$home_q10 + main_unique$home_q11+ main_unique$home_q15
main_unique$home_langstim <- main_unique$home_q12 + main_unique$home_q16 + main_unique$home_q17 + main_unique$home_q23+ main_unique$home_q24 
main_unique$home_academicstim <- main_unique$home_q18+main_unique$home_q19+main_unique$home_q20+main_unique$home_q21+main_unique$home_q22
main_unique$home_variety <- main_unique$home_q13+main_unique$home_q30+main_unique$home_q32+main_unique$home_q33+main_unique$home_q31+main_unique$home_q14+main_unique$home_q5
main_unique$home_modeling  #this requires coding of written-in responses
main_unique$home_total <- main_unique$home_learningmats +main_unique$home_langstim+main_unique$home_academicstim+main_unique$home_variety

#the HOME has only dichotomous variables, need to use IRT on this?
cogstim_variables <- main_unique %>% dplyr::select(.,ID,home_learningmats:home_variety, litnum_freq_num_avg:litnum_freq_finemotor_avg, litnum_child_books,litnum_hrs_childreadto) %>% data.frame()
#cogstim_variables <- main %>% dplyr::select(.,childaces_q1:childaces_q10b_notreversed, pss_q1:pss_q10,wlbq_q1:wlbq_q12)
#remove variables with no variance
#cogstim_variables <- cogstim_variables%>% dplyr::select(.,-c(childaces_q8,childaces_q9,childaces_q6b,childaces_q6a,childaces_q6d,  pss_q4_notreversed,pss_q5_notreversed,pss_q8_notreversed))
#cogstim_variables <- cogstim_variables%>% dplyr::select(.,-c(pss_q4_notreversed,pss_q5_notreversed,pss_q8_notreversed, pss_q7_notreversed))
cogstim_variables <- cogstim_variables[complete.cases(cogstim_variables),]
cogstim_variables_id <- cogstim_variables$ID
cogstim_variables <- dplyr::select(cogstim_variables, -ID)
view(dfSummary(cogstim_variables))
#since this is a mix of continuous and dichotomous variables, need to get a correlation matrix in some other way
#convert variables with < 4 categories to factor? Easier if they're continuous!
# col_names <- sapply(cogstim_variables, function(col) length(unique(col)) <= 4)
# cogstim_variables[ , col_names] <- lapply(cogstim_variables[ , col_names] , factor)

#one factor analysis with factanal()
fit <- factanal(cogstim_variables, factors = 3, rotation = "varimax", scores = "regression")
scree.plot(fit$correlation)
#a different scree plot
ev <- eigen(cor(cogstim_variables)) # get eigenvalues
ap <- parallel(subject=nrow(cogstim_variables),var=ncol(cogstim_variables), rep=100, cent=.05)
nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)
plotnScree(nS) #looks like should be extracting 3 or 4 factors
firstcogstimfactor <- fit$scores[,1]
temp <- data.frame(cogstim_variables_id,firstcogstimfactor)
colnames(temp) <- c("ID","firstcogstimfactor")
main_unique<- right_join(temp,main_unique,by="ID")
main_replicate_unique<- right_join(temp,main_replicate_unique,by="ID")

#another factor analysis with fa()
#main_unique <- main_unique %>% group_by(ID) %>% filter(row_number() == 1)
fa.2 <- fa(cogstim_variables, nfactors = 3, rotate = "varimax")
print(fa.2,cut = .1)
factorscores <- factor.scores(cogstim_variables, fa.2)$scores[,1]
temp <- data.frame(cogstim_variables_id,factorscores)
colnames(temp) <- c("ID","firstcogstimfactor1")
main_unique<- right_join(temp,main_unique,by="ID")
main_replicate_unique<- right_join(temp,main_replicate_unique,by="ID")