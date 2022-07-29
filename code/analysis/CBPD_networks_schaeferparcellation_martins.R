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
subjdata_dir="~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/subjectData/"
outdir=paste0("~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/",pipeline)
cluster_outdir=paste0("/cbica/home/tooleyu/projects/in_progress/within_between_network_conn_CBPD/data/imageData/",pipeline)
netdata_dir=paste0("~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/",pipeline)
cluster_mounted_data="/cbica/home/tooleyu/projects/in_progress/within_between_network_conn_CBPD/data/imageData/"

#ages <- read.csv(paste0(subjdata_dir,"CBPD_data_DMD_2020.02.24_revised.csv")) #find most recent CBPD data wherever it is
main <- read.csv(paste0(subjdata_dir,"CBPD_data_DMD_2020.05.22_UNOFFICIAL_mergedPITMRIdropout.csv")) #find most recent CBPD data wherever it is
withavgweight <- read.csv(paste0(netdata_dir,"/n150_long_inc_within_between_Yeo7_avgruns_",parcellation,pipeline,"_withmodulpartcoef_with_QA.csv"))

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
#main <- main %>% filter(.,ID!="sub-CBPD0020")

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
main_excluded <- main %>% filter(., fd_mean_avg > 1  | pctSpikesFD_avg > 0.3)

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
cluster_mounted_data="~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/imageData/"

main_replicate <- read.csv(paste0(subjdata_dir,"CBPD_data_DMD_2020.05.22_UNOFFICIAL_mergedPITMRIdropout.csv")) #find most recent CBPD data wherever it is
withavgweight <- read.csv(paste0(netdata_dir,"/n150_long_inc_within_between_Yeo7_avgruns_",parcellation,pipeline,"_withmodulpartcoef_with_QA.csv"))

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

#filter out only kids under 8
main_replicate_under9 <- filter(main_replicate_filt, age_scan <=9)

#take only the first data we have from a subject
main_replicate_filt$base_ID <- stri_split(main_replicate_filt$record_id.y,fixed="_", simplify=T)[,1]
main_replicate_filt$longitudinal_visit_num[is.na(main_replicate_filt$longitudinal_visit_num)] <- 1

#only take one set of brain data from each participant, no matter the timepoint.
main_replicate_unique <- main_replicate_filt  %>% group_by(base_ID) %>% arrange(longitudinal) %>% filter(row_number() == 1)
dim(main_replicate_unique)

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
  outfile=paste0(outdir,"/CBPD_n91_schaefer400_avgruns.Rdata")
  outfile2=paste0(cluster_outdir,"/CBPD_n91_schaefer400_avgruns.Rdata")
}else {outfile=paste0(outdir,"/CBPD_n92_schaefer400_allruns.Rdata")
outfile2=paste0(cluster_outdir,"/CBPD_n92_schaefer400_allruns.Rdata")
}
save(main, main_filt, main_unique,file=outfile)
save(main, main_filt, main_unique,file=outfile2)

save(main, main_filt, main_unique, networks_age_pvals_fdr,networks_ses_pvals_fdr, networks_age_pvals_fdr_replicate,networks_ses_pvals_fdr_replicate, main_replicate, main_replicate_filt, main_replicate_unique, file=outfile)
save(main, main_filt, main_unique, networks_age_pvals_fdr,networks_ses_pvals_fdr, networks_age_pvals_fdr_replicate,networks_ses_pvals_fdr_replicate, main_replicate, main_replicate_filt, main_replicate_unique, file=outfile2)
}

