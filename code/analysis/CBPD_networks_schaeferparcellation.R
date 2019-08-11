library(dplyr)
library(stats)
library(parallel)
library(lm.beta)
library(packrat)
library(summarytools)
library(psych)
library(polycor)
library(psy)
library(nFactors)
#library(MASS)
# Load Data ---------------------------------------------------------------
subjdata_dir="~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/subjectData/"
outdir="~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_censor_5contig_fd0.5dvars1.75_drpvls/"
netdata_dir="~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_censor_5contig_fd0.5dvars1.75_drpvls/"

ages <- read.csv(paste0(subjdata_dir,"CBPD_data_190729_age_ses_ccti.csv")) #find most recent CBPD data wherever it is
main <- read.csv(paste0(subjdata_dir,"CBPD_parent_questionnaires_upto_177_080819.csv")) #find most recent CBPD data wherever it is
withavgweight <- read.csv(paste0(netdata_dir,"n75_within_between_Yeo7_Schaefer400_gsr_censor_5contig_fd0.5dvars1.75_drpvls_withmodulpartcoef_with_QA.csv"))

#filter out extra variables in average weight
withavgweight$record_id <- withavgweight$ID
#merge in main CBPD data
main$ID <- paste0("sub-",main$record_id)
main <- merge(withavgweight,main, by="ID")
#keep only ages, merge those in
ages$ID <- paste0("sub-",ages$record_id)
ages <- ages %>% dplyr::select(., ID,dob_entered:age_ques)
main <- merge(main,ages, by="ID")

#take out the participant with the glitter in her hair and artifact in rest?
main <- main %>% filter(.,ID!="sub-CBPD0020")

# Data Cleaning -----------------------------------------------------------
#make factor variables as factors
main$male <- factor(main$male, levels=c(0,1), labels=c("Female", "Male"))
main$part_coef <- (main$part_coef_neg+main$part_coef_pos)/2

#make sure NA in nVolsCensored doesn't turn into missing data
main$nVolsCensored[is.na(main$nVolsCensored)]<- 0
main$nVolCensored[is.na(main$nVolCensored)]<- 0

#make race simpler
main$race2 <- ifelse(main$race_americanindian==1, 3, ifelse(main$race_asian == 1, 3, ifelse(main$race_hawaiian==1, 3, ifelse(main$race_black==1, 2, ifelse(main$race_white==1, 1, ifelse(main$race_americanindian+main$race_asian+main$race_black+main$race_hawaiian+main$race_white > 1, 4, NA))))))
main$race2 <- factor(main$race2, labels=c("White", "Black", "Other"))
main$ethnicity <- factor(main$ethnicity, labels=c("Not Hispanic or Latino","Hispanic or Latino"))
#filter out some unneeded variables
main <- main %>% dplyr::select(., -c(colorado_child_temperament_index_timestamp:ccti_sum))
#main <- main %>% select(.,-c(cbcl_18mo_admin:cbcl_6yr_complete))
#main <- main %>% select(.,-c(has_diagnoses:letterword_identification_comple))

main$ses_composite <- as.numeric(scale(main$parent1_edu)+scale(main$income_median))

#make categorical child aces
main$aces3category <- ifelse(main$childaces_sum_ignorenan <= 1, 0, ifelse(main$childaces_sum_ignorenan == 2, 1, ifelse(main$childaces_sum_ignorenan>=3, 2, NA)))
main$aces3category <- factor(main$aces3category, labels=c("None or One", "Two", "Three+"))

#Filter out runs from participants with < 50% of frames remaining after the 0.5mm and 1.75 DVARS thresholds
main$pctVolsCensored <- main$nVolsCensored/main$size_t
main$pctVolsCensored <- main$nVolCensored/main$size_t
#Take out people with mean motion over 1, or more than 50% of frames censored (pctSpikesFD), or masys motion > 10 mm
main_filt <- main %>% filter(., pctSpikesFD < 0.5 & relMaxRMSMotion < 10 & relMeanRMSMotion < 1)
#filter out only kids under 8
main_under9 <- filter(main_filt, age_scan <=9)
main_unique <- main_filt %>% group_by(ID) %>% filter(row_number() == 1)

#then melt it into wide format and average across participants with more than one run?
view(dfSummary(main))

# Load Replicate in Schaefer200 Data ---------------------------------------------------------------
subjdata_dir="~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/subjectData/"
outdir="~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_censor_5contig_fd0.5dvars1.75_drpvls/"
netdata_dir="~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_censor_5contig_fd0.5dvars1.75_drpvls/"

ages <- read.csv(paste0(subjdata_dir,"CBPD_data_190729_age_ses_ccti.csv")) #find most recent CBPD data wherever it is
main_replicate <- read.csv(paste0(subjdata_dir,"CBPD_parent_questionnaires_upto_177_080819.csv")) #find most recent CBPD data wherever it is
withavgweight <- read.csv(paste0(netdata_dir,"n75_within_between_Yeo7_Schaefer400_gsr_censor_5contig_fd0.5dvars1.75_drpvls_withmodulpartcoef_with_QA.csv"))

#filter out extra variables in average weight
withavgweight$record_id <- withavgweight$ID
#merge in main_replicate CBPD data
main_replicate$ID <- paste0("sub-",main_replicate$record_id)
main_replicate <- merge(withavgweight,main_replicate, by="ID")
#keep only ages, merge those in
ages$ID <- paste0("sub-",ages$record_id)
ages <- ages %>% dplyr::select(., ID,dob_entered:age_ques)
main_replicate <- merge(main_replicate,ages, by="ID")


#take out the participant with the glitter in her hair and artifact in rest?
main_replicate <- main_replicate %>% filter(.,ID!="sub-CBPD0020")

# Data Cleaning -----------------------------------------------------------
#make factor variables as factors
main_replicate$male <- factor(main_replicate$male, levels=c(0,1), labels=c("Female", "Male"))
main_replicate$part_coef <- (main_replicate$part_coef_neg+main_replicate$part_coef_pos)/2

#make sure NA in nVolsCensored doesn't turn into missing data
main_replicate$nVolsCensored[is.na(main_replicate$nVolsCensored)]<- 0
main_replicate$nVolCensored[is.na(main_replicate$nVolCensored)]<- 0

#make race simpler
main_replicate$race2 <- ifelse(main_replicate$race_americanindian==1, 3, ifelse(main_replicate$race_asian == 1, 3, ifelse(main_replicate$race_hawaiian==1, 3, ifelse(main_replicate$race_black==1, 2, ifelse(main_replicate$race_white==1, 1, ifelse(main_replicate$race_americanindian+main_replicate$race_asian+main_replicate$race_black+main_replicate$race_hawaiian+main_replicate$race_white > 1, 4, NA))))))
main_replicate$race2 <- factor(main_replicate$race2, labels=c("White", "Black", "Other"))
main_replicate$ethnicity <- factor(main_replicate$ethnicity, labels=c("Not Hispanic or Latino","Hispanic or Latino"))
#filter out some unneeded variables
main_replicate <- main_replicate %>% dplyr::select(., -c(colorado_child_temperament_index_timestamp:ccti_sum))
#main_replicate <- main_replicate %>% select(.,-c(cbcl_18mo_admin:cbcl_6yr_complete))
#main_replicate <- main_replicate %>% select(.,-c(has_diagnoses:letterword_identification_comple))

main_replicate$ses_composite <- as.numeric(scale(main_replicate$parent1_edu)+scale(main_replicate$income_median))

#make categorical child aces
main_replicate$aces3category <- ifelse(main_replicate$childaces_sum_ignorenan <= 1, 0, ifelse(main_replicate$childaces_sum_ignorenan == 2, 1, ifelse(main_replicate$childaces_sum_ignorenan>=3, 2, NA)))
main_replicate$aces3category <- factor(main_replicate$aces3category, labels=c("None or One", "Two", "Three+"))

#Filter out runs from participants with < 50% of frames remain_replicateing after the 0.5mm and 1.75 DVARS thresholds
main_replicate$pctVolsCensored <- main_replicate$nVolsCensored/main_replicate$size_t
main_replicate$pctVolsCensored <- main_replicate$nVolCensored/main_replicate$size_t
#Take out people with mean motion over 1, or more than 50% of frames censored (pctSpikesFD), or masys motion > 10 mm
main_replicate_filt <- main_replicate %>% filter(., pctSpikesFD < 0.5 & relMaxRMSMotion < 10 & relMeanRMSMotion < 1)
#filter out only kids under 8
main_replicate_under9 <- filter(main_replicate_filt, age_scan <=9)
main_replicate_unique <- main_replicate_filt %>% group_by(ID) %>% filter(row_number() == 1)

# Look at each column and correct for multiple comparisons
covariates="~ age_scan+male+fd_mean+avgweight+pctSpikesFD+size_t"

#make a dataframe with no repeats of net comparisons
main_unique <- dplyr::select(main_unique, -c(sys2to1,sys3to1,sys3to2,sys4to1,sys4to2,sys4to3,sys5to1,sys5to2,sys5to3,sys5to4,sys6to1,sys6to2,sys6to3,sys6to4,sys6to5,sys7to1,sys7to2,sys7to3,sys7to4,sys7to5,sys7to6))
#run and compare for multiple comparisons again
m <- mclapply(names(main_unique[,50:77]), function(sys) {as.formula(paste(sys, covariates, sep=""))},mc.cores=2)
networks_Age_pvals <- mclapply(m, function(sys) { summary(lm(formula = sys,data=main_unique))$coefficients[2,4]},mc.cores=1)
networks_Age_pvals <- as.data.frame(networks_Age_pvals)
networks_Age_pvals <- t(networks_Age_pvals)
networks_Age_pvals <- as.data.frame(networks_Age_pvals)
#bonferroni correct
networks_Age_pvals_fdr <- p.adjust(networks_Age_pvals$V1, method="fdr")
networks_age_pvals_fdr <- data.frame(networks_Age_pvals_fdr,names(main_unique[,50:77]))
colnames(networks_age_pvals_fdr) <- c("pvalue", "network")
#FDR correction shows DMN to attentional networks and visual to dorsal attention is marginal.

#make a dataframe with no repeats of net comparisons
main_replicate_unique <- dplyr::select(main_replicate_unique, -c(sys2to1,sys3to1,sys3to2,sys4to1,sys4to2,sys4to3,sys5to1,sys5to2,sys5to3,sys5to4,sys6to1,sys6to2,sys6to3,sys6to4,sys6to5,sys7to1,sys7to2,sys7to3,sys7to4,sys7to5,sys7to6))
#run and compare for multiple comparisons again
m <- mclapply(names(main_replicate_unique[,50:77]), function(sys) {as.formula(paste(sys, covariates, sep=""))},mc.cores=2)
networks_Age_pvals <- mclapply(m, function(sys) { summary(lm(formula = sys,data=main_replicate_unique))$coefficients[2,4]},mc.cores=1)
networks_Age_pvals <- as.data.frame(networks_Age_pvals)
networks_Age_pvals <- t(networks_Age_pvals)
networks_Age_pvals <- as.data.frame(networks_Age_pvals)
#bonferroni correct
networks_Age_pvals_fdr_replicate <- p.adjust(networks_Age_pvals$V1, method="fdr")
networks_age_pvals_fdr_replicate <- data.frame(networks_Age_pvals_fdr_replicate,names(main_replicate_unique[,78:105]))
#FDR correction shows DMN to attentional networks and visual to dorsal attention is marginal.
colnames(networks_age_pvals_fdr_replicate) <- c("pvalue", "network")

# Separate Networks and Effects of SES ------------------------------------
covariates="~ age_scan+male+fd_mean+avgweight+pctSpikesFD+size_t+aces3category"
#make a dataframe with no repeats of net comparisons
main_unique <- dplyr::select(main_unique, -c(sys2to1,sys3to1,sys3to2,sys4to1,sys4to2,sys4to3,sys5to1,sys5to2,sys5to3,sys5to4,sys6to1,sys6to2,sys6to3,sys6to4,sys6to5,sys7to1,sys7to2,sys7to3,sys7to4,sys7to5,sys7to6))
#run and compare for multiple comparisons again
m <- mclapply(names(main_unique[,54:81]), function(sys) {as.formula(paste(sys, covariates, sep=""))},mc.cores=2)
networks_ses_pvals <- mclapply(m, function(sys) { summary(lm(formula = sys,data=main_unique))$coefficients[8,4]},mc.cores=1)
networks_ses_pvals <- as.data.frame(networks_ses_pvals)
networks_ses_pvals <- t(networks_ses_pvals)
networks_ses_pvals <- as.data.frame(networks_ses_pvals)
#bonferroni correct
networks_ses_pvals_fdr <- p.adjust(networks_ses_pvals$V1, method="fdr")
networks_ses_pvals_fdr <- data.frame(networks_ses_pvals_fdr,names(main_unique[,50:77]))
colnames(networks_ses_pvals_fdr) <- c("pvalue", "network")
networks_ses_pvals_fdr

#make a dataframe with no repeats of net comparisons
main_replicate_unique <- dplyr::select(main_replicate_unique, -c(sys2to1,sys3to1,sys3to2,sys4to1,sys4to2,sys4to3,sys5to1,sys5to2,sys5to3,sys5to4,sys6to1,sys6to2,sys6to3,sys6to4,sys6to5,sys7to1,sys7to2,sys7to3,sys7to4,sys7to5,sys7to6))
#run and compare for multiple comparisons again
m <- mclapply(names(main_replicate_unique[,54:81]), function(sys) {as.formula(paste(sys, covariates, sep=""))},mc.cores=2)
networks_ses_pvals <- mclapply(m, function(sys) { summary(lm(formula = sys,data=main_replicate_unique))$coefficients[8,4]},mc.cores=1)
networks_ses_pvals <- as.data.frame(networks_ses_pvals)
networks_ses_pvals <- t(networks_ses_pvals)
networks_ses_pvals <- as.data.frame(networks_ses_pvals)
#bonferroni correct
networks_ses_pvals_fdr_replicate <- p.adjust(networks_ses_pvals$V1, method="fdr")
networks_ses_pvals_fdr_replicate <- data.frame(networks_ses_pvals_fdr_replicate,names(main_unique[,50:77]))
colnames(networks_ses_pvals_fdr_replicate) <- c("pvalue", "network")
networks_ses_pvals_fdr_replicate
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
cogstim_variables <- cogstim_variables%>% dplyr::select(.,-c(pss_q4_notreversed,pss_q5_notreversed,pss_q8_notreversed, pss_q7_notreversed))
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

# Save models for use in markdown file ------------------------------------
save(main, main_filt, main_unique, networks_age_pvals_fdr,networks_ses_pvals_fdr, networks_age_pvals_fdr_replicate,networks_ses_pvals_fdr_replicate, main_replicate, main_replicate_filt, main_replicate_unique, file=paste0(outdir,"CBPD_n75_schaefer400.Rdata"))

# Age and measures of segregation -----------------------------------------
#Need to control for the amount of good data a participant has, size_t

#look at mean within and between with age
lm_within_sys_age <- lm(mean_within_sys~age_scan+male+fd_mean+avgweight+pctSpikesFD+size_t, data=main_unique)
summary(lm_within_sys_age)
lm.beta(l)
lm_between_sys_age <- lm(mean_between_sys~age_scan+male+fd_mean+avgweight+pctSpikesFD+size_t, data=main_unique)
summary(lm_between_sys_age)
lm.beta(l)
lm_segreg_age<- lm(system_segreg~age_scan+male+fd_mean+avgweight+pctSpikesFD+size_t, data=main_unique)
summary(lm_segreg_age)
lm_modul_age <- lm(modul~age_scan+male+fd_mean+avgweight+pctSpikesFD+size_t, data=main_unique)
summary(lm_modul_age)
lm_part_coef_age <- lm(part_coef~age_scan+male+fd_mean+avgweight+pctSpikesFD+size_t, data=main_unique)
summary(lm_part_coef_age)

#Does number of communities detected with modul change with age?
lm_num_comms_age <- lm(num_comms_modul_avg~age_scan+male+fd_mean+avgweight+pctSpikesFD+size_t, data=main_unique)
summary(l) #Decreases slightly with age
#look at average weight with age
l <- lm(avgweight~age_scan+male+fd_mean+pctSpikesFD+size_t, data=main_unique)
summary(l)
lm.beta(l)

# Non-linear effects of age? ----------------------------------------------
gam_part_coef_age <- gam(part_coef~s(age_scan)+male+fd_mean+avgweight+pctSpikesFD+size_t, data=main_unique)
l <- gam(modul~s(ageAtScan1cent)+sex+avgweight+envSES, data=master, method = "REML")

#This is in the RMarkdown document.

# Separate Networks and Effects of Age ------------------------------------
networks <- select(main_unique, sys1to1:sys7to7)
nets <- colnames(networks)
#look at non-linear interaction between age and envSES 
for (net in nets){
  name<-paste0("lm_",net)
  formula<-formula(paste0(net, '~age_scan+male+fd_mean+avgweight+pctSpikesFD+size_t'))
  assign(name, lm(formula, data=main_unique))
  #p_val[net] <- summary(name)$coefficients[2,4]
}
summary(lm_sys1to1) #increase with age
summary(lm_sys1to2)
summary(lm_sys1to3) #increase with age
lm.beta(lm_sys1to3)
summary(lm_sys1to4)
summary(lm_sys1to5)
summary(lm_sys1to6)
summary(lm_sys1to7)
summary(lm_sys2to2)
summary(lm_sys2to3)
summary(lm_sys2to4)
summary(lm_sys2to5)
summary(lm_sys2to6) #effect of SES composite
summary(lm_sys2to7)
summary(lm_sys6to6)
summary(lm_sys6to3)
summary(lm_sys6to4)
summary(lm_sys6to5)
summary(lm_sys6to6)
summary(lm_sys6to7)
summary(lm_sys7to7)#increase with age
summary(lm_sys7to3) #marginal decrease with age
lm.beta(lm_sys7to3)
summary(lm_sys7to4)
visreg(lm_sys7to4) #strong decrease with age
summary(lm_sys7to5)

# Environmental effects on networks -------------------------------------------------

#look at segreg measures with SES
lm_within_sys_age <- lm(mean_within_sys~age_scan+male+fd_mean+avgweight+pctSpikesFD+size_t+neighborhood_perceived_quality, data=main_unique)
summary(lm_within_sys_age)
lm.beta(l)
lm_between_sys_age <- lm(mean_between_sys~age_scan+male+fd_mean+avgweight+pctSpikesFD+size_t+neighborhood_perceived_quality, data=main_unique)
summary(lm_between_sys_age)
lm.beta(l)
lm_segreg_age<- lm(system_segreg~age_scan+male+fd_mean+avgweight+pctSpikesFD+size_t+neighborhood_perceived_quality, data=main_unique)
summary(lm_segreg_age)
lm_modul_age <- lm(modul_avg~age_scan+male+fd_mean+avgweight+pctSpikesFD+size_t+neighborhood_perceived_quality, data=main_unique)
summary(lm_modul_age)
lm_part_coef_age <- lm(part_coef~age_scan+male+fd_mean+avgweight+pctSpikesFD+size_t+neighborhood_perceived_quality, data=main_unique)
summary(lm_part_coef_age)

visreg(lm_between_sys_age)
visreg(lm_within_sys_age)
visreg(lm_segreg_age)
visreg(lm_modul_age)
visreg(lm_part_coef_age)
visreg(lm_num_comms_age)

#look at whether there is an interaction
lm_within_sys_age_income <- lm(mean_within_sys~age_scan*neighborhood_perceived_quality+male+fd_mean+avgweight+pctSpikesFD+size_t, data=main_unique)
summary(lm_within_sys_age_income)
lm.beta(l)
lm_between_sys_age_income <- lm(mean_between_sys~age_scan*neighborhood_perceived_quality+male+fd_mean+avgweight+pctSpikesFD+size_t, data=main_unique)
summary(lm_between_sys_age_income)
lm.beta(l)
lm_segreg_age<- lm(system_segreg~age_scan*neighborhood_perceived_quality+male+fd_mean+avgweight+pctSpikesFD+size_t, data=main_unique)
summary(lm_segreg_age)
lm_modul_age <- lm(modul_avg~age_scan*neighborhood_perceived_quality+male+fd_mean+avgweight+pctSpikesFD+size_t, data=main_unique)
summary(lm_modul_age)
lm_part_coef_age <- lm(part_coef~age_scan*neighborhood_perceived_quality+male+fd_mean+avgweight+pctSpikesFD+size_t, data=main_unique)
summary(lm_part_coef_age)

visreg(lm_within_sys_age, "age_scan", by="income_median")
visreg(lm_between_sys_age, "age_scan", by="income_median")

# Specific Networks and Environment ---------------------------------------
nets=c("sys1to3", "sys3to7", "sys4to7")
for (net in nets){
  name<-paste0("lm_",net)
  formula<-formula(paste0(net, '~age_scan+male+fd_mean+avgweight+pctSpikesFD+size_t'))
  assign(name, lm(formula, data=main_unique))
  #p_val[net] <- summary(name)$coefficients[2,4]
}
summary(lm_sys1to3)
summary(lm_sys3to7)
summary(lm_sys4to7)
visreg(lm_sys1to3)
visreg(lm_sys3to7)

for (net in nets){
  name<-paste0("lm_",net)
  formula<-formula(paste0(net, '~age_scan*aces3category+male+fd_mean+avgweight+pctSpikesFD+size_t+ses_composite'))
  assign(name, lm(formula, data=main_unique))
  #p_val[net] <- summary(name)$coefficients[2,4]
}
summary(lm_sys1to3)
summary(lm_sys3to7)
summary(lm_sys4to7)


