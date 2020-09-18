library(dplyr)
library(stats)
library(parallel)
library(lm.beta)
library(packrat)
library(summarytools)
library(psych)
library(polycor)
library(psy)


# Load data ---------------------------------------------------------------
question_data_dir="~/Box/Mackey_Lab/CBPD (825656)/Data/PCIT_Questionasking_Ursula/"
#2 header rows
qa_data <- read.csv(paste0(question_data_dir, "question_asking_three_coders_parent_qs_and_longitudinal_082020.csv"), skip = 1)


# Check correlations ------------------------------------------------------
colnames(qa_data)
#Take only a few relevant variables
subset <- qa_data %>% dplyr::select(.,contains("_name") | contains("Total_q") | contains("Total_pq"))


# Look at Emily-Ortal (n=96) ----------------------------------------------
subset %>% filter(Coder1_name=="Emily" & Coder2_name=="Ortal" ) %>% 
  dplyr::select(contains("coder1") | contains("coder2")) %>% dplyr::select_if(is.numeric) %>% 
  cor(., use = "pairwise.complete.obs")

subset_eo <- subset %>% filter(Coder1_name=="Emily" & Coder2_name=="Ortal" ) %>% 
  dplyr::select(contains("coder1") | contains("coder2"))
icc(select(subset_eo, contains("Total_q_kid_allpcit_coder")))
ICC(select(subset_eo, contains("Total_q_kid_allpcit_coder")))

#Correlations range from 92-95.
#ICC for ICC2k (average, random raters) is 0.96 for kid q's, 0.98 for parent q's, 0.98 for parent pedagogical questions 
icc(select(subset_eo, contains("Total_q_parent_bookreading_coder")))
ICC(select(subset_eo, contains("Total_q_parent_bookreading_coder")))

icc(select(subset_eo, contains("Total_pq_parent_bookreading_coder")))
ICC(select(subset_eo, contains("Total_pq_parent_bookreading_coder")))

# Look at Christina-Ortal (n=33) ----------------------------------------------
subset %>% filter(Coder1_name=="Christina" & Coder2_name=="Ortal" ) %>% 
  dplyr::select(contains("coder1") | contains("coder2")) %>% dplyr::select_if(is.numeric) %>% 
  cor(., use = "pairwise.complete.obs")

subset_co <- subset %>% filter(Coder1_name=="Christina" & Coder2_name=="Ortal" ) %>% 
  dplyr::select(contains("coder1") | contains("coder2"))
icc(select(subset_co, contains("Total_q_kid_allpcit_coder")))
ICC(select(subset_co, contains("Total_q_kid_allpcit_coder")))

#Correlations range from 0.86-0.95-0.999.
#ICC for ICC2k (average, random raters) is 0.83 for kid q's, 0.99 for parent q's, 0.96 for parent pedagogical questions 
icc(select(subset_co, contains("Total_q_parent_bookreading_coder")))
ICC(select(subset_co, contains("Total_q_parent_bookreading_coder")))

icc(select(subset_co, contains("Total_pq_parent_bookreading_coder")))
ICC(select(subset_co, contains("Total_pq_parent_bookreading_coder")))

# Look at Christina-Adrian (n=33)------------------------------------------------
subset%>% filter(Coder1_name=="Christina" & Coder3_name=="Adrian" ) %>% 
  dplyr::select(contains("coder1") | contains("coder3")) %>% dplyr::select_if(is.numeric) %>% 
  cor(., use = "pairwise.complete.obs")

subset_ca <- subset %>% filter(Coder1_name=="Christina" & Coder3_name=="Adrian" ) %>% 
  dplyr::select(contains("coder1") | contains("coder3"))
icc(select(subset_ca, contains("Total_q_kid_allpcit_coder")))
ICC(select(subset_ca, contains("Total_q_kid_allpcit_coder")))

#Correlations range from 0.95-0.98-0.99.
#ICC for ICC2k (average, random raters) is 0.95 for kid q's, 0.99 for parent q's, 0.99 for parent pedagogical questions 
icc(select(subset_ca, contains("Total_q_parent_bookreading_coder")))
ICC(select(subset_ca, contains("Total_q_parent_bookreading_coder")))

icc(select(subset_ca, contains("Total_pq_parent_bookreading_coder")))
ICC(select(subset_ca, contains("Total_pq_parent_bookreading_coder")))

# Look at Ortal-Adrian (n=33)------------------------------------------------
subset%>% filter(Coder2_name=="Ortal" & Coder3_name=="Adrian" ) %>% 
  dplyr::select(contains("coder2") | contains("coder3")) %>% dplyr::select_if(is.numeric) %>% 
  cor(., use = "pairwise.complete.obs")

subset_oa <- subset %>% filter(Coder2_name=="Ortal" & Coder3_name=="Adrian" ) %>% 
  dplyr::select(contains("coder2")  | contains("coder3"))
icc(select(subset_oa, contains("Total_q_kid_allpcit_coder")))
ICC(select(subset_oa, contains("Total_q_kid_allpcit_coder")))

#Correlations range from 0.85-0.99-0.94.
#ICC for ICC2k (average, random raters) is 0.90 for kid q's, 0.99 for parent q's, 0.96 for parent pedagogical questions 
icc(select(subset_oa, contains("Total_q_parent_bookreading_coder")))
ICC(select(subset_oa, contains("Total_q_parent_bookreading_coder")))

icc(select(subset_oa, contains("Total_pq_parent_bookreading_coder")))
ICC(select(subset_oa, contains("Total_pq_parent_bookreading_coder")))

# Average across any and all existing raters ------------------------------
qa_data$Total_q_kid_allpcit_averaged_all <-rowMeans(cbind(qa_data$Total_q_kid_allpcit_coder1, qa_data$Total_q_kid_allpcit_coder2, qa_data$Total_q_kid_allpcit_coder3), na.rm = T)

qa_data$Total_q_parent_bookreading_averaged_all <-rowMeans(cbind(qa_data$Total_q_parent_bookreading_coder1, qa_data$Total_q_parent_bookreading_coder2, qa_data$Total_q_parent_bookreading_coder3), na.rm = T)
qa_data$Total_pq_parent_bookreading_averaged_all <-rowMeans(cbind(qa_data$Total_pq_parent_bookreading_coder1, qa_data$Total_pq_parent_bookreading_coder2, qa_data$Total_pq_parent_bookreading_coder3), na.rm = T)
qa_data$Total_ddpi_parent_bookreading_averaged_all <-rowMeans(cbind(qa_data$Total_ddpi_parent_bookreading_coder1, qa_data$Total_ddpi_parent_bookreading_coder2, qa_data$Total_ddpi_parent_bookreading_coder3), na.rm = T)
qa_data$Total_rq_parent_bookreading_averaged_all <-rowMeans(cbind(qa_data$Total_rq_parent_bookreading_coder1, qa_data$Total_rq_parent_bookreading_coder2, qa_data$Total_rq_parent_bookreading_coder3), na.rm = T)

# Make Exclude column and exclude on other factors ------------------------
summary(qa_data)
str(qa_data)
#convert to time formats
qa_data$Total_time_coded_all <- as.POSIXct(strptime(qa_data$Total_time_coded_all, format = "%H:%M:%S"))
qa_data$Total_time_coded_all_coder2 <- as.POSIXct(strptime(qa_data$Total_time_coded_all_coder2, format = "%H:%M:%S"))
qa_data$Total_time_coded_all_coder3 <- as.POSIXct(strptime(qa_data$Total_time_coded_all_coder3, format = "%H:%M:%S"))
#average total times and book-reading times
qa_data$Total_time_coded_averaged_all <-as.POSIXct(rowMeans(cbind(qa_data$Total_time_coded_all, qa_data$Total_time_coded_all_coder2, qa_data$Total_time_coded_all_coder3), na.rm = T), origin = "1970-01-01")
qa_data$Total_time_coded_book_averaged_all <-as.POSIXct(rowMeans(cbind(as.POSIXct(strptime(qa_data$Total_time_coded, format = "%H:%M:%S")), 
                                                                       as.POSIXct(strptime(qa_data$Total_time_coded_coder2, format = "%H:%M:%S")), 
                                                                       as.POSIXct(strptime(qa_data$Total_time_coded_coder3, format = "%H:%M:%S"))), na.rm = T), origin = "1970-01-01")


#exclude too short coding times
qa_data$exclude <- ifelse(qa_data$Total_time_coded_all < as.POSIXct("2020-09-18 00:16:00"), 1, 0)
qa_data$exclude2 <- ifelse(qa_data$Total_time_coded_all_coder2 < as.POSIXct("2020-09-18 00:16:00"), 1, 0)

#check inconsistencies
data.frame(qa_data$record_id, qa_data$exclude==qa_data$exclude2)

#Subjects to exclude bc didn't follow timers such that time is < 16 min or other issues with timers going off
exclude_from_kid_qasking=c("CBPD0114","CBPD0138", "CBPD0115", "CBPD0122", "CBPD0131", "CBPD0137", "CBPD0166", "CBPD0170", "CBPD0172", 
"CBPD0173", "CBPD0181", "CBPD0184", "CBPD0202", "CBPD0203", "CBPD0205", "CBPD0064_2", "CBPD0142_2", "CBPD0037_2",
"CBPD0096_2")

qa_data$exclude <- NA
qa_data$exclude <- ifelse(qa_data$record_id %in% exclude_from_kid_qasking, 1, qa_data$exclude)
maybes=c("CBPD0156", "CBPD0157", "CBPD0162")
#CBPD0195 has read the book before!

#to control for time spent in regressions need to convert to minutes
qa_data$Total_time_coded_averaged_all <- minute(qa_data$Total_time_coded_averaged_all)+(second(qa_data$Total_time_coded_averaged_all)/60)
qa_data$Total_time_coded_book_averaged_all <- minute(qa_data$Total_time_coded_book_averaged_all)+(second(qa_data$Total_time_coded_book_averaged_all)/60)
qa_data$exclude <- ifelse(qa_data$Total_time_coded_averaged_all < 16, 1, qa_data$exclude)

# Merge with demographic and brain data, exclude CBPD subjects with diagnoses  --------
demo <- read.csv("~/Box/Mackey_Lab/CBPD (825656)/Data/Data Maintenance Day/2020-06-Data/CBPD_data_DMD_2020.06.29.csv") %>% select(-c( cbcl_18mo_admin:totalDuration.adult_p_play_avg)) %>% rename(.,age_behav=age_behav.x)
#CHECK THAT THIS IS HOW WE CHOOSE TO DO RACE
demo$race2 <- ifelse(demo$race_americanindian+demo$race_asian+demo$race_black+demo$race_hawaiian+demo$race_white > 1, 4, ifelse(demo$race_americanindian==1, 3, ifelse(demo$race_other==1,3, ifelse(demo$race_asian == 1, 3, ifelse(demo$race_hawaiian==1, 3, ifelse(demo$race_black==1, 2, ifelse(demo$race_white==1, 1, NA)))))))
demo$race2 <- factor(demo$race2, labels=c("White", "Black", "Other", "Multiracial"))
demo$male <- factor(demo$male, labels=c("Female", "Male"))
#Check demo$potential_exclusion_reasons, Anne says reported diagnoses are CBPD0148, CBPD0164, CBPD0165, CBPD0210, and CBPD0211 
demo$exclude_demo <- NA
demo <- demo %>% rowwise() %>% mutate(exclude_demo=if_else(record_id %in% c("CBPD0148","CBPD0164", "CBPD0165","CBPD0210","CBPD0211"), 1, 0)) %>% ungroup()
demo <- demo %>% filter(exclude_demo!=1)
full <- left_join(qa_data,demo,by="record_id") 

# Kid Questions ---------------------------------------------------
library(PerformanceAnalytics)
library(lubridate)
#kid q's only, filter out exclusions for not enough time or full-PCIT timing issues
full_kidq <- full %>% filter(is.na(exclude))

full_kidq %>% select(age_behav,parent1_edu, income_median, epcuriosity_dtype_sum, epcuriosity_itype_sum, wisc_vocab_raw, 
                     Total_q_kid_allpcit_averaged_all, Total_q_parent_bookreading_averaged_all, Total_pq_parent_bookreading_averaged_all) %>% 
chart.Correlation()
#SES and I-D type
summary(lm(Total_q_kid_allpcit_averaged_all~age_behav+Total_time_coded_averaged_all+parent1_edu+income_median,data=full_kidq))
visreg(lm(Total_q_kid_allpcit_averaged_all~age_behav+parent1_edu+income_median,data=full_kidq))

summary(lm(Total_q_kid_allpcit_averaged_all~age_behav+Total_time_coded_averaged_all+parent1_edu+income_median+epcuriosity_itype_sum,data=full_kidq))

#Weschler-hacky combination across measures
full_kidq$full_mr_both <- ifelse(!is.na(full_kidq$wppsi_matrix_raw), full_kidq$wppsi_matrix_raw, full_kidq$wisc_matrix_raw)
full_kidq$full_vocabinfo_both <- ifelse(!is.na(full_kidq$wppsi_info_raw), full_kidq$wppsi_info_raw, full_kidq$wisc_vocab_raw)

summary(lm(Total_q_kid_allpcit_averaged_all~age_behav+Total_time_coded_averaged_all+full_vocabinfo_both,data=full_kidq))
visreg(lm(Total_q_kid_allpcit_averaged_all~age_behav+Total_time_coded_averaged_all+full_vocabinfo_both,data=full_kidq))

#relation to other measures->parent questions
summary(lm(Total_q_kid_allpcit_averaged_all~age_behav+Total_time_coded_averaged_all+Total_q_parent_bookreading_averaged_all,data=full_kidq))
visreg(lm(Total_q_kid_allpcit_averaged_all~age_behav+Total_time_coded_averaged_all+Total_q_parent_bookreading_averaged_all,data=full_kidq))

# Parent q's --------------------------------------------------------------
#should exclude folks with < 5 min?
#full_kidq1 <- filter(full_kidq, Total_time_coded_book_averaged_all > 5)

#control for time on book-reading
summary(lm(Total_q_parent_bookreading_averaged_all~age_behav+Total_time_coded_book_averaged_all,data=full_kidq))
visreg(lm(Total_q_parent_bookreading_averaged_all~age_behav+Total_time_coded_book_averaged_all,data=full_kidq))

#SES
summary(lm(Total_q_parent_bookreading_averaged_all~age_behav+parent1_edu+Total_time_coded_book_averaged_all,data=full_kidq))
summary(lm(Total_q_parent_bookreading_averaged_all~age_behav+income_median+parent1_edu+Total_time_coded_book_averaged_all,data=full_kidq))
visreg(lm(Total_q_parent_bookreading_averaged_all~age_behav+income_median+parent1_edu+Total_time_coded_book_averaged_all,data=full_kidq))

#I-D type
summary(lm(Total_q_parent_bookreading_averaged_all~age_behav+epcuriosity_dtype_sum+Total_time_coded_book_averaged_all,data=full_kidq))
summary(lm(Total_q_parent_bookreading_averaged_all~age_behav+income_median+parent1_edu+epcuriosity_itype_sum+Total_time_coded_book_averaged_all,data=full_kidq))

#Weschler
summary(lm(Total_q_parent_bookreading_averaged_all~age_behav+Total_time_coded_book_averaged_all+full_vocabinfo_both,data=full_kidq))
summary(lm(Total_q_parent_bookreading_averaged_all~age_behav+Total_time_coded_book_averaged_all+full_mr_both,data=full_kidq))
visreg(lm(Total_q_parent_bookreading_averaged_all~age_behav+Total_time_coded_book_averaged_all+full_mr_both,data=full_kidq))

#relation to other measures->parent pedagogical questions
summary(lm(Total_q_parent_bookreading_averaged_all~age_behav+Total_time_coded_book_averaged_all+Total_pq_parent_bookreading_averaged_all,data=full_kidq))
summary(lm(Total_q_parent_bookreading_averaged_all~age_behav+parent1_edu+Total_time_coded_book_averaged_all+Total_pq_parent_bookreading_averaged_all,data=full_kidq))

# Pedagogical questions ---------------------------------------------------
#control for time on book-reading
summary(lm(Total_pq_parent_bookreading_averaged_all~age_behav+Total_time_coded_book_averaged_all,data=full_kidq))
summary(lm(Total_pq_parent_bookreading_averaged_all~age_behav+Total_time_coded_book_averaged_all+Total_q_parent_bookreading_averaged_all,data=full_kidq))

#SES
summary(lm(Total_pq_parent_bookreading_averaged_all~age_behav+income_median+Total_time_coded_book_averaged_all,data=full_kidq))
summary(lm(Total_pq_parent_bookreading_averaged_all~age_behav+parent1_edu+Total_time_coded_book_averaged_all,data=full_kidq))
summary(lm(Total_pq_parent_bookreading_averaged_all~age_behav+parent1_edu+income_median+Total_time_coded_book_averaged_all,data=full_kidq))


#I-D types
summary(lm(Total_pq_parent_bookreading_averaged_all~age_behav+epcuriosity_dtype_sum+Total_time_coded_book_averaged_all,data=full_kidq))
summary(lm(Total_pq_parent_bookreading_averaged_all~age_behav+epcuriosity_dtype_sum+income_median+Total_time_coded_book_averaged_all,data=full_kidq))

#Weschler MR/Info-vocab
summary(lm(Total_pq_parent_bookreading_averaged_all~age_behav+full_mr_both+Total_time_coded_book_averaged_all,data=full_kidq))
summary(lm(Total_pq_parent_bookreading_averaged_all~age_behav+full_mr_both+parent1_edu+Total_time_coded_book_averaged_all,data=full_kidq))
summary(lm(Total_pq_parent_bookreading_averaged_all~age_behav+full_vocabinfo_both+Total_time_coded_book_averaged_all,data=full_kidq))

# Relationship with brain networks ----------------------------------------


