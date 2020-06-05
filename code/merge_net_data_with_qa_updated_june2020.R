library(dplyr)
library(psych)
library(mgcv)
library(stringi)
library(stringr)
library(R.matlab)

# Loop through each parcellation----------------------------------------------
run="both"
parcellations=c("schaefer400_","schaefer200_")
for (parcellation in parcellations){
  #parcellation="schaefer200_"
  
# Loop through each pipeline ----------------------------------------------
#pipelines=c('nogsr_spkreg_fd0.5dvars1.75_drpvls','gsr_spkreg_fd0.5dvars1.75_drpvls', "gsr_censor_5contig_fd0.5dvars1.75_drpvls", "gsr_censor_5contig_fd1.25dvars2_drpvls", "nogsr_spkreg_fd1.25dvars2_drpvls")
pipelines=c('gsr_spkreg_fd0.5dvars1.75_drpvls', 'nogsr_spkreg_fd0.5dvars1.75_drpvls', "nogsr_spkreg_fd1.25dvars2_drpvls")
for (pipeline in pipelines){
#pipeline= "nogsr_spkreg_fd1.25dvars2_drpvls"

# SETUP -------------------------------------------------------------------
#Cluster mounted locally on personal computer
netdatadir=paste0("/data/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/imageData/", pipeline)
localnetdatadir=paste0("/Users/utooley/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/", pipeline)
sublistdir="/data/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/subjectLists/"
qadir="/data/picsl/mackey_group/BPD/CBPD_bids/derivatives/mriqc_fd_2_mm/"
xcpdir=paste0("/data/picsl/mackey_group/BPD/CBPD_bids/derivatives/xcpEngine_", pipeline)
analysis_dir="~/Documents/bassett_lab/tooleyEnviNetworks/analyses/"

# Read in files -----------------------------------------------------------
#net data
#file1<-read.csv(paste0(netdatadir, "n47_within_between_Yeo7_Schaefer400.csv"))
# if (run=="averaged"){ 
file1 <- read.csv(paste0(netdatadir,"/n125_long_inc_within_between_Yeo7_avgruns_",parcellation,"withmodulpartcoef.csv"))
# } else {
# file1 <- read.csv(paste0(netdatadir,"/n74_within_between_Yeo7_",parcellation,"withmodulpartcoef.csv"))
# }
#MRIQC Data
#file2<-read.table(paste0(qadir, "group_bold.tsv"), sep = '\t', header = TRUE)
qafile<- read.csv(paste0(netdatadir, "/n125_within_between_Yeo7_", parcellation, pipeline, "_QA.csv"))
qa1mm <- read.table("/data/picsl/mackey_group/CBPD_OLD/CBPD_bids/derivatives/mriqc_fd_1_mm/group_bold.tsv", header=TRUE)
#subject list
subjlist <- read.csv(paste0(sublistdir, "n125_cross_sect_one_or_more_nonsleep_rest_10mm_max_RMS_perrun_at_least_130_vols.csv"), header = TRUE)
#xcp quality data, if it's the older fixed one or in a newer folder.
# if (file.exists(paste0(xcpdir, "/XCP_QAVARS_fixed.csv"))) {
  #qa2 <- read.csv(paste0(xcpdir, "/XCP_QAVARS_fixed.csv")) 
# } else {
#   qa2 <- read.csv(paste0(xcpdir, "/XCP_QAVARS.csv")) 
# }

# Data Cleaning -----------------------------------------------------------
#make ID a character vector
subjlist$ID <- as.character(subjlist$ID)
file1$ID <- as.character(file1$ID)
file1$ID <- trimws(file1$ID) #take off any extra whitespace that might impede

#split QA file ID name on underscores and extract run number
# file2$ID <-  stri_split_fixed(file2$bids_name,"_", simplify = TRUE)[,1]
# file2$scan_type <-stri_split_fixed(file2$bids_name,"_", simplify = TRUE)[,2] 
# file2$run <- stri_split_fixed(file2$bids_name,"_", simplify = TRUE)[,3]
# #file2$run <- as.numeric(str_replace(file2$run,"run-0",""))
# file2 <- moveMe(file2, c("ID", "run","scan_type"), "after", "bids_name")

#Split MRIQC 1 mm file on underscores, filter out non-rest, extract run number
#split QA file ID name on underscores and extract run number
qa1mm$ID <-  stri_split_fixed(qa1mm$bids_name,"_", simplify = TRUE)[,1]
#make new longitudinal IDS
qa1mm$sess <- stri_split_fixed(qa1mm$bids_name,"_", simplify = TRUE)[,2]
qa1mm$ses <- substr(qa1mm$sess, nchar(qa1mm$sess), nchar(qa1mm$sess))
qa1mm$ID <- ifelse(qa1mm$ses=="1",qa1mm$ID,paste0(qa1mm$ID,qa1mm$ses))

qa1mm$scan_type <-stri_split_fixed(qa1mm$bids_name,"_", simplify = TRUE)[,3] 
qa1mm$run <- stri_split_fixed(qa1mm$bids_name,"_", simplify = TRUE)[,4]

qa1mm <- moveMe(qa1mm, c("ID","run","scan_type"), "after", "bids_name")
qa1mm <- filter(qa1mm,scan_type=="task-rest") %>% dplyr::select(c("bids_name","ID","run","scan_type", "fd_num", "fd_perc","size_t"))
qa1mm<-dplyr::rename(qa1mm, fd_num_1mm=fd_num, fd_perc_1mm=fd_perc, size_t_1mm=size_t)

# Match up Datafiles ------------------------------------------------------
#Filter QA file for rest only
qafile <- qafile %>% filter(.,scan_type=="task-rest")

#add 'sub' prefix to the subject list so it matches
#subjlist$ID <- paste0("sub-",subjlist$ID)
#merge the network data with the subject list with the run that was used to calculate it
#master<-right_join(subjlist,file1, by=c("ID", "run"))
subjlist$ID <- trimws(subjlist$ID, which = "both")
file1$ID <- trimws(as.character(file1$ID),which = "both") #silly whitespace popping up when factor converted to character when merging
master<-right_join(subjlist,file1, by=c("ID")) #now that we've averaged both runs together, just merge on ID
#merge in the 2 mm MRIQC and XCP QA data, ignoring runs that were not used for network calculations
master <- right_join(qafile,master, by=c("ID", "run"))
#merge in the 1 mm MRIQC data
master <- right_join(qa1mm,master, by=c("ID", "run")) %>% dplyr::select(-X)

# Summarise any run-wise statistics for MRIQC 1 mm and 2 mm files---------------------------------------
master$nVolCensored[is.na(master$nVolCensored)]<- 0

#filter out extraneous QA variables and make summary variables of volumes and censored volumes-unnecessary we've already done below
# master <- master %>% dplyr::select(., -c(aor:fber)) %>% dplyr::select(.,-c(spacing_tr:summary_fg_stdv)) %>% group_by(ID) %>% 
#   mutate(totalnVolCensored=sum(nVolCensored), totalSizet=sum(size_t)) %>% ungroup()

#average motion and outliers across the two runs, weighted by the length of each run as a percentage of the total, same for percent spikes FD.
master <- master %>% mutate(perc_vols=size_t/totalSizet, fd_mean_weight=fd_mean*perc_vols, pctSpikesFD_weight=pctSpikesFD*perc_vols, 
                            fd_perc_1mm_weight=fd_perc_1mm*perc_vols, fd_perc_2mm_weight=fd_perc*perc_vols) %>% group_by(ID) %>% 
  mutate(fd_mean_avg=sum(fd_mean_weight), pctVolsCensored=(totalnVolCensored/totalSizet), pctSpikesFD_avg=sum(pctSpikesFD_weight), 
         fd_perc_1mm_avg=sum(fd_perc_1mm_weight), fd_perc_2mm_avg=sum(fd_perc_2mm_weight))

# Make a second rest run a second column? ----------------------------------------------------------
## Include number of volumes and the number of bad vols/outliers/censored vols in each run 

## Include baseline number of vols in each run from MRIQC

## Include Jaccard and other indices, make sure no outliers?

# Write out Data ----------------------------------------------------------
#create directory if it doesn't already exist
dir.create(localnetdatadir)
#write the network data file back into the output folder
write.csv(master,paste0("~/Downloads/n125_long_inc_within_between_Yeo7_avgruns_",parcellation, pipeline,"_withmodulpartcoef_with_QA.csv"))
write.csv(master,paste0(netdatadir,"/n125_long_inc_within_between_Yeo7_avgruns_",parcellation,pipeline,"_withmodulpartcoef_with_QA.csv"))
write.csv(master,paste0(localnetdatadir,"/n125_long_inc_within_between_Yeo7_avgruns_",parcellation,pipeline,"_withmodulpartcoef_with_QA.csv"))
}
}
# Just QA data for MATLAB, do this first -------------------------------------------------
parcellations=c("schaefer400_","schaefer200_")
for (parcellation in parcellations){
  pipeline="nogsr_spkreg_fd1.25dvars2_drpvls"
  #pipelines=c('nogsr_spkreg_fd0.5dvars1.75_drpvls','gsr_spkreg_fd0.5dvars1.75_drpvls', "gsr_censor_5contig_fd0.5dvars1.75_drpvls", "gsr_censor_5contig_fd1.25dvars2_drpvls", "nogsr_spkreg_fd1.25dvars2_drpvls")
  pipelines=c("nogsr_spkreg_fd1.25dvars2_drpvls",'gsr_spkreg_fd0.5dvars1.75_drpvls','nogsr_spkreg_fd0.5dvars1.75_drpvls')
  for (pipeline in pipelines){
    netdatadir=paste0("/data/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/imageData/", pipeline)
    localnetdatadir=paste0("/Users/utooley/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/", pipeline)
    sublistdir="/data/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/subjectLists/"
    qadir="/data/picsl/mackey_group/BPD/CBPD_bids/derivatives/mriqc_fd_2_mm/"
    xcpdir=paste0("/data/picsl/mackey_group/BPD/CBPD_bids/derivatives/xcpEngine_", pipeline)
    #MRIQC Data
    file2<-read.table(paste0(qadir, "group_bold.tsv"), sep = '\t', header = TRUE)
    #subject list
    subjlist <- read.csv(paste0(sublistdir, "n150_cross_sect_one_or_more_nonsleep_rest_at_least_130_vols.csv"), header = TRUE)
    #xcp quality data, if it's the older fixed one or in a newer folder.
    if (file.exists(paste0(xcpdir, "/XCP_QAVARS_fixed.csv"))) {
      qa2 <- read.csv(paste0(xcpdir, "/XCP_QAVARS_fixed.csv")) 
    } else {
      qa2 <- read.csv(paste0(xcpdir, "/XCP_QAVARS.csv")) 
    }
    #subjlist<-dplyr::rename(subjlist, run=id1)
    #subjlist<-dplyr::rename(subjlist, ID=id0)
    qa2<-dplyr::rename(qa2, run=id1)
    qa2<-dplyr::rename(qa2, ID=id0)
    #subjlist <- dplyr::select(subjlist, -img)
    file2<-dplyr::rename(file2, fd_num_2mm=fd_num)
    #make ID a character vector
    subjlist$ID <- as.character(subjlist$ID)
    #split QA file ID name on underscores and extract run number
    file2$ID <-  stri_split_fixed(file2$bids_name,"_", simplify = TRUE)[,1]
    file2$scan_type <-stri_split_fixed(file2$bids_name,"_", simplify = TRUE)[,2] 
    file2$run <- stri_split_fixed(file2$bids_name,"_", simplify = TRUE)[,3]
    #file2$run <- as.numeric(str_replace(file2$run,"run-0",""))
    file2 <- moveMe(file2, c("ID", "run","scan_type"), "after", "bids_name")
    #Filter QA file for rest only
    file2 <- file2 %>% filter(.,scan_type=="task-rest")
    #add 'sub' prefix to the subject list so it matches
    #subjlist$ID <- paste0("sub-",subjlist$ID)
    #merge the network data with the subject list with the run that was used to calculate it
    master<-right_join(file2,subjlist, by=c("ID", "run"))
    #merge in the xcp quality data
    master <- right_join(qa2,master, by=c("ID", "run"))
    #master$nVolCensored[is.na(master$nVolsCensored)]<- 0
    #filter out runs with relMaxRMS over 10 mm and then exclude people who have less than 130 vols after that
    master <- filter(master, relMaxRMSMotion< 10)
    #filter out extraneous QA variables and make summary variables of volumes and censored volumes
    master <- master %>% dplyr::select(., -c(aor:fber)) %>% dplyr::select(.,-c(spacing_tr:summary_fg_stdv)) %>% group_by(ID) %>% 
      mutate(totalnVolCensored=sum(nVolCensored), totalSizet=sum(size_t)) %>% ungroup()
    #exclude people who have less than 130 vols total after this
    master <- filter(master, totalSizet >= 130)
    #write it out on the cluster
    #create directory if it doesn't already exist
    dir.create(localnetdatadir)
    dir.create(netdatadir)
    write.csv(master,paste0(netdatadir,"/n125_within_between_Yeo7_",parcellation,pipeline,"_QA.csv"))
    write.csv(master,paste0(localnetdatadir,"/n125_within_between_Yeo7_",parcellation,pipeline,"_QA.csv"))
  }
}

# MoveMe Function ---------------------------------------------------------
moveMe <- function(data, tomove, where = "last", ba = NULL) {
  temp <- setdiff(names(data), tomove)
  x <- switch(
    where,
    first = data[c(tomove, temp)],
    last = data[c(temp, tomove)],
    before = {
      if (is.null(ba)) stop("must specify ba column")
      if (length(ba) > 1) stop("ba must be a single character string")
      data[append(temp, values = tomove, after = (match(ba, temp)-1))]
    },
    after = {
      if (is.null(ba)) stop("must specify ba column")
      if (length(ba) > 1) stop("ba must be a single character string")
      data[append(temp, values = tomove, after = (match(ba, temp)))]
    })
  x
}
