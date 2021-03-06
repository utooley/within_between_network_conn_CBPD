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
  
# Loop through each pipeline ----------------------------------------------
#pipelines=c('nogsr_spkreg_fd0.5dvars1.75_drpvls','gsr_spkreg_fd0.5dvars1.75_drpvls', "gsr_censor_5contig_fd0.5dvars1.75_drpvls", "gsr_censor_5contig_fd1.25dvars2_drpvls", "nogsr_spkreg_fd1.25dvars2_drpvls")
pipelines=c("gsr_censor_5contig_fd0.5dvars1.75_drpvls", "gsr_censor_5contig_fd1.25dvars2_drpvls", "nogsr_spkreg_fd1.25dvars2_drpvls")
for (pipeline in pipelines){

# SETUP -------------------------------------------------------------------
#Cluster mounted locally on personal computer
netdatadir=paste0("~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/imageData/", pipeline)
localnetdatadir=paste0("/Users/utooley/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/", pipeline)
sublistdir="~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/subjectLists/"
qadir="~/Desktop/cluster/picsl/mackey_group/BPD/CBPD_bids/derivatives/mriqc_fd_2_mm/"
xcpdir=paste0("~/Desktop/cluster/picsl/mackey_group/BPD/CBPD_bids/derivatives/xcpEngine_", pipeline)
analysis_dir="~/Documents/bassett_lab/tooleyEnviNetworks/analyses/"

# Read in files -----------------------------------------------------------
#net data
#file1<-read.csv(paste0(netdatadir, "n47_within_between_Yeo7_Schaefer400.csv"))
if (run=="averaged"){
file1 <- read.csv(paste0(netdatadir,"/n74_within_between_Yeo7_avgruns_",parcellation,"withmodulpartcoef.csv"))
} else {
file1 <- read.csv(paste0(netdatadir,"/n74_within_between_Yeo7_",parcellation,"withmodulpartcoef.csv"))
}
#MRIQC Data
file2<-read.table(paste0(qadir, "group_bold.tsv"), sep = '\t', header = TRUE)
#subject list
subjlist <- read.csv(paste0(sublistdir, "n74_cohort_mult_runs_usable_t1_rest_1mm_outliers_10_2mm_80119.csv"), header = TRUE)
#xcp quality data, if it's the older fixed one or in a newer folder.
if (file.exists(paste0(xcpdir, "/XCP_QAVARS.csv"))) {
  qa2 <- read.csv(paste0(xcpdir, "/XCP_QAVARS.csv")) 
} else {
  qa2 <- read.csv(paste0(xcpdir, "/XCP_QAVARS.csv")) 
}

# Data Cleaning -----------------------------------------------------------
#file1<-dplyr::rename(file1, ID=subjlist)
subjlist<-dplyr::rename(subjlist, run=id1)
subjlist<-dplyr::rename(subjlist, ID=id0)
qa2<-dplyr::rename(qa2, run=id1)
qa2<-dplyr::rename(qa2, ID=id0)
subjlist <- dplyr::select(subjlist, -img)
file2<-dplyr::rename(file2, fd_num_2mm=fd_num)

#make ID a character vector
subjlist$ID <- as.character(subjlist$ID)
file1$ID <- as.character(file1$ID)
file1$ID <- trimws(file1$ID) #take off any extra whitespace that might impede

#split QA file ID name on underscores and extract run number
file2$ID <-  stri_split_fixed(file2$bids_name,"_", simplify = TRUE)[,1]
file2$scan_type <-stri_split_fixed(file2$bids_name,"_", simplify = TRUE)[,2] 
file2$run <- stri_split_fixed(file2$bids_name,"_", simplify = TRUE)[,3]
#file2$run <- as.numeric(str_replace(file2$run,"run-0",""))
file2 <- moveMe(file2, c("ID", "run","scan_type"), "after", "bids_name")

#write the QA file back out into the MRIQC folder for future use
write.csv(file2, paste0(qadir, "group_bold_with_names.csv"))

# Match up Datafiles ------------------------------------------------------
#Filter QA file for rest only
file2 <- file2 %>% filter(.,scan_type=="task-rest")

#add 'sub' prefix to the subject list so it matches
#subjlist$ID <- paste0("sub-",subjlist$ID)
#merge the network data with the subject list with the run that was used to calculate it
#master<-right_join(subjlist,file1, by=c("ID", "run"))
if (run=="averaged"){
  master<-right_join(subjlist,file1, by=c("ID")) #now that we've averaged both runs together, just merge on ID
} else if (run == "both"){
master<-right_join(subjlist,file1, by=c("ID", "run")) 
}
#merge in the QA data, ignoring runs that were not used for network calculations
master <- right_join(file2,master, by=c("ID", "run"))
#merge in the xcp quality data
master <- right_join(qa2,master, by=c("ID", "run"))


# Summarise any run-wise statistics ---------------------------------------
master$nVolCensored[is.na(master$nVolCensored)]<- 0
#filter out extraneous QA variables and make summary variables of volumes and censored volumes
master <- master %>% dplyr::select(., -c(aor:fber)) %>% dplyr::select(.,-c(spacing_tr:summary_fg_stdv)) %>% group_by(ID) %>% 
  mutate(totalnVolCensored=sum(nVolCensored), totalSizet=sum(size_t)) %>% ungroup()

master <- master %>% mutate(perc_vols=size_t/totalSizet, fd_mean_weight=fd_mean*perc_vols, pctSpikesFD_weight=pctSpikesFD*perc_vols) %>% group_by(ID) %>% 
  mutate(fd_mean_avg=sum(fd_mean_weight), pctVolsCensored=(totalnVolCensored/totalSizet), pctSpikesFD_avg=sum(pctSpikesFD_weight))
#average motion across the two runs, weighted by the length of each run as a percentage of the total, same for percent spikes FD.
# Make a second rest run a second column? ----------------------------------------------------------
## Include number of volumes and the number of bad vols/outliers/censored vols in each run 

## Include baseline number of vols in each run from MRIQC

## Include Jaccard and other indices, make sure no outliers?

# Write out Data ----------------------------------------------------------
#create directory if it doesn't already exist
dir.create(localnetdatadir)
#write the network data file back into the output folder
if (run=="averaged"){
  write.csv(master,paste0("~/Downloads/n74_within_between_Yeo7_avgruns_",parcellation, pipeline,"_withmodulpartcoef_with_QA.csv"))
  write.csv(master,paste0(netdatadir,"/n74_within_between_Yeo7_avgruns_",parcellation,pipeline,"_withmodulpartcoef_with_QA.csv"))
  write.csv(master,paste0(localnetdatadir,"/n74_within_between_Yeo7_avgruns_",parcellation,pipeline,"_withmodulpartcoef_with_QA.csv"))
} else {
  write.csv(master,paste0("~/Downloads/n74_within_between_Yeo7_",parcellation, pipeline,"_withmodulpartcoef_with_QA.csv"))
  write.csv(master,paste0(netdatadir,"/n74_within_between_Yeo7_",parcellation,pipeline,"_withmodulpartcoef_with_QA.csv"))
  write.csv(master,paste0(localnetdatadir,"/n74_within_between_Yeo7_",parcellation,pipeline,"_withmodulpartcoef_with_QA.csv"))
}
}
}


# Just QA data for MATLAB -------------------------------------------------
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
    #filter out extraneous QA variables and make summary variables of volumes and censored volumes
    master <- master %>% dplyr::select(., -c(aor:fber)) %>% dplyr::select(.,-c(spacing_tr:summary_fg_stdv)) %>% group_by(ID) %>% 
      mutate(totalnVolCensored=sum(nVolCensored), totalSizet=sum(size_t)) %>% ungroup()
    #write it out on the cluster
    #create directory if it doesn't already exist
    dir.create(localnetdatadir)
    dir.create(netdatadir)
    write.csv(master,paste0(netdatadir,"/n150_within_between_Yeo7_",parcellation,pipeline,"_QA.csv"))
    write.csv(master,paste0(localnetdatadir,"/n150_within_between_Yeo7_",parcellation,pipeline,"_QA.csv"))
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
