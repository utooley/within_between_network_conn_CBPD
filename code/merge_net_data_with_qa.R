library(dplyr)
library(psych)
library(mgcv)
library(stringi)
library(stringr)
library(R.matlab)

# SETUP -------------------------------------------------------------------
#Cluster mounted locally on personal computer
netdatadir="~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_censor_5contig_fd1.25dvars2_drpvls/"
sublistdir="~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/subjectLists/"
qadir="~/Desktop/cluster/picsl/mackey_group/BPD/CBPD_bids/derivatives/mriqc_fd_2_mm/"
xcpdir="~/Desktop/cluster/picsl/mackey_group/BPD/CBPD_bids/derivatives/xcpEngine_gsr_censor_5contig_fd1.25dvars2_drpvls/"
analysis_dir="~/Documents/bassett_lab/tooleyEnviNetworks/analyses/"


# Read in files -----------------------------------------------------------
#net data
#file1<-read.csv(paste0(netdatadir, "n47_within_between_Yeo7_Schaefer400.csv"))
file1 <- read.csv(paste0(netdatadir, "n64_within_between_Yeo7_Schaefer400_withmodul.csv"))
#MRIQC Data
file2<-read.table(paste0(qadir, "group_bold.tsv"), sep = '\t', header = TRUE)
#subject list
subjlist <- read.csv(paste0(sublistdir, "n64_cohort_mult_runs_usable_t1_rest_1mm_outliers_10_2mm_060419.csv"), header = TRUE)
#xcp quality data
qa2 <- read.csv(paste0(xcpdir, "XCP_QAVARS.csv"))

# Data Cleaning -----------------------------------------------------------
file1<-dplyr::rename(file1, ID=Var1)
file1<-dplyr::rename(file1, run=Var2)
subjlist<-dplyr::rename(subjlist, run=id1)
subjlist<-dplyr::rename(subjlist, ID=id0)
qa2<-dplyr::rename(qa2, run=id1)
qa2<-dplyr::rename(qa2, ID=id0)
subjlist <- dplyr::select(subjlist, -img)

#make ID a character vector
subjlist$ID <- as.character(subjlist$ID)
file1$ID <- as.character(file1$ID)

#split QA file ID name on underscores and extract run number
file2$ID <-  stri_split_fixed(file2$bids_name,"_", simplify = TRUE)[,1]
file2$scan_type <-stri_split_fixed(file2$bids_name,"_", simplify = TRUE)[,2] 
file2$run <- stri_split_fixed(file2$bids_name,"_", simplify = TRUE)[,3]
#file2$run <- as.numeric(str_replace(file2$run,"run-0",""))
file2 <- moveMe(file2, c("ID", "run","scan_type"), "after", "bids_name")

#write the QA file back out into the MRIQC folder for future use
write.table(file2, paste0(qadir, "group_bold_with_names.tsv"), sep = '\t')

# Match up Datafiles ------------------------------------------------------
#Filter QA file for rest only
file2 <- file2 %>% filter(.,scan_type=="task-rest")

#add 'sub' prefix to the subject list so it matches
subjlist$ID <- paste0("sub-",subjlist$ID)
#merge the network data with the subject list with the run that was used to calculate it
master<-right_join(subjlist,file1, by=c("ID", "run"))
#merge in the QA data, ignoring runs that were not used for network calculations
master <- right_join(file2,master, by=c("ID", "run"))
#merge in the xcp quality data
master <- right_join(qa2,master, by=c("ID", "run"))

#filter out extraneous QA variables 
master <- master %>% select(., -c(aor:fber)) %>% select(.,-c(spacing_tr:summary_fg_stdv))

# Make a second rest run a second column? ----------------------------------------------------------
## Include number of volumes and the number of bad vols/outliers/censored vols in each run 

## Include baseline number of vols in each run from MRIQC

## Include Jaccard and other indices, make sure no outliers?

# Write out Data ----------------------------------------------------------

#write the network data file back into the output folder
write.csv(master,"~/Downloads/n64_within_between_Yeo_Schaefer400_gsr_censor_5contig_fd1.25dvars2_drpvls_with_qa.csv")
write.csv(master,paste0(netdatadir,"n64_within_between_Yeo_Schaefer400_gsr_censor_5contig_fd1.25dvars2_drpvls_with_QA.csv"))

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
