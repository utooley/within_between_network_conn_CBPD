library(dplyr)
library(psych)
library(mgcv)
library(stringi)
library(stringr)
library(R.matlab)

# Just QA data for MATLAB, do this first -------------------------------------------------
parcellations=c("schaefer400_","schaefer200_")
for (parcellation in parcellations){
  pipeline="gsr_spkreg_fd0.25dvars1.75_drpvls"
  #pipelines=c('nogsr_spkreg_fd0.5dvars1.75_drpvls','gsr_spkreg_fd0.5dvars1.75_drpvls', "gsr_censor_5contig_fd0.5dvars1.75_drpvls", "gsr_censor_5contig_fd1.25dvars2_drpvls", "nogsr_spkreg_fd1.25dvars2_drpvls")
  #pipelines=c("nogsr_spkreg_fd1.25dvars2_drpvls",'gsr_spkreg_fd0.5dvars1.75_drpvls','nogsr_spkreg_fd0.5dvars1.75_drpvls')
  #for (pipeline in pipelines){
  netdatadir=paste0("~/cbica/home/tooleyu//projects/in_progress/within_between_network_conn_CBPD/data/imageData/", pipeline)
  localnetdatadir=paste0("/Users/utooley/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/", pipeline)
  sublistdir="~/cbica/home/tooleyu//projects/in_progress/within_between_network_conn_CBPD/data/subjectLists/"
  qadir="~/cbica/projects/cbpd_main_data/CBPD_bids_crosssectional//derivatives/mriqc_fd_2_mm/"
  xcpdir=paste0("~/cbica/projects/cbpd_main_data/CBPD_bids_crosssectional//derivatives/xcpEngine_", pipeline)
  #MRIQC Data
  file2<-read.table(paste0(qadir, "group_bold.tsv"), sep = '\t', header = TRUE)
  #subject list
  subjlist <- read.csv(paste0(sublistdir, "n125_fd0.25_pipeline_reprocess_revisions_for_xcpEngine.csv"), header = TRUE)
  #xcp quality data, if it's the older fixed one or in a newer folder.
  if (file.exists(paste0(xcpdir, "/XCP_QAVARS_xcpEngine_", pipeline, "_fixed.csv"))) {
    qa2 <- read.csv(paste0(xcpdir, "/XCP_QAVARS_xcpEngine_",pipeline, "_fixed.csv")) #this is the group xcp output csv that I looked at and fixed before this
  } else {
    qa2 <- read.csv(paste0(xcpdir, "/XCP_QAVARS_xcpEngine_", pipeline, ".csv")) 
  }
  #subjlist<-dplyr::rename(subjlist, run=id1)
  #subjlist<-dplyr::rename(subjlist, ID=id0)
  qa2<-dplyr::rename(qa2, run=id1)
  qa2<-dplyr::rename(qa2, ID=id0)
  #subjlist <- dplyr::select(subjlist, -img)
  file2<-dplyr::rename(file2, fd_num_2mm=fd_num)
  #make ID a character vector
  subjlist$ID <- as.character(subjlist$id0)
  subjlist$run <- as.character(subjlist$id1)
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
  #master <- filter(master, relMaxRMSMotion< 10)
  #filter out extraneous QA variables and make summary variables of volumes and censored volumes
  master <- master %>% dplyr::select(., -c(aor:fber)) %>% dplyr::select(.,-c(spacing_tr:summary_fg_stdv)) %>% group_by(ID) %>% 
    mutate(totalnVolCensored=sum(nVolCensored), totalSizet=sum(size_t)) %>% ungroup()
  #exclude people who have less than 130 vols total after this
  master <- filter(master, totalSizet >= 130)
  #write it out on the cluster
  #create directory if it doesn't already exist
  dir.create(localnetdatadir)
  dir.create(netdatadir)
  write.csv(master,paste0(netdatadir,"/n125_within_between_Yeo7_",parcellation,pipeline,"_QA_revisions.csv"))
  write.csv(master,paste0(localnetdatadir,"/n125_within_between_Yeo7_",parcellation,pipeline,"_QA_revisions.csv"))
}
}

# MoveMe Function -- run first ---------------------------------------------------------
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