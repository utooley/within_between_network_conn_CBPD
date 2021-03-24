library(fsbrain) #this may not work if you're editing the script directly on the cluster...
library(freesurferformats)
library(dplyr)
library(R.matlab)
library(stringr)
library(ggplot2)
library(tidyr)
library(visreg)

#rearrange the order of the brains in the T9 view of fsbrain
source("~/Documents/tools/fsbrain_fix_t9.R")
environment(brainview.t9) <- asNamespace('fsbrain')
assignInNamespace("brainview.t9", brainview.t9, ns = "fsbrain")

# SETUP -------------------------------------------------------------------
#CUBIC cluster filepaths, mounted locally
subjects_dir = "/Users/utooley/Documents/tools/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/";

#make yeo colors
yeo_colors=colorRampPalette(c("#000000", "#7B287E", "#5CA1C8", "#0A9045","#C33AF8", "#dcf8a4", "#EF9C23", "#E34A53")) #these are corr
subject_id = 'fsaverage';       # for function which use one subject only
atlas='Schaefer2018_400Parcels_7Networks_order'

#Set RGL defaults, larger window and how to make snapshots!
rgloptions=list("windowRect"=c(50,50,1000,1000));

#FIGURE OUT HOW TO ROTATE ON A DIFFERENT AXIS FOR RGL
#This edits the hidden function vis.coloredmeshes.rotating, figure out a way to edit this permanently. 
#Change x=0, y=0, z=1, maybe rpm
trace(fsbrain:::brainview.sr, edit=TRUE)
trace(fsbrain::brainviews, edit=TRUE)


# Brains ------------------------------------------------------------------
subjects_dir = "/Users/utooley/Documents/tools/parcellations/Yeo_from_freesurfer/";

rglactions=list("snapshot_png"="~/Downloads/yeo_sys3_yellow.png")
vis.subject.annot(subjects_dir, 'fsaverage','Yeo2011_7Networks_N1000','both', "inflated", views="t4",
                  rgloptions = rgloptions,rglactions = rglactions)

lh <- subject.annot(subjects_dir, 'fsaverage6', 'lh','Yeo2011_7Networks_N1000' )
yeo7_lh <- as.numeric(str_split(lh$label_names, "_", simplify = TRUE)[,2])
yeo7_lh[is.na(yeo7_lh)] <- 0
rh <- subject.annot(subjects_dir, 'fsaverage6', 'rh','Yeo2011_7Networks_N1000' )
yeo7_rh <- as.numeric(str_split(rh$label_names, "_", simplify = TRUE)[,2])
yeo7_rh[is.na(yeo7_rh)] <- 0

yeo7 <- c(yeo7_lh, yeo7_rh)
yeo_sys3 <- ifelse(yeo7==3,1,0)

makecmap_options=list('colFn'=colorRampPalette(c("lightgray",  "#0A9045")))
rglactions=list("snapshot_png"="~/Downloads/yeo_sys3_green_small.png")
vis.data.on.subject(subjects_dir, 'fsaverage6',yeo_sys3[1:40962], yeo_sys3[40963:81924], 
                    "inflated", makecmap_options = makecmap_options, views="t4", rgloptions = rgloptions, rglactions = rglactions)

#GIF
rglactions=list("movie"="yeo_sys3_green")
rgloptions=list("windowRect"=c(50, 50, 50, 50));  
vis.data.on.subject(subjects_dir, 'fsaverage6',yeo_sys3[1:40962], NULL, 
                    "inflated", makecmap_options = makecmap_options, views="sr",rglactions = rglactions, rgloptions = rgloptions)


# Plots -------------------------------------------------------------------
pipeline='gsr_spkreg_fd0.5dvars1.75_drpvls'
pipeline='nogsr_spkreg_fd0.5dvars1.75_drpvls'
pipeline="nogsr_spkreg_fd1.25dvars2_drpvls"
network_dir=paste0("~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/",pipeline)

load(paste0(network_dir,"/CBPD_n90_schaefer400_allruns.Rdata"))
main_unique_lowmotion <- dplyr::filter(main_unique, fd_mean_avg < 0.5)

# Aces effects on attn nets 

#3to3 and 4to4, effects of SES or child ACES
sys3to3_aces <- lm(sys3to3~age_scan+male+fd_mean_avg+avgweight+totalSizet+childaces_sum, data=main_unique_lowmotion)
summary(sys3to3_aces)

sys4to4_aces <- lm(sys4to4~age_scan+male+fd_mean_avg+avgweight+totalSizet+childaces_sum, data=main_unique_lowmotion)
summary(sys4to4_aces)
visreg(sys3to3_aces)
visreg(sys4to4_aces)

sys3to3_aces <- lm(sys3to3~age_scan+male+fd_mean_avg+avgweight+totalSizet+childaces_sum+ses_composite, data=main_unique_lowmotion)
summary(sys3to3_aces)

#Plot
visreg(sys3to3_aces, "childaces_sum", alpha = NULL, main="", xlab="Number of ACES", 
       ylab="Adj dorsal attention network connectivity",line.par = list(col = "orange"))

#Interaction
sys3to3_aces <- lm(sys3to3~age_scan*childaces_sum++male+fd_mean_avg+avgweight+totalSizet+ses_composite, data=main_unique)
summary(sys3to3_aces)
