library(dplyr)
library(stats)
library(parallel)
library(lm.beta)
library(packrat)
library(summarytools)
library(visreg)
library(mgcv)
library(RLRsim)
library(GGally)
library(stringr)
library(tidyverse)
library(data.table)
library(bigmemory)
library(biglm)
library(effectsize)
library(car)
library(lsa)
library(papaja)
library(psych)
options(scipen = 999)

# Setup -------------------------------------------------------------------
library(fsbrain) #this may not work if you're editing the script directly on the cluster...
library(freesurferformats)
#rearrange the order of the brains in the T9 view of fsbrain
source("~/Documents/tools/fsbrain_fix_t9.R")
environment(brainview.t9) <- asNamespace('fsbrain')
assignInNamespace("brainview.t9", brainview.t9, ns = "fsbrain")
subjects_dir = "/Users/utooley/Documents/tools/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/";
subject_id = 'fsaverage';       # for functions which use one subject only
atlas='Schaefer2018_400Parcels_7Networks_order'

#Get the Schaefer atlas region names
schaefer_atlas_region_names_lh = get.atlas.region.names(atlas, template_subjects_dir = subjects_dir,template_subject=subject_id, hemi="lh");
schaefer_atlas_region_names_rh = get.atlas.region.names(atlas, template_subjects_dir = subjects_dir,template_subject=subject_id, hemi="rh");

output_image_directory="~/Documents/projects/in_progress/within_between_network_conn_CBPD/output/figures/brains/"
output_directory="~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/ms_data/"

#Set RGL defaults, larger window and how to make snapshots!
rgloptions=list("windowRect"=c(50,50,1000,1000));

# Main analyses -----------------------------------------------------------
pipeline='gsr_spkreg_fd0.5dvars1.75_drpvls'
network_dir=paste0("~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/",pipeline)
load(paste0(network_dir,"/CBPD_n92_schaefer400_allruns.Rdata"))

## Demographics
summary(main_unique$age_scan)
sum(main_unique$race_black)/dim(main_unique)[1]
sum(main_unique$race_white)/dim(main_unique)[1]
sum(main_unique$race_asian)/dim(main_unique)[1]
(sum(main_unique$race_other)+sum(main_unique$race_americanindian)+sum(main_unique$race_hawaiian))/dim(main_unique)[1]
sum(main_unique$ethnicity=="Hispanic or Latino", na.rm = T)/dim(main_unique)[1]
table(main_unique$male)/dim(main_unique)[1]
table(main_unique$avg_parentedu)/dim(main_unique)[1]
describe(main_unique$avg_parentedu)
table(main_unique$income_median)/dim(main_unique)[1]
describe(main_unique$income_median)
#Scan characteristics 
x <- left_join(main_unique,main_filt, by = c("ID"))
x %>% group_by(ID) %>% summarise(count1=n()) %>% count(count1)

#Figure 1
#Segregation and age
outcomes <- main_unique %>% ungroup() %>% dplyr::select(mean_within_sys,mean_between_sys,system_segreg, modul_avg, part_coef,avgclustco_both) %>% names()
for (outcome in outcomes){
  model<-paste0("lm_",outcome,"_age")
  formula<-formula(paste0(outcome, '~age_scan+male+fd_mean_avg+avgweight+totalSizet'))
  assign(model, lm(formula, data=main_unique))
  output <- paste0("summary_",model)
  assign(output, apa_print(Anova(get(model)), es="pes",mse=FALSE))
  eta <- paste0("eta_",model)
  assign(eta, eta_squared(Anova(get(model)), ci=0.95))
  cohensf <- paste0("cohensf_",model)
  assign(cohensf, cohens_f(Anova(get(model)), ci=0.95))
  beta <- paste0("beta_",model)
  assign(beta, lm.beta(get(model)))
}

## Edge-level analyses
load("~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/edgewise_age_effects_all_cov_n92.Rdata")
edge_age_pvals_mat <- matrix(nrow = 400, ncol=400)
edge_age_pvals_mat[lower.tri(edge_age_pvals_mat, diag=FALSE)] <- edgewise_age_pvals_fdr[,1] #take uncorrected p-values for now
#copy lower triangle to upper triangle
edge_age_pvals_mat[upper.tri(edge_age_pvals_mat)] <- t(edge_age_pvals_mat)[upper.tri(edge_age_pvals_mat)]
edge_age_betas_mat <- matrix(nrow = 400, ncol=400)
edge_age_betas_mat[lower.tri(edge_age_betas_mat, diag=FALSE)] <- edgewise_age_betas
edge_age_betas_mat[upper.tri( edge_age_betas_mat)] <- t(edge_age_betas_mat)[upper.tri( edge_age_betas_mat)]

## Plot on brain
pvalues=c(0.01,0.001,0.0001,0.00001)
pvalue=0.001
#Use different pval thresholds and plot parcels that have age-nonlin-sig edges at that threshold
pos_edges <- edge_age_pvals_mat;pos_edges[edge_age_betas_mat<0] <- NA
neg_edges <- edge_age_pvals_mat;neg_edges[edge_age_betas_mat>0] <- NA
for (pvalue in pvalues){
  indices <- which(neg_edges<pvalue, arr.ind = T) #indices of edges that have age effects  below a pval
  l <- data.frame(table(indices)) #for each parcel, how many age-significant edges does it have at a given pval thresh
  values <- vector(mode = "double", 400)
  for (i in 1:400){
    values[as.numeric(as.character(l$indices[i]))] <- l$Freq[i]#replace the indices in values with the num of sig edges from that parcel in the table
  }
  print(max(values))
  values <- values/2 #because we indexed the full matrix, there are duplicates for each edge
  num_of_sig_age_edges_lh=as.list(setNames(c(0, values[1:200]), schaefer_atlas_region_names_lh))
  num_of_sig_age_edges_rh=as.list(setNames(c(0, values[201:400]), schaefer_atlas_region_names_rh))
  #colormap
  colFn_diverging = colorRampPalette(c("white","#3602D9"));makecmap_options=list('colFn'=colFn_diverging) #range= c(0,10) 
  rglactions=list("snapshot_png"=paste0(output_image_directory,"age_pvals_neg_", pvalue,"regions_colorbar_varying.png")) #purple="#7502E3", red=#E0011C, blue=#3602D9
  vis.region.values.on.subject(subjects_dir, 'fsaverage6', 'Schaefer2018_400Parcels_7Networks_order',  num_of_sig_age_edges_lh, 
                               num_of_sig_age_edges_rh, makecmap_options = makecmap_options, "inflated", views="t4", draw_colorbar = T, rgloptions = rgloptions, rglactions = rglactions)
}

#Nodal age effects
n92_all_edges<- readRDS("~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/n92_all_edges.RData");n92_all_edges <- n92_all_edges[,-1]
#get nodal stregth for each node for each subject
n92_nodal_strength <- apply(n92_all_edges,1, function (x){
  goback <- matrix(nrow = 400, ncol=400)
  goback[lower.tri(goback, diag=FALSE)] <-x
  goback[upper.tri(goback)] <- t(goback)[upper.tri(goback)]
  return(colMeans(goback, na.rm = T))
})
n92_nodal_strength <- t(n92_nodal_strength)
n92_nodal_strength_data <- data.frame(main_unique,n92_nodal_strength)
parcel_sd_pvals<- lapply(names(n92_nodal_strength_data[,1646:2045]), function(x) { summary(lm(as.formula(paste0(x,"~age_scan+male+fd_mean_avg+avgweight+totalSizet")), data=n92_nodal_strength_data))$coef[2,4]})
parcel_sd_pvals <- unlist(parcel_sd_pvals)
parcel_sd_pvals_fdr <- cbind(parcel_sd_pvals,p.adjust(parcel_sd_pvals,method = "fdr"))
#get age betas
parcel_sd_betas<- lapply(names(n92_nodal_strength_data[,1646:2045]), function(x) { lm.beta(lm(as.formula(paste0(x,"~age_scan+male+fd_mean_avg+avgweight+totalSizet")), data=n92_nodal_strength_data))$standardized.coefficients[[2]]})
parcel_sd_betas <- unlist(parcel_sd_betas)
to_plot=ifelse(parcel_sd_pvals_fdr[,2] <0.05,parcel_sd_betas,NA)
lh=as.list(setNames(c(NA,to_plot[1:200]), schaefer_atlas_region_names_lh)) #this has the medial wall in it.
rh=as.list(setNames(c(NA,to_plot[201:400]), schaefer_atlas_region_names_rh))
colFn_diverging = colorRampPalette(rev(RColorBrewer::brewer.pal(9, name="RdBu")));makecmap_options=list('colFn'=colFn_diverging, symm=T)
rglactions=list("snapshot_png"=paste0(output_image_directory,"age_nodal_strength_one_hemi.png"), 'shift_hemis_apart'=TRUE)
#To plot a single hemisphere at an angle, see this gist: https://gist.github.com/dfsp-spirit/d63a6b56f2c11c92086a81911138b453
#must call hook and add the option first (first number is angle to rotate, then x-y-z of axis of rotation. Z is up-down, y is front-back, x is L-R)
hook = function(cm){ for(mesh_idx in seq(length(cm))) { cat(sprintf("Rotating mesh %d.\n", mesh_idx)); cm[[mesh_idx]]$mesh = rgl::rotate3d(cm[[mesh_idx]]$mesh, -pi*0.89, 0,0.3,1);  }; return(cm); }
options('fsbrain.callback_hook_interactive_setup_mesh_before_render' = hook);
vis.region.values.on.subject(subjects_dir, 'fsaverage6', 'Schaefer2018_400Parcels_7Networks_order',  lh, 
                             rh, makecmap_options = makecmap_options, "inflated", views="t4", draw_colorbar = T, rgloptions = rgloptions, rglactions = rglactions)
#change it back
options('fsbrain.callback_hook_interactive_setup_mesh_before_render' = NULL);
nodes <- which(parcel_sd_pvals_fdr[,2] <0.05);nodes <- names(n92_nodal_strength_data[,1645+nodes])
##1 is OFC, number 2 is mPFC
summary(lm(get(nodes[1])~age_scan+male+fd_mean_avg+avgweight+totalSizet, data=n92_nodal_strength_data))
visreg(lm(get(nodes[1])~age_scan+male+fd_mean_avg+avgweight+totalSizet, data=n92_nodal_strength_data))
lm.beta(lm(get(nodes[1])~age_scan+male+fd_mean_avg+avgweight+totalSizet, data=n92_nodal_strength_data))
summary(lm(get(nodes[2])~age_scan+male+fd_mean_avg+avgweight+totalSizet, data=n92_nodal_strength_data))
visreg(lm(get(nodes[2])~age_scan+male+fd_mean_avg+avgweight+totalSizet, data=n92_nodal_strength_data))
lm.beta(lm(get(nodes[2])~age_scan+male+fd_mean_avg+avgweight+totalSizet, data=n92_nodal_strength_data))

#Figure 2
#System-specific effects
networks_age_pvals_fdr #FDR-corrected
outcomes <- main_unique %>% ungroup() %>% dplyr::select(matches("sys.to.")) %>% names()
for (outcome in outcomes){
  model<-paste0("lm_",outcome)
  formula<-formula(paste0(outcome, '~age_scan+male+fd_mean_avg+avgweight+totalSizet'))
  assign(model, lm(formula, data=main_unique))
  output <- paste0("summary_",model)
  assign(output, apa_print(Anova(get(model)),es="pes", mse=FALSE))
  eta <- paste0("eta_",model)
  assign(eta, eta_squared(Anova(get(model)), ci=0.95))
  cohensf <- paste0("cohensf_",model)
  assign(cohensf, cohens_f(Anova(get(model)), ci=0.95))
  beta <- paste0("beta_", model)
  assign(beta, lm.beta(get(model)))
}
summary(lm_sys1to3);visreg(lm_sys1to3,"age_scan");lm.beta(lm_sys1to3)
summary(lm_sys3to7);visreg(lm_sys3to7,"age_scan");lm.beta(lm_sys3to7)
summary(lm_sys4to7);visreg(lm_sys4to7,"age_scan");lm.beta(lm_sys4to7)
summary(lm_sys3to4);visreg(lm_sys3to4,"age_scan");lm.beta(lm_sys3to4)

## Cosine similarity analysis

#read in average matrix
average_FC_mat <- as.matrix(read.csv("~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/n92_average_zscored_FC_matrix.csv"))
#load back in edge age effects
load("~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/edgewise_age_effects_all_cov.Rdata")
#make a matrix
edge_age_pvals_mat <- matrix(nrow = 400, ncol=400);edge_age_pvals_mat[lower.tri(edge_age_pvals_mat, diag=FALSE)] <- edgewise_age_pvals_fdr[,1] #take uncorrected p-values for now
edge_age_pvals_mat[upper.tri(edge_age_pvals_mat)] <- t(edge_age_pvals_mat)[upper.tri(edge_age_pvals_mat)]
#put betas in matrix
edge_age_betas_mat <- matrix(nrow = 400, ncol=400);edge_age_betas_mat[lower.tri(edge_age_betas_mat, diag=FALSE)] <- edgewise_age_betas #take uncorrected p-values for now
edge_age_betas_mat[upper.tri(edge_age_betas_mat)] <- t(edge_age_betas_mat)[upper.tri(edge_age_betas_mat)]

#Use different pval thresholds and plot parcels that have age-nonlin-sig edges at that threshold
pos_edges <- edge_age_pvals_mat;pos_edges[edge_age_betas_mat<0] <- NA
neg_edges <- edge_age_pvals_mat;neg_edges[edge_age_betas_mat>0] <- NA

#Use this code to do edges that have positive sig age effects and edges that have negative sig age effects separately
sig_parcel_cosine=vector()
indices <- which(neg_edges<0.001, arr.ind = T) #indices of edges that have age effects below a pval
l <- data.frame(table(indices)) #for each parcel, how many age-significant edges does it have at a given pval thresh--indices gives duplicates!
values <- vector(mode = "double", 400)
for (i in 1:400){
  values[as.numeric(as.character(l$indices[i]))] <- l$Freq[i]#replace the indices in values with the num of sig edges from that parcel in the table
}
for (x in 1:dim(indices)[1]){
  sig_parcel_cosine[x] <- cosine(average_FC_mat[,indices[x,1]], average_FC_mat[,indices[x,2]])#take both columns
}

null_parcel_cosine=vector()
indices <- which(neg_edges<0.001, arr.ind = T) #indices of edges that have age effects  below a pval
indices[,1] <- sample(1:400,dim(indices)[1], replace = T)#replace with random integers 1-400
indices[,2] <- sample(1:400,dim(indices)[1], replace = T)#replace with random integers 1-400
for (x in 1:dim(indices)[1]){
  null_parcel_cosine[x] <- cosine(average_FC_mat[,indices[x,1]], average_FC_mat[,indices[x,2]])
}
save(sig_parcel_cosine, null_parcel_cosine, file="~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/cosine_similarity_neg_edges_0.001_data.Rdata")

load(file="~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/cosine_similarity_pos_edges_0.001_data.Rdata")
#plot the difference
plt <- ggplot() + theme_classic()+
  geom_density(aes(x=null_parcel_cosine, ..count..),fill="purple", alpha=0.2) + 
  geom_density(aes(x=sig_parcel_cosine,..count..),fill="red", alpha=0.2)+ ggtitle(paste0("Connectivity similarity of parcels w/ pos age FC effects at p <",pvalue)) +
  labs(y= "Number of pairs of parcels", x = "Cosine similarity between parcel connectivity")
print(plt)
print(t.test(sig_parcel_cosine,null_parcel_cosine))

# Main: Reasoning analyses ------------------------------------------------
#Figure 3
main_unique$matrix_reasoning_both <- ifelse(is.na(main_unique$wppsi_matrix_valid),main_unique$wisc_matrix_scaled,ifelse(main_unique$wppsi_matrix_valid==0,NA,main_unique$wppsi_matrix_scaled))
main_unique$matrix_reasoning_both_raw <- ifelse(is.na(main_unique$wppsi_matrix_valid),main_unique$wisc_matrix_raw,ifelse(main_unique$wppsi_matrix_valid==0,NA,main_unique$wppsi_matrix_raw))
main_unique$matrix_reasoning_both;main_unique$matrix_reasoning_both_raw

## Whole-brain positive and negative to reasoning
n92_all_edges<- readRDS("~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/n92_all_edges.RData");n92_all_edges <- n92_all_edges[,-1]

# edgewise_mr_pvals<- mclapply(1:79800, function(x) { summary(lm(n92_all_edges[,x]~main_unique$matrix_reasoning_both+main_unique$age_scan+main_unique$male+main_unique$fd_mean_avg+main_unique$avgweight+main_unique$totalSizet))$coef[2,4]}, mc.cores = 4)
# edgewise_mr_pvals <- unlist(edgewise_mr_pvals)
# edgewise_mr_pvals_fdr <- cbind(edgewise_mr_pvals,p.adjust(edgewise_mr_pvals,method = "fdr"))
# #get mr betas
# edgewise_mr_betas<- mclapply(1:79800, function(x) { lm.beta(lm(n92_all_edges[,x]~main_unique$matrix_reasoning_both+main_unique$age_scan+main_unique$male+main_unique$fd_mean_avg+main_unique$avgweight+main_unique$totalSizet))$standardized.coefficients[[2]]}, mc.cores = 4)
# edgewise_mr_betas <- unlist(edgewise_mr_betas)
# #save for reloading in future
# save(edgewise_mr_pvals_fdr, edgewise_mr_betas,file= "~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/edgewise_mr_effects_all_cov_n92.Rdata")

#load back in edge MR effects
load("~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/edgewise_mr_effects_all_cov_n92.Rdata")

#code to read back into a 400 x 400 matrix! This works correctly.
# vec<- mat[lower.tri(mat)]
edge_mr_pvals_mat <- matrix(nrow = 400, ncol=400)
edge_mr_pvals_mat[lower.tri(edge_mr_pvals_mat, diag=FALSE)] <- edgewise_mr_pvals_fdr[,1] #take uncorrected p-values for now
#copy lower triangle to upper triangle
edge_mr_pvals_mat[upper.tri(edge_mr_pvals_mat)] <- t(edge_mr_pvals_mat)[upper.tri(edge_mr_pvals_mat)]
edge_mr_betas_mat <- matrix(nrow = 400, ncol=400)
edge_mr_betas_mat[lower.tri(edge_mr_betas_mat, diag=FALSE)] <- edgewise_mr_betas
edge_mr_betas_mat[upper.tri( edge_mr_betas_mat)] <- t(edge_mr_betas_mat)[upper.tri( edge_mr_betas_mat)]

#pvalues=c(0.01,0.001,0.0001,0.00001)
pvalue=0.001
#Use different pval thresholds and plot parcels that have age-nonlin-sig edges at that threshold
pos_edges <- edge_mr_pvals_mat;pos_edges[edge_mr_betas_mat<0] <- NA
neg_edges <- edge_mr_pvals_mat;neg_edges[edge_mr_betas_mat>0] <- NA
#for (pvalue in pvalues){
  indices <- which(neg_edges<pvalue, arr.ind = T) #indices of edges that have age effects  below a pval
  l <- data.frame(table(indices)) #for each parcel, how many age-significant edges does it have at a given pval thresh
  values <- vector(mode = "double", 400)
  for (i in 1:400){
    values[as.numeric(as.character(l$indices[i]))] <- l$Freq[i]#replace the indices in values with the num of sig edges from that parcel in the table
  }
  print(max(values))
  values <- values/2 #because we indexed the full matrix, there are duplicates for each edge
  num_of_sig_age_edges_lh=as.list(setNames(c(0, values[1:200]), schaefer_atlas_region_names_lh))
  num_of_sig_age_edges_rh=as.list(setNames(c(0, values[201:400]), schaefer_atlas_region_names_rh))
  #colormap
  colFn_diverging = colorRampPalette(c("white","#3602D9"));makecmap_options=list('colFn'=colFn_diverging, range= c(0,10)) 
  rglactions=list("snapshot_png"=paste0(output_image_directory,"mr_neg_pvals_", pvalue,"regions.png"), 'shift_hemis_apart'= T) #purple="#7502E3", red=#E0011C, blue=#3602D9
  view=vis.region.values.on.subject(subjects_dir, 'fsaverage6', 'Schaefer2018_400Parcels_7Networks_order',  num_of_sig_age_edges_lh, 
                               num_of_sig_age_edges_rh, makecmap_options = makecmap_options, "inflated", views="t9", draw_colorbar = T, rgloptions = rgloptions, rglactions = rglactions)
#}  
output_brain_img = "fsbrain_arranged.png";
vislayout.from.coloredmeshes(view, view_angles = 'sd_medial_lh', rglactions = rglactions, output_img = output_brain_img);

#Raincloud plots
source("~/Documents/tools/raincloud.R")
community_assign <- read.delim("~/Documents/tools/parcellations/schaefer400x7CommunityAffiliation.1D", header = F)
communities <- paste0("comm",community_assign$V1) %>% car::recode(.,"'comm1'='Visual'; 'comm2'='Somatomotor' ; 'comm3'='Dorsal Attention';'comm4'='Ventral Attention';'comm5'='Limbic';'comm6'='Frontoparietal';'comm7'='Default'") %>% factor(levels=c("Visual","Somatomotor","Dorsal Attention", "Ventral Attention", "Limbic", "Frontoparietal", "Default"))
edge_mr_pvals_mat <- matrix(nrow = 400, ncol=400);edge_mr_pvals_mat[lower.tri(edge_mr_pvals_mat, diag=FALSE)] <- edgewise_mr_pvals_fdr[,1];edge_mr_pvals_mat[upper.tri(edge_mr_pvals_mat)] <- t(edge_mr_pvals_mat)[upper.tri(edge_mr_pvals_mat)]
edge_mr_betas_mat <- matrix(nrow = 400, ncol=400);edge_mr_betas_mat[lower.tri(edge_mr_betas_mat, diag=FALSE)] <- edgewise_mr_betas; edge_mr_betas_mat[upper.tri( edge_mr_betas_mat)] <- t(edge_mr_betas_mat)[upper.tri( edge_mr_betas_mat)]
#Use different pval thresholds and plot parcels that have age-nonlin-sig edges at that threshold
pos_edges <- edge_mr_pvals_mat;pos_edges[edge_mr_betas_mat<0] <- NA
neg_edges <- edge_mr_pvals_mat;neg_edges[edge_mr_betas_mat>0] <- NA
pvalue=0.001
indices <- which(neg_edges<pvalue, arr.ind = T) #indices of edges that have age effects  below a pval
l <- data.frame(table(indices)) #for each parcel, how many age-significant edges does it have at a given pval thresh
values <- vector(mode = "double", 400)
for (i in 1:400){
  values[as.numeric(as.character(l$indices[i]))] <- l$Freq[i]#replace the indices in values with the num of sig edges from that parcel in the table
}
values <- values/2 #because we indexed twice
num_edges <- data.frame(factor(communities, ),values);colnames(num_edges) <- c("community","data"); 
longdata = num_edges %>% dplyr::group_by(community) %>% mutate(med = median(data))

g <- ggplot(data = longdata, aes(y = data, x = community, fill=community)) +
  #geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .5) +
  geom_point(aes(y = data, color = community), position = position_jitter(width = .15), size = .3, alpha = 0.5) +geom_boxplot(width = .3, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  ylim(c(0,10)) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  # scale_color_manual(values=c("#4669dc","pink","purple")) +
  #scale_color_distiller(palette="Blues", direction = 1) +
  #scale_fill_distiller(palette="Blues", direction = 1) +
  scale_color_manual(values=c("#7B287E", "#5CA1C8", "#0A9045","#C33AF8", "#B6C988", "#EF9C23", "#E34A53")) +
  scale_fill_manual(values= c("#7B287E", "#5CA1C8", "#0A9045","#C33AF8", "#B6C988", "#EF9C23", "#E34A53")) +
  #scale_color_gradient2(low="#FFFFFF",mid= "#FB6A4A", high="#67000D", midpoint= 0.475, aesthetics = c("color", "fill")) +
  theme_bw() +
  raincloud_theme
print(g)

kruskal.test(data~community, data = num_edges)
pairwise.wilcox.test(num_edges$data, num_edges$community,
                     p.adjust.method = "fdr")
num_edges$community <- relevel(num_edges$community, "Visual")
anova(lm(data~community,data=num_edges))
TukeyHSD(aov(lm(data~community,data=num_edges)))

#Examine betas instead
betas <- list()
edge_mr_betas_mat_thresh <- edge_mr_betas_mat
links_to_show=edge_mr_pvals_mat<0.001;edge_mr_betas_mat_thresh[links_to_show==FALSE] <- NA
rownames(edge_mr_betas_mat_thresh) <- communities;colnames(edge_mr_betas_mat_thresh) <- communities
for (j in 1:7){
  set=unique(communities)[j]
  #take the rows from vis and the columns from SM (and vice-versa?)
  betas[[set]] <- na.omit(c(edge_mr_betas_mat_thresh[rownames(edge_mr_betas_mat_thresh)==set,]))
}
betas<- lapply(betas, function (x) {length (x) <- 18953;x}) #make the all the same length and pad with NAS
betas <- do.call("cbind", betas) #can also rbind if want rows
colnames(betas) <- unique(communities);longdata <- melt(betas);colnames(longdata) <- c("x","community","value")
kruskal.test(value~community, data = longdata)
pairwise.wilcox.test(longdata$value, longdata$community,
                     p.adjust.method = "fdr")
summary(lm(value~community,data=longdata))
TukeyHSD(aov(lm(value~community,data=longdata)))

#Test measures of segregation with reasoning
outcomes <- main_unique %>% ungroup() %>% dplyr::select(mean_within_sys,mean_between_sys,system_segreg, modul_avg, part_coef,avgclustco_both) %>% names()
for (outcome in outcomes){
  model<-paste0("lm_",outcome,"_age")
  formula<-formula(paste0('matrix_reasoning_both~age_scan+male+fd_mean_avg+avgweight+totalSizet+',outcome))
  assign(model, lm(formula, data=main_unique))
  output <- paste0("summary_",model)
  assign(output, apa_print(Anova(get(model)), es="pes",mse=FALSE))
  eta <- paste0("eta_",model)
  assign(eta, eta_squared(Anova(get(model)), ci=0.95))
  cohensf <- paste0("cohensf_",model)
  assign(cohensf, cohens_f(Anova(get(model)), ci=0.95))
  beta <- paste0("beta_",model)
  assign(beta, lm.beta(get(model)))
  print(outcome)
  print(summary(get(model)))
}
#Test all 4 age-sig systems with Reasoning
lm_reasoning_sys1to3<- lm(matrix_reasoning_both~age_scan+male+fd_mean_avg+avgweight+totalSizet+sys1to3, data=main_unique)
lm_reasoning_sys3to7<- lm(matrix_reasoning_both~age_scan+male+fd_mean_avg+avgweight+totalSizet+sys3to7, data=main_unique)
lm_reasoning_sys4to7<- lm(matrix_reasoning_both~age_scan+male+fd_mean_avg+avgweight+totalSizet+sys4to7, data=main_unique)
summary_lm_reasoning_sys1to3 <- apa_print(Anova(lm_reasoning_sys1to3),es="pes",mse=F);beta_lm_reasoning_sys1to3 <- lm.beta(lm_reasoning_sys1to3)
eta_lm_reasoning_sys1to3 <- eta_squared(Anova(lm_reasoning_sys1to3), ci=0.95);cohensf_lm_reasoning_sys1to3 <- cohens_f(Anova(lm_reasoning_sys1to3), ci=0.95)
summary_lm_reasoning_sys3to7 <- apa_print(Anova(lm_reasoning_sys3to7),es="pes",mse=F);beta_lm_reasoning_sys3to7 <- lm.beta(lm_reasoning_sys3to7)
eta_lm_reasoning_sys3to7 <- eta_squared(Anova(lm_reasoning_sys3to7), ci=0.95);cohensf_lm_reasoning_sys3to7 <- cohens_f(Anova(lm_reasoning_sys3to7), ci=0.95)
summary_lm_reasoning_sys4to7 <- apa_print(Anova(lm_reasoning_sys4to7),es="pes",mse=F);beta_lm_reasoning_sys4to7 <- lm.beta(lm_reasoning_sys4to7)
eta_lm_reasoning_sys4to7 <- eta_squared(Anova(lm_reasoning_sys4to7), ci=0.95);cohensf_lm_reasoning_sys4to7 <- cohens_f(Anova(lm_reasoning_sys4to7), ci=0.95)

visreg(lm_sys1to3, "sys1to3", main= "Matrix Reasoning (scaled score)", ylab="Scaled score", xlab="VIS to DA connectivity")
visreg(lm_sys3to7, "sys3to7", main= "Matrix Reasoning (scaled score)", ylab="Scaled score", xlab="DM to DA connectivity")
visreg(lm_sys4to7)

#Test reasoning with FP system
lm_reasoning_sys6to6<- lm(matrix_reasoning_both~age_scan+male+fd_mean_avg+avgweight+totalSizet+sys6to6, data=main_unique)
summary_lm_reasoning_sys6to6 <- apa_print(Anova(lm_reasoning_sys6to6),es="pes",mse=F);beta_lm_reasoning_sys6to6 <- lm.beta(lm_reasoning_sys6to6)
eta_lm_reasoning_sys6to6 <- eta_squared(Anova(lm_reasoning_sys6to6), ci=0.95);cohensf_lm_reasoning_sys6to6 <- cohens_f(Anova(lm_reasoning_sys6to6), ci=0.95)
summary_lm_reasoning_sys6to6
summary(lm_reasoning_sys6to6)

## Play
lm_sys3to7<- lm(litnum_freq_play_avg~age_scan+male+fd_mean_avg+avgweight+totalSizet+sys3to7, data=main_unique)
lm_sys4to7<- lm(sys4to7~age_scan+male+fd_mean_avg+avgweight+totalSizet+litnum_freq_play_avg, data=main_unique)
summary(lm(litnum_freq_play_avg~age_behav.x+male+ses_composite, data=main_unique));visreg(lm(litnum_freq_play_avg~age_behav.x+male+ses_composite, data=main_unique))
summary(lm_sys3to7);lm.beta(lm_sys3to7) #Not significant
summary(lm_sys4to7)
visreg(lm_sys3to7, "sys3to7", main="Frequency of play",ylab="", xlab="DMN to DAN")
visreg(lm_sys4to7, "litnum_freq_play_avg", main="DMN to VAN",ylab="", xlab="litnum_freq_play_avg")

# Replication in alternate parcellation --------
pipeline='gsr_spkreg_fd0.5dvars1.75_drpvls'
network_dir=paste0("~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/",pipeline)

#trace and edit this function to remove MSE
trace(papaja:::apa_print.anova, edit=TRUE)

load(paste0(network_dir,"/CBPD_n92_schaefer400_allruns.Rdata"))
#Use main_replicate_unique instead
lm_replicate_schaefer_within_sys_age <- lm(mean_within_sys~age_scan+male+fd_mean_avg+avgweight+totalSizet, data=main_replicate_unique)
summary_replicate_schaefer_within_sys_age <- apa_print(Anova(lm_replicate_schaefer_within_sys_age));lm_eta_replicate_schaefer_within_sys_age <- eta_squared(Anova(lm_replicate_schaefer_within_sys_age), ci=0.95)
outcomes <- main_replicate_unique %>% ungroup() %>% dplyr::select(mean_within_sys,mean_between_sys,system_segreg, modul_avg, part_coef,avgclustco_both) %>% names()
for (outcome in outcomes){
  model<-paste0("lm_replicate_schaefer_",outcome,"_age")
  formula<-formula(paste0(outcome, '~age_scan+male+fd_mean_avg+avgweight+totalSizet'))
  assign(model, lm(formula, data=main_replicate_unique))
  output <- paste0("summary_",model)
  assign(output, apa_print(Anova(get(model)), es="pes",mse=F))
  eta <- paste0("eta_",model)
  assign(eta, eta_squared(Anova(get(model)), ci=0.95))
  cohensf <- paste0("cohensf_",model)
  assign(cohensf, cohens_f(Anova(get(model)), ci=0.95))
  beta <- paste0("beta",model)
  assign(beta, lm.beta(get(model)))
  print(outcome)
  print(get(beta))
  print(summary(get(model)))
}
visreg(lm_replicate_schaefer_between_sys_age,"age_scan", main="Mean between-system connectivity", ylab="", xlab="Age at scan (years)")
visreg(lm_replicate_schaefer_within_sys_age,"age_scan", main= "Mean within-system connectivity", ylab="", xlab="Age at scan (years)")
visreg(lm_replicate_schaefer_segreg_age,"age_scan", main="System segregation", ylab="", xlab="Age at scan (years)")
visreg(lm_replicate_schaefer_modul_avg_age,"age_scan", main="Modularity quality index (Q)", ylab="", xlab="Age at scan (years)")
visreg(lm_replicate_schaefer_part_coef_age,"age_scan", main="Average participation coefficient", ylab="", xlab="Age at scan (years)")
visreg(lm_replicate_schaefer_clust_co_age,"age_scan", main="Average clustering coefficient", ylab="", xlab="Age at scan (years)")

networks_age_pvals_fdr_replicate #FDR-corrected
outcomes <- main_replicate_unique %>% ungroup() %>% dplyr::select(mean_within_sys,sys1to3,sys3to4,sys3to7,sys4to7) %>% names()
for (outcome in outcomes){
  model<-paste0("lm_replicate_schaefer_",outcome)
  formula<-formula(paste0(outcome, '~age_scan+male+fd_mean_avg+avgweight+totalSizet'))
  assign(model, lm(formula, data=main_replicate_unique))
  output <- paste0("summary_",model)
  assign(output, apa_print(Anova(get(model)),es="pes",mse=F))
  eta <- paste0("eta_",model)
  assign(eta, eta_squared(Anova(get(model)), ci=0.95))
  cohensf <- paste0("cohensf_",model)
  assign(cohensf, cohens_f(Anova(get(model)), ci=0.95))
  beta <- paste0("beta_", model)
  assign(beta, lm.beta(get(model)))
}
summary(lm_replicate_schaefer_sys1to3);visreg(lm_replicate_schaefer_sys1to3,"age_scan");lm.beta(lm_replicate_schaefer_sys1to3)
summary(lm_replicate_schaefer_sys3to7);visreg(lm_replicate_schaefer_sys3to7,"age_scan");lm.beta(lm_replicate_schaefer_sys3to7)
summary(lm_replicate_schaefer_sys4to7);visreg(lm_replicate_schaefer_sys4to7,"age_scan");lm.beta(lm_replicate_schaefer_sys4to7)
summary(lm_replicate_schaefer_sys3to4);visreg(lm_replicate_schaefer_sys3to4,"age_scan");lm.beta(lm_replicate_schaefer_sys3to4)

## Edge-level analyses
pipeline='gsr_spkreg_fd0.5dvars1.75_drpvls'
edges_dir=paste0("/cbica/home/tooleyu/projects/in_progress/within_between_network_conn_CBPD/data/imageData/",pipeline,"/schaefer200zNetworks_avg/")
#make a huge dataframe of subjects x edges

# files <- data.frame(list.files(edges_dir)) %>% filter( list.files.edges_dir. %in% paste0(main_unique$ID, "_schaefer200MNI_zavgnetwork.txt")); files$list.files.edges_dir. <- paste0(edges_dir,files$list.files.edges_dir.)
# vect <- matrix(nrow = dim(main_unique)[1], ncol = 19900)
# for (i in 1:dim(main_unique)[1]){
# mat <- data.frame(read_csv(file =paste0(edges_dir,main_unique$ID[i], "_schaefer200MNI_zavgnetwork.txt"), col_names = F))
# vect[i,1:19900]<-mat[lower.tri(mat)] #read out the lower triangle from matrix by row
# }
# 
# #code to read back into a 200 x 200 matrix! This works correctly.
# vec<- mat[lower.tri(mat)]
# goback <- matrix(nrow = 200, ncol=200)
# goback[lower.tri(goback, diag=FALSE)] <- vec
# goback <- t(goback)
# goback
# 
saveRDS(vect,"~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/n92_schaefer200_replicate_all_edges.RData")

#Read back in the matrix of edges, so you don't create every time.
#n90_all_edges<- readRDS("~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/n90_all_edges.Rds");n90_all_edges <- n90_all_edges[,-1]
n92_all_edges<- readRDS("~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/n92_schaefer200_replicate_all_edges.RData")

#multi-core apply the linear model across the matrix of edges, get age pvals 
edgewise_age_pvals<- mclapply(1:19900, function(x) { summary(lm(n92_all_edges[,x]~main_unique$age_scan+main_unique$male+main_unique$fd_mean_avg+main_unique$avgweight+main_unique$totalSizet))$coef[2,4]}, mc.cores = 4)
edgewise_age_pvals <- unlist(edgewise_age_pvals)
edgewise_age_pvals_fdr <- cbind(edgewise_age_pvals,p.adjust(edgewise_age_pvals,method = "fdr"))
#get age betas
edgewise_age_betas<- mclapply(1:19900, function(x) { lm.beta(lm(n92_all_edges[,x]~main_unique$age_scan+main_unique$male+main_unique$fd_mean_avg+main_unique$avgweight+main_unique$totalSizet))$standardized.coefficients[[2]]}, mc.cores = 4)
edgewise_age_betas <- unlist(edgewise_age_betas)
#save for reloading in future
save(edgewise_age_pvals_fdr, edgewise_age_betas,file= "~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/schaefer200_edgewise_age_effects_all_cov_n92.Rdata")

#load back in edge age effects
load("~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/schaefer200_edgewise_age_effects_all_cov_n92.Rdata")

#code to read back into a 200 x 200 matrix! This works correctly.
# vec<- mat[lower.tri(mat)]
edge_age_pvals_mat <- matrix(nrow = 200, ncol=200)
edge_age_pvals_mat[lower.tri(edge_age_pvals_mat, diag=FALSE)] <- edgewise_age_pvals_fdr[,1] #take uncorrected p-values for now
#copy lower triangle to upper triangle
edge_age_pvals_mat[upper.tri(edge_age_pvals_mat)] <- t(edge_age_pvals_mat)[upper.tri(edge_age_pvals_mat)]
edge_age_betas_mat <- matrix(nrow = 200, ncol=200)
edge_age_betas_mat[lower.tri(edge_age_betas_mat, diag=FALSE)] <- edgewise_age_betas
edge_age_betas_mat[upper.tri( edge_age_betas_mat)] <- t(edge_age_betas_mat)[upper.tri( edge_age_betas_mat)]

## Plot on brain
schaefer_atlas_region_names_lh = get.atlas.region.names('Schaefer2018_200Parcels_7Networks_order', template_subjects_dir = subjects_dir,template_subject=subject_id, hemi="lh");
schaefer_atlas_region_names_rh = get.atlas.region.names('Schaefer2018_200Parcels_7Networks_order', template_subjects_dir = subjects_dir,template_subject=subject_id, hemi="rh");
output_image_directory="/Users/utooley/Dropbox/projects/in_progress/within_between_network_conn_CBPD/output/figures/replicate_schaefer200/"
pvalues=c(0.01,0.001,0.0001,0.00001)
#Use different pval thresholds and plot parcels that have age-nonlin-sig edges at that threshold
pos_edges <- edge_age_pvals_mat;pos_edges[edge_age_betas_mat<0] <- NA
neg_edges <- edge_age_pvals_mat;neg_edges[edge_age_betas_mat>0] <- NA
for (pvalue in pvalues){
  indices <- which(edge_age_pvals_mat<pvalue, arr.ind = T) #indices of edges that have age effects  below a pval
  l <- data.frame(table(indices)) #for each parcel, how many age-significant edges does it have at a given pval thresh
  values <- vector(mode = "double", 400)
  for (i in 1:200){
    values[as.numeric(as.character(l$indices[i]))] <- l$Freq[i]#replace the indices in values with the num of sig edges from that parcel in the table
  }
  print(max(values))
  #values <- values/2 #because we indexed the full matrix, there are duplicates for each edge
  num_of_sig_age_edges_lh=as.list(setNames(c(0, values[1:100]), schaefer_atlas_region_names_lh))
  num_of_sig_age_edges_rh=as.list(setNames(c(0, values[101:200]), schaefer_atlas_region_names_rh))
  #colormap
  colFn_diverging = colorRampPalette(c("white","#7502E3"));makecmap_options=list('colFn'=colFn_diverging, range= c(0,10)) 
  rglactions=list("snapshot_png"=paste0(output_image_directory,"age_pvals_", pvalue,"regions.png")) #purple="#7502E3", red=#E0011C, blue=#3602D9
  vis.region.values.on.subject(subjects_dir, 'fsaverage6', 'Schaefer2018_400Parcels_7Networks_order',  num_of_sig_age_edges_lh, 
                               num_of_sig_age_edges_rh, makecmap_options = makecmap_options, "inflated", views="t4", draw_colorbar = T, rgloptions = rgloptions, rglactions = rglactions)
}
#get edge matrix, make average nodal strength
n92_all_edges<- readRDS("~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/n92_schaefer200_replicate_all_edges.RData")

#get nodal stregth for each node for each subject
n92_nodal_strength <- apply(n92_all_edges,1, function (x){
  goback <- matrix(nrow = 200, ncol=200)
  goback[lower.tri(goback, diag=FALSE)] <-x
  goback[upper.tri(goback)] <- t(goback)[upper.tri(goback)]
  return(colMeans(goback, na.rm = T))
})
n92_nodal_strength <- t(n92_nodal_strength)
n92_nodal_strength_data <- data.frame(main_replicate_unique,n92_nodal_strength)

subjectwise_metric <- data.frame(n92_nodal_strength,rowMeans(n92_nodal_strength, na.rm = T)) %>% rename(global_mean_metric=rowMeans.n92_nodal_strength..na.rm...T.)
hist(subjectwise_metric$global_mean_metric)
#take out the 2 outliers
#subjectwise_metric <- subjectwise_metric %>% filter(.,global_mean_metric>1)
subjectwise_metric$ID <- main_replicate_unique$ID
subjectwise_metric<- left_join(main_replicate_unique,subjectwise_metric, by= "ID")
summary(lm(global_mean_metric~age_scan+male+fd_mean_avg+totalSizet, data=subjectwise_metric))
visreg(lm(global_mean_metric~age_scan+male+fd_mean_avg+avgweight+totalSizet, data=subjectwise_metric))
exactRLRT(gamm(global_mean_metric~s(age_scan)+male+fd_mean_avg+avgweight+totalSizet, data=subjectwise_metric, method = "REML")$lme)
visreg(gam(global_mean_metric~s(age_scan)+male+fd_mean_avg+avgweight+totalSizet, data=subjectwise_metric, method = "REML"))
# #take out rows with NAs/Inf at the parcel level
# subjectwise_metric %>% ungroup() %>%  select(V1:V100) %>% summary() mutate(global_mean_metric=rowMeans(as.matrix(.))) %>%  select(global_mean_metric)
# mutate(global_mean_metric=rowMedians(as.matrix(.),na.rm=T))
#parcelwise
parcel_sd_pvals<- lapply(names(subjectwise_metric[,1601:1800]), function(x) { summary(lm(as.formula(paste0(x,"~age_scan+male+fd_mean_avg+avgweight+totalSizet")), data=subjectwise_metric))$coef[2,4]})
parcel_sd_pvals <- unlist(parcel_sd_pvals)
parcel_sd_pvals_fdr <- cbind(parcel_sd_pvals,p.adjust(parcel_sd_pvals,method = "fdr"))
#get age betas
parcel_sd_betas<- lapply(names(subjectwise_metric[,1601:1800]), function(x) { lm.beta(lm(as.formula(paste0(x,"~age_scan+male+fd_mean_avg+avgweight+totalSizet")), data=subjectwise_metric))$standardized.coefficients[[2]]})
parcel_sd_betas <- unlist(parcel_sd_betas)
name="parcel_age_node_strength_pvals"
to_plot=parcel_sd_pvals_fdr[,1] #Change this to plot different things
to_plot=ifelse(parcel_sd_pvals_fdr[,1] <0.01,parcel_sd_betas,NA)

lh=as.list(setNames(c(NA,to_plot[1:100]), schaefer_atlas_region_names_lh)) #this has the medial wall in it.
rh=as.list(setNames(c(NA,to_plot[101:200]), schaefer_atlas_region_names_rh))
#colormap
#colFn_diverging = colorRampPalette(rev(c("white","white","pink","red")));makecmap_options=list('colFn'=colFn_diverging)
colFn_diverging = colorRampPalette(rev(RColorBrewer::brewer.pal(9, name="RdBu")));makecmap_options=list('colFn'=colFn_diverging, symm=T)
rglactions=list("snapshot_png"=paste0(output_image_directory,name,".png"), 'shift_hemis_apart'=TRUE)
vis.region.values.on.subject(subjects_dir, 'fsaverage6', 'Schaefer2018_400Parcels_7Networks_order',  lh, 
                             rh, makecmap_options = makecmap_options, "inflated", views="t9", draw_colorbar = T, rgloptions = rgloptions, rglactions = rglactions)
## Reasoning
#Test all 3 systems with Reasoning
main_replicate_unique$matrix_reasoning_both <- ifelse(is.na(main_replicate_unique$wppsi_matrix_valid),main_replicate_unique$wisc_matrix_scaled,ifelse(main_replicate_unique$wppsi_matrix_valid==0,NA,main_replicate_unique$wppsi_matrix_scaled))
main_replicate_unique$matrix_reasoning_both_raw <- ifelse(is.na(main_replicate_unique$wppsi_matrix_valid),main_replicate_unique$wisc_matrix_raw,ifelse(main_replicate_unique$wppsi_matrix_valid==0,NA,main_replicate_unique$wppsi_matrix_raw))
main_replicate_unique$matrix_reasoning_both;main_replicate_unique$matrix_reasoning_both_raw
lm_replicate_schaefer_reasoning_sys1to3<- lm(matrix_reasoning_both~age_scan+male+fd_mean_avg+avgweight+totalSizet+sys1to3, data=main_replicate_unique)
lm_replicate_schaefer_reasoning_sys3to7<- lm(matrix_reasoning_both~age_scan+male+fd_mean_avg+avgweight+totalSizet+sys3to7, data=main_replicate_unique)
lm_replicate_schaefer_reasoning_sys4to7<- lm(matrix_reasoning_both~age_scan+male+fd_mean_avg+avgweight+totalSizet+sys4to7, data=main_replicate_unique)
summary_lm_replicate_schaefer_reasoning_sys1to3 <- apa_print(Anova(lm_replicate_schaefer_reasoning_sys1to3),es="pes",mse=F);beta_lm_replicate_schaefer_reasoning_sys1to3 <- lm.beta(lm_replicate_schaefer_reasoning_sys1to3)
eta_lm_replicate_schaefer_reasoning_sys1to3 <- eta_squared(Anova(lm_replicate_schaefer_reasoning_sys1to3), ci=0.95);cohensf_lm_replicate_schaefer_reasoning_sys1to3 <- cohens_f(Anova(lm_replicate_schaefer_reasoning_sys1to3), ci=0.95)
summary_lm_replicate_schaefer_reasoning_sys3to7 <- apa_print(Anova(lm_replicate_schaefer_reasoning_sys3to7),es="pes",mse=F);beta_lm_replicate_schaefer_reasoning_sys3to7 <- lm.beta(lm_replicate_schaefer_reasoning_sys3to7)
eta_lm_replicate_schaefer_reasoning_sys3to7 <- eta_squared(Anova(lm_replicate_schaefer_reasoning_sys3to7), ci=0.95);cohensf_lm_replicate_schaefer_reasoning_sys3to7 <- cohens_f(Anova(lm_replicate_schaefer_reasoning_sys3to7), ci=0.95)
summary_lm_replicate_schaefer_reasoning_sys4to7 <- apa_print(Anova(lm_replicate_schaefer_reasoning_sys4to7),es="pes",mse=F);beta_lm_replicate_schaefer_reasoning_sys4to7 <- lm.beta(lm_replicate_schaefer_reasoning_sys4to7)
eta_lm_replicate_schaefer_reasoning_sys4to7 <- eta_squared(Anova(lm_replicate_schaefer_reasoning_sys4to7), ci=0.95);cohensf_lm_replicate_schaefer_reasoning_sys4to7 <- cohens_f(Anova(lm_replicate_schaefer_reasoning_sys4to7), ci=0.95)

visreg(lm_replicate_schaefer_reasoning_sys1to3, "sys1to3", main= "Matrix Reasoning (scaled score)", ylab="Scaled score", xlab="VIS to DA connectivity")
visreg(lm_replicate_schaefer_reasoning_sys3to7, "sys3to7", main= "Matrix Reasoning (scaled score)", ylab="Scaled score", xlab="DM to DA connectivity")
visreg(lm_replicate_schaefer_reasoning_sys4to7)

## Play
lm_replicate_schaefer_play_sys3to7<- lm(litnum_freq_play_avg~age_scan+male+fd_mean_avg+avgweight+totalSizet+sys3to7, data=main_replicate_unique)
lm_replicate_schaefer_play_sys4to7<- lm(litnum_freq_play_avg~age_scan+male+fd_mean_avg+avgweight+totalSizet+sys4to7, data=main_replicate_unique)
summary(lm(litnum_freq_play_avg~age_behav.x+male+ses_composite, data=main_replicate_unique))
#visreg(lm(litnum_freq_play_avg~age_behav.x+male+ses_composite, data=main_replicate_unique))
summary(lm_replicate_schaefer_play_sys3to7);lm.beta(lm_replicate_schaefer_play_sys3to7) #Not significant
summary_lm_replicate_schaefer_reasoning_sys4to7 <- apa_print(Anova(lm_replicate_schaefer_play_sys3to7),es="pes",mse=F);beta_lm_replicate_schaefer_play_sys3to7 <- lm.beta(lm_replicate_schaefer_play_sys3to7)
eta_lm_replicate_schaefer_play_sys3to7 <- eta_squared(Anova(lm_replicate_schaefer_play_sys3to7), ci=0.95);cohensf_lm_replicate_schaefer_play_sys3to7 <- cohens_f(Anova(lm_replicate_schaefer_play_sys3to7), ci=0.95)
summary(lm_replicate_schaefer_play_sys4to7)
visreg(lm_replicate_schaefer_play_sys3to7, "sys3to7", main="Frequency of play",ylab="", xlab="DMN to DAN")
visreg(lm_replicate_schaefer_play_sys4to7, "litnum_freq_play_avg", main="DMN to VAN",ylab="", xlab="litnum_freq_play_avg")

save(list=c("main_replicate_unique", "networks_age_pvals_fdr_replicate", ls(pattern = "summary_lm_replicate_*"), 
            ls(pattern = "eta_lm_replicate_schaefer_*"), ls(pattern = "cohensf_lm_replicate_schaefer_*"), ls(pattern = "beta_lm_replicate_schaefer_*")),
     file = paste0(output_directory,"replicate_schaefer_data.Rdata"))

# Replication in pipeline without GSR --------

pipeline='nogsr_spkreg_fd0.5dvars1.75_drpvls'
pipeline="nogsr_spkreg_fd1.25dvars2_drpvls"
network_dir=paste0("~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/",pipeline)
load(paste0(network_dir,"/CBPD_n92_schaefer400_allruns.Rdata"))

#there is an extra participant here who gets filtered out in the other pipeline
main_unique <- filter(main_unique,ID != "sub-CBPD0018")
lm_within_sys_age <- lm(mean_within_sys~age_scan+male+fd_mean_avg+avgweight+totalSizet, data=main_unique)
summary(lm_within_sys_age)
lm_between_sys_age <- lm(mean_between_sys~age_scan+male+fd_mean_avg+avgweight+totalSizet, data=main_unique)
summary(lm_between_sys_age)
lm_segreg_age<- lm(system_segreg~age_scan+male+fd_mean_avg+avgweight+totalSizet, data=main_unique)
summary(lm_segreg_age)
lm_modul_avg_age <- lm(modul_avg~age_scan+male+fd_mean_avg+avgweight+totalSizet, data=main_unique)
summary(lm_modul_avg_age)
lm_part_coef_age <- lm(part_coef~age_scan+male+fd_mean_avg+avgweight+totalSizet, data=main_unique)
summary(lm_part_coef_age)
lm_clust_co_age <- lm(avgclustco_both~age_scan+male+fd_mean_avg+avgweight+totalSizet, data=main_unique)
summary(lm_clust_co_age)
visreg(lm_between_sys_age,"age_scan", main="Mean between-system connectivity", ylab="", xlab="Age at scan (years)")
visreg(lm_within_sys_age,"age_scan", main= "Mean within-system connectivity", ylab="", xlab="Age at scan (years)")
visreg(lm_segreg_age,"age_scan", main="System segregation", ylab="", xlab="Age at scan (years)")
visreg(lm_modul_avg_age,"age_scan", main="Modularity quality index (Q)", ylab="", xlab="Age at scan (years)")
visreg(lm_part_coef_age,"age_scan", main="Average participation coefficient", ylab="", xlab="Age at scan (years)")
visreg(lm_clust_co_age,"age_scan", main="Average clustering coefficient", ylab="", xlab="Age at scan (years)")

#System-level analyses
networks_age_pvals_fdr #FDR-corrected
lm_sys1to3<- lm(sys1to3~age_scan+male+fd_mean_avg+avgweight+totalSizet, data=main_unique)
lm_sys3to7<- lm(sys3to7~age_scan+male+fd_mean_avg+avgweight+totalSizet, data=main_unique)
lm_sys4to7<- lm(sys4to7~age_scan+male+fd_mean_avg+avgweight+totalSizet, data=main_unique)
summary(lm_sys1to3);visreg(lm_sys1to3,"age_scan")
summary(lm_sys3to7);visreg(lm_sys3to7,"age_scan")
summary(lm_sys4to7);visreg(lm_sys4to7,"age_scan")

#Edge-level analyses
edges_dir=paste0("/cbica/home/tooleyu/projects/in_progress/within_between_network_conn_CBPD/data/imageData/",pipeline,"/schaefer400zNetworks_avg/")
#make a huge dataframe of subjects x edges

files <- data.frame(list.files(edges_dir)) %>% filter( list.files.edges_dir. %in% paste0(main_unique$ID, "_schaefer400MNI_zavgnetwork.txt")); files$list.files.edges_dir. <- paste0(edges_dir,files$list.files.edges_dir.)
# vect <- matrix(nrow = dim(main_unique)[1], ncol = 79801)
# for (i in 1:dim(main_unique)[1]){
# mat <- data.frame(read_csv(file =paste0(edges_dir,main_unique$ID[i], "_schaefer400MNI_zavgnetwork.txt"), col_names = F))
# vect[i,2:79801]<-mat[lower.tri(mat)] #read out the lower triangle from matrix by row
# }
# 
# #code to read back into a 400 x 400 matrix! This works correctly.
# vec<- mat[lower.tri(mat)]
# goback <- matrix(nrow = 400, ncol=400)
# goback[lower.tri(goback, diag=FALSE)] <- vec
# goback <- t(goback)
# goback
# 
# saveRDS(vect,"~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/n92_all_edges.Rdata")

#Read back in the matrix of edges, so you don't create every time.
#n90_all_edges<- readRDS("~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/n90_all_edges.Rds");n90_all_edges <- n90_all_edges[,-1]
n92_all_edges<- readRDS("~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/n92_all_edges.RData");n92_all_edges <- n92_all_edges[,-1]


## Run edgwise linear models
#Run linear models on edges, controlling for age, sex and average motion.
# edgewise_age_pvals<- apply(n90_all_edges,2, function(x) { summary(lm(x~main_unique$age_scan+main_unique$male+main_unique$fd_mean_avg+main_unique$avgweight+main_unique$pctVolsCensored+main_unique$totalSizet))$coef[2,4]}) #without multicore

#multi-core apply the linear model across the matrix of edges, get age pvals 
edgewise_age_pvals<- mclapply(1:79800, function(x) { summary(lm(n92_all_edges[,x]~main_unique$age_scan+main_unique$male+main_unique$fd_mean_avg+main_unique$avgweight+main_unique$totalSizet))$coef[2,4]}, mc.cores = 4)
edgewise_age_pvals <- unlist(edgewise_age_pvals)
edgewise_age_pvals_fdr <- cbind(edgewise_age_pvals,p.adjust(edgewise_age_pvals,method = "fdr"))
#get age betas
edgewise_age_betas<- mclapply(1:79800, function(x) { lm.beta(lm(n92_all_edges[,x]~main_unique$age_scan+main_unique$male+main_unique$fd_mean_avg+main_unique$avgweight+main_unique$totalSizet))$standardized.coefficients[[2]]}, mc.cores = 4)
edgewise_age_betas <- unlist(edgewise_age_betas)
#save for reloading in future
save(edgewise_age_pvals_fdr, edgewise_age_betas,file= "~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/edgewise_age_effects_all_cov_n92.Rdata")

#load back in edge age effects
load("~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/edgewise_age_effects_all_cov_n92.Rdata")

#code to read back into a 400 x 400 matrix! This works correctly.
# vec<- mat[lower.tri(mat)]
edge_age_pvals_mat <- matrix(nrow = 400, ncol=400)
edge_age_pvals_mat[lower.tri(edge_age_pvals_mat, diag=FALSE)] <- edgewise_age_pvals_fdr[,1] #take uncorrected p-values for now
#copy lower triangle to upper triangle
edge_age_pvals_mat[upper.tri(edge_age_pvals_mat)] <- t(edge_age_pvals_mat)[upper.tri(edge_age_pvals_mat)]
edge_age_betas_mat <- matrix(nrow = 400, ncol=400)
edge_age_betas_mat[lower.tri(edge_age_betas_mat, diag=FALSE)] <- edgewise_age_betas
edge_age_betas_mat[upper.tri( edge_age_betas_mat)] <- t(edge_age_betas_mat)[upper.tri( edge_age_betas_mat)]

# Cognition associations
## Play
lm_sys3to7<- lm(sys3to7~age_scan+male+fd_mean_avg+avgweight+totalSizet+litnum_freq_play_avg, data=main_unique)
lm_sys4to7<- lm(sys4to7~age_scan+male+fd_mean_avg+avgweight+totalSizet+litnum_freq_play_avg, data=main_unique)
summary(lm(litnum_freq_play_avg~age_behav.x+male+ses_composite, data=main_unique));visreg(lm(litnum_freq_play_avg~age_behav.x+male+ses_composite, data=main_unique))
summary(lm_sys3to7) #Not significant
summary(lm_sys4to7)
visreg(lm_sys3to7, "sys3to7", main="Frequency of play",ylab="", xlab="DMN to DAN")
visreg(lm_sys4to7, "litnum_freq_play_avg", main="DMN to VAN",ylab="", xlab="litnum_freq_play_avg")

## Reasoning

#Test all 3 systems with Reasoning
main_unique$matrix_reasoning_both <- ifelse(is.na(main_unique$wppsi_matrix_valid),main_unique$wisc_matrix_scaled,ifelse(main_unique$wppsi_matrix_valid==0,NA,main_unique$wppsi_matrix_scaled))
main_unique$matrix_reasoning_both_raw <- ifelse(is.na(main_unique$wppsi_matrix_valid),main_unique$wisc_matrix_raw,ifelse(main_unique$wppsi_matrix_valid==0,NA,main_unique$wppsi_matrix_raw))
main_unique$matrix_reasoning_both;main_unique$matrix_reasoning_both_raw
lm_sys1to3<- lm(matrix_reasoning_both~age_scan+male+fd_mean_avg+avgweight+totalSizet+sys1to3, data=main_unique)
lm_sys3to7<- lm(matrix_reasoning_both~age_scan+male+fd_mean_avg+avgweight+totalSizet+sys3to7, data=main_unique)
lm_sys4to7<- lm(matrix_reasoning_both~age_scan+male+fd_mean_avg+avgweight+totalSizet+sys4to7, data=main_unique)
summary(lm_sys1to3)
summary(lm_sys3to7)
summary(lm_sys4to7)
cohens_f(lm_sys1to3, ci = 0.95)
cohens_f(lm_sys3to7, ci = 0.95)
cohens_f(lm_sys4to7, ci = 0.95)
visreg(lm_sys1to3, "sys1to3", main= "Matrix Reasoning (scaled score)", ylab="Scaled score", xlab="VIS to DA connectivity")
visreg(lm_sys3to7, "sys3to7", main= "Matrix Reasoning (scaled score)", ylab="Scaled score", xlab="DM to DA connectivity")
visreg(lm_sys4to7)


# Replication including longitudinal data ---------------------------------
pipeline='gsr_spkreg_fd0.5dvars1.75_drpvls'
network_dir=paste0("~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/",pipeline)

load(paste0(network_dir,"/CBPD_n92_schaefer400_allruns.Rdata"))
#melt data to get longitudinal data only
main_filt_long <- dplyr::select(main_filt, c(ID:record_id.y,date_scan, date_behav, age_scan,male,fd_mean_avg,avgweight,pctVolsCensored,totalSizet,race2,ethnicity, longitudinal_visit_num,base_ID, part_coef))
#need to filter out only the first datapoint from each timepoint, since the network data is repeated across multiple runs at each timepoint
#to see why -> main_filt %>% select(base_ID,longitudinal_visit_num,sys1to1)
main_long_only <- main_filt %>% group_by(base_ID, longitudinal_visit_num) %>% filter(row_number() == 1)
main_long_only %>% select(base_ID,longitudinal_visit_num,sys1to1) #now has one set of network stats from each timepoint for each subject
table(main_long_only$longitudinal_visit_num)
table(main_unique$longitudinal_visit_num) #but we had a few longitudinal visits in the original dataset

#get gaps between visits
main_long_only<- main_long_only %>% group_by(base_ID) %>% mutate(diffageT1toT2=nth(age_scan,2)-first(age_scan),
                                                                 diffageT2toT3=nth(age_scan,3)-nth(age_scan,2)) %>% ungroup()
#Summary in months
main_long_only %>% select(diffageT1toT2) %>%  summarise(.,min=min(diffageT1toT2*12, na.rm=T), max=max(diffageT1toT2*12, na.rm=T), mean=mean(diffageT1toT2*12, na.rm=T),med=median(diffageT1toT2*12, na.rm=T), sd=sd(diffageT1toT2*12, na.rm=T))
main_long_only %>% select(diffageT2toT3) %>%  summarise(.,min=min(diffageT2toT3*12, na.rm=T), max=max(diffageT2toT3*12, na.rm=T),mean=mean(diffageT2toT3*12, na.rm=T),med=median(diffageT2toT3*12, na.rm=T), sd=sd(diffageT2toT3*12, na.rm=T))

#Global measures analyses
lm_within_sys_age <- lm(mean_within_sys~age_scan+male+fd_mean_avg+avgweight+totalSizet, data=main_long_only)
summary(lm_within_sys_age)
lm_between_sys_age <- lm(mean_between_sys~age_scan+male+fd_mean_avg+avgweight+totalSizet, data=main_long_only)
summary(lm_between_sys_age)
lm_segreg_age<- lm(system_segreg~age_scan+male+fd_mean_avg+avgweight+totalSizet, data=main_long_only)
summary(lm_segreg_age)
lm_modul_avg_age <- lm(modul_avg~age_scan+male+fd_mean_avg+avgweight+totalSizet, data=main_long_only)
summary(lm_modul_avg_age)
lm_part_coef_age <- lm(part_coef~age_scan+male+fd_mean_avg+avgweight+totalSizet, data=main_long_only)
summary(lm_part_coef_age)
lm_clust_co_age <- lm(avgclustco_both~age_scan+male+fd_mean_avg+avgweight+totalSizet, data=main_long_only)
summary(lm_clust_co_age)
visreg(lm_between_sys_age,"age_scan", main="Mean between-system connectivity", ylab="", xlab="Age at scan (years)")
visreg(lm_within_sys_age,"age_scan", main= "Mean within-system connectivity", ylab="", xlab="Age at scan (years)")
visreg(lm_segreg_age,"age_scan", main="System segregation", ylab="", xlab="Age at scan (years)")
visreg(lm_modul_avg_age,"age_scan", main="Modularity quality index (Q)", ylab="", xlab="Age at scan (years)")
visreg(lm_part_coef_age,"age_scan", main="Average participation coefficient", ylab="", xlab="Age at scan (years)")
visreg(lm_clust_co_age,"age_scan", main="Average clustering coefficient", ylab="", xlab="Age at scan (years)")

#Segregation and age
outcomes <- main_long_only %>% ungroup() %>% dplyr::select(mean_within_sys,mean_between_sys,system_segreg, modul_avg, part_coef,avgclustco_both) %>% names()
for (outcome in outcomes){
  model<-paste0("lm_",outcome,"_age")
  formula<-formula(paste0(outcome, '~age_scan+male+fd_mean_avg+avgweight+totalSizet'))
  assign(model, lm(formula, data=main_long_only))
  output <- paste0("summary_",model)
  assign(output, apa_print(Anova(get(model)), es="pes",mse=FALSE))
  eta <- paste0("eta_",model)
  assign(eta, eta_squared(Anova(get(model)), ci=0.95))
  cohensf <- paste0("cohensf_",model)
  assign(cohensf, cohens_f(Anova(get(model)), ci=0.95))
  beta <- paste0("beta_",model)
  assign(beta, lm.beta(get(model)))
  print(outcome)
  print(get(beta))
  print(summary(get(model)))
}

## FDR-correct across systems
# Look at each column and correct for multiple comparisons
covariates="~ age_scan+male+fd_mean_avg+avgweight+totalSizet"
#make a dataframe with no repeats of net comparisons
main_long_only <- dplyr::select(main_long_only, -c(sys2to1,sys3to1,sys3to2,sys4to1,sys4to2,sys4to3,sys5to1,sys5to2,sys5to3,sys5to4,sys6to1,sys6to2,sys6to3,sys6to4,sys6to5,sys7to1,sys7to2,sys7to3,sys7to4,sys7to5,sys7to6))
#run and compare for multiple comparisons again
m <- mclapply(names((dplyr::select(ungroup(main_long_only),sys1to1:sys7to7))), function(sys) {as.formula(paste(sys, covariates, sep=""))},mc.cores=2)
networks_Age_pvals <- mclapply(m, function(sys) { summary(lm(formula = sys,data=main_long_only))$coefficients[2,4]},mc.cores=1)
networks_Age_pvals <- as.data.frame(networks_Age_pvals)
networks_Age_pvals <- t(networks_Age_pvals)
networks_Age_pvals <- as.data.frame(networks_Age_pvals)
#bonferroni correct
networks_Age_pvals_fdr <- p.adjust(networks_Age_pvals$V1, method="fdr")
networks_age_pvals_fdr <- data.frame(networks_Age_pvals_fdr,names((dplyr::select(ungroup(main_long_only),sys1to1:sys7to7))))
colnames(networks_age_pvals_fdr) <- c("pvalue", "network")
#FDR correction shows DMN to attentional networks and visual to dorsal attention is marginal.
print(networks_age_pvals_fdr)
networks_age_pvals_fdr
lm_sys1to3<- lm(sys1to3~age_scan+male+fd_mean_avg+avgweight+totalSizet, data=main_long_only)
lm_sys3to7<- lm(sys3to7~age_scan+male+fd_mean_avg+avgweight+totalSizet, data=main_long_only)
lm_sys4to7<- lm(sys4to7~age_scan+male+fd_mean_avg+avgweight+totalSizet, data=main_long_only)
lm_sys3to4<- lm(sys3to4~age_scan+male+fd_mean_avg+avgweight+totalSizet, data=main_long_only)
lm_sys1to1<- lm(sys1to1~age_scan+male+fd_mean_avg+avgweight+totalSizet, data=main_long_only)
lm_sys1to7<- lm(sys1to7~age_scan+male+fd_mean_avg+avgweight+totalSizet, data=main_long_only)
lm_sys3to6<- lm(sys3to6~age_scan+male+fd_mean_avg+avgweight+totalSizet, data=main_long_only)
lm_sys7to7<- lm(sys7to7~age_scan+male+fd_mean_avg+avgweight+totalSizet, data=main_long_only)
lm_sys4to4<- lm(sys4to4~age_scan+male+fd_mean_avg+avgweight+totalSizet, data=main_long_only)
summary(lm_sys1to3);visreg(lm_sys1to3,"age_scan");lm.beta(lm_sys1to3)
summary(lm_sys3to7);visreg(lm_sys3to7,"age_scan");lm.beta(lm_sys3to7)
summary(lm_sys4to7);visreg(lm_sys4to7,"age_scan");lm.beta(lm_sys4to7)
summary(lm_sys3to4);visreg(lm_sys3to4,"age_scan");lm.beta(lm_sys3to4)
summary(lm_sys1to1);visreg(lm_sys1to1,"age_scan");lm.beta(lm_sys1to1)
summary(lm_sys1to7);visreg(lm_sys1to7,"age_scan");lm.beta(lm_sys1to7)
summary(lm_sys3to6);visreg(lm_sys3to6,"age_scan");lm.beta(lm_sys3to6)
summary(lm_sys7to7);visreg(lm_sys7to7,"age_scan");lm.beta(lm_sys7to7)
summary(lm_sys4to4);visreg(lm_sys4to4,"age_scan");lm.beta(lm_sys7to7)

## Edge-level analyses
pipeline='gsr_spkreg_fd0.5dvars1.75_drpvls'
edges_dir=paste0("/cbica/home/tooleyu/projects/in_progress/within_between_network_conn_CBPD/data/imageData/",pipeline,"/schaefer400zNetworks_avg/")
#make a huge dataframe of subjects x edges

# files <- data.frame(list.files(edges_dir)) %>% filter( list.files.edges_dir. %in% paste0(main_long_only$ID, "_schaefer400MNI_zavgnetwork.txt")); files$list.files.edges_dir. <- paste0(edges_dir,files$list.files.edges_dir.)
# vect <- matrix(nrow = dim(main_long_only)[1], ncol = 79801)
# for (i in 1:dim(main_long_only)[1]){
# mat <- data.frame(read_csv(file =paste0(edges_dir,main_long_only$ID[i], "_schaefer400MNI_zavgnetwork.txt"), col_names = F))
# vect[i,2:79801]<-mat[lower.tri(mat)] #read out the lower triangle from matrix by row
# }
# 
# #code to read back into a 400 x 400 matrix! This works correctly.
# vec<- mat[lower.tri(mat)]
# goback <- matrix(nrow = 400, ncol=400)
# goback[lower.tri(goback, diag=FALSE)] <- vec
# goback <- t(goback)
# goback
# saveRDS(vect,"~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/n118_all_edges.Rdata")

#Read back in the matrix of edges, so you don't create every time.
n118_all_edges<- readRDS("~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/n118_all_edges.RData");
#multi-core apply the linear model across the matrix of edges, get age pvals 
# edgewise_age_pvals<- mclapply(1:79800, function(x) { summary(lm(n118_all_edges[,x]~main_long_only$age_scan+main_long_only$male+main_long_only$fd_mean_avg+main_long_only$avgweight+main_long_only$totalSizet))$coef[2,4]}, mc.cores = 4)
# edgewise_age_pvals <- unlist(edgewise_age_pvals)
# edgewise_age_pvals_fdr <- cbind(edgewise_age_pvals,p.adjust(edgewise_age_pvals,method = "fdr"))
# #get age betas
# edgewise_age_betas<- mclapply(1:79800, function(x) { lm.beta(lm(n118_all_edges[,x]~main_long_only$age_scan+main_long_only$male+main_long_only$fd_mean_avg+main_long_only$avgweight+main_long_only$totalSizet))$standardized.coefficients[[2]]}, mc.cores = 4)
# edgewise_age_betas <- unlist(edgewise_age_betas)
# #save for reloading in future
# save(edgewise_age_pvals_fdr, edgewise_age_betas,file= "~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/replicate_longitudinal/edgewise_age_effects_all_cov_n118.Rdata")

#load back in edge age effects
load("~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/replicate_longitudinal/edgewise_age_effects_all_cov_n118.Rdata")

#code to read back into a matrix! This works correctly.
edge_age_pvals_mat <- matrix(nrow = 400, ncol=400)
edge_age_pvals_mat[lower.tri(edge_age_pvals_mat, diag=FALSE)] <- edgewise_age_pvals_fdr[,1] #take uncorrected p-values for now
#copy lower triangle to upper triangle
edge_age_pvals_mat[upper.tri(edge_age_pvals_mat)] <- t(edge_age_pvals_mat)[upper.tri(edge_age_pvals_mat)]
edge_age_betas_mat <- matrix(nrow = 400, ncol=400)
edge_age_betas_mat[lower.tri(edge_age_betas_mat, diag=FALSE)] <- edgewise_age_betas
edge_age_betas_mat[upper.tri( edge_age_betas_mat)] <- t(edge_age_betas_mat)[upper.tri( edge_age_betas_mat)]

## Plot on brain
schaefer_atlas_region_names_lh = get.atlas.region.names('Schaefer2018_400Parcels_7Networks_order', template_subjects_dir = subjects_dir,template_subject=subject_id, hemi="lh");
schaefer_atlas_region_names_rh = get.atlas.region.names('Schaefer2018_400Parcels_7Networks_order', template_subjects_dir = subjects_dir,template_subject=subject_id, hemi="rh");
output_image_directory="/Users/utooley/Dropbox/projects/in_progress/within_between_network_conn_CBPD/output/figures/replicate_longitudinal/"
pvalues=c(0.01,0.001,0.0001,0.00001)
#Use different pval thresholds and plot parcels that have age-nonlin-sig edges at that threshold
pos_edges <- edge_age_pvals_mat;pos_edges[edge_age_betas_mat<0] <- NA
neg_edges <- edge_age_pvals_mat;neg_edges[edge_age_betas_mat>0] <- NA
for (pvalue in pvalues){
  indices <- which(pos_edges<pvalue, arr.ind = T) #indices of edges that have age effects  below a pval
  l <- data.frame(table(indices)) #for each parcel, how many age-significant edges does it have at a given pval thresh
  values <- vector(mode = "double", 400)
  for (i in 1:400){
    values[as.numeric(as.character(l$indices[i]))] <- l$Freq[i]#replace the indices in values with the num of sig edges from that parcel in the table
  }
  values <- values/2 #because we indexed the full matrix, there are duplicates for each edge
  print(max(values))
  num_of_sig_age_edges_lh=as.list(setNames(c(0, values[1:200]), schaefer_atlas_region_names_lh))
  num_of_sig_age_edges_rh=as.list(setNames(c(0, values[201:400]), schaefer_atlas_region_names_rh))
  #colormap
  colFn_diverging = colorRampPalette(c("white","#E0011C"));makecmap_options=list('colFn'=colFn_diverging, range= c(0,20)) 
  rglactions=list("snapshot_png"=paste0(output_image_directory,"age_pos_pvals_", pvalue,"regions.png")) #purple="#7502E3", red=#E0011C, blue=#3602D9
  vis.region.values.on.subject(subjects_dir, 'fsaverage6', 'Schaefer2018_400Parcels_7Networks_order',  num_of_sig_age_edges_lh, 
                               num_of_sig_age_edges_rh, makecmap_options = makecmap_options, "inflated", views="t4", draw_colorbar = T, rgloptions = rgloptions, rglactions = rglactions)
}
n118_all_edges<- readRDS("~/Documents/projects/in_progress/within_between_network_conn_CBPD/data/imageData/gsr_spkreg_fd0.5dvars1.75_drpvls/n118_all_edges.RData");

#get nodal stregth for each node for each subject
n118_nodal_strength <- apply(n118_all_edges,1, function (x){
  goback <- matrix(nrow = 400, ncol=400)
  goback[lower.tri(goback, diag=FALSE)] <-x
  goback[upper.tri(goback)] <- t(goback)[upper.tri(goback)]
  return(colMeans(goback, na.rm = T))
})
n118_nodal_strength <- t(n118_nodal_strength)
n118_nodal_strength_data <- data.frame(main_long_only,n118_nodal_strength)
#parcelwise
parcel_sd_pvals<- lapply(names(n118_nodal_strength_data[,1669:2068]), function(x) { summary(lm(as.formula(paste0(x,"~age_scan+male+fd_mean_avg+avgweight+totalSizet")), data=n118_nodal_strength_data))$coef[2,4]})
parcel_sd_pvals <- unlist(parcel_sd_pvals)
parcel_sd_pvals_fdr <- cbind(parcel_sd_pvals,p.adjust(parcel_sd_pvals,method = "fdr"))
#get age betas
parcel_sd_betas<- lapply(names(n118_nodal_strength_data[,1669:2068]), function(x) { lm.beta(lm(as.formula(paste0(x,"~age_scan+male+fd_mean_avg+avgweight+totalSizet")), data=n118_nodal_strength_data))$standardized.coefficients[[2]]})
parcel_sd_betas <- unlist(parcel_sd_betas)
name="parcel_age_node_strength_pvals"
to_plot=parcel_sd_pvals_fdr[,1] #Change this to plot different things
to_plot=ifelse(parcel_sd_pvals_fdr[,2] <0.05,parcel_sd_betas,NA)

lh=as.list(setNames(c(NA,to_plot[1:200]), schaefer_atlas_region_names_lh)) #this has the medial wall in it.
rh=as.list(setNames(c(NA,to_plot[201:400]), schaefer_atlas_region_names_rh))
#colormap
#colFn_diverging = colorRampPalette(rev(c("white","white","pink","red")));makecmap_options=list('colFn'=colFn_diverging)
colFn_diverging = colorRampPalette(rev(RColorBrewer::brewer.pal(9, name="RdBu")));makecmap_options=list('colFn'=colFn_diverging, symm=T)
rglactions=list("snapshot_png"=paste0(output_image_directory,name,".png"), 'shift_hemis_apart'=TRUE)
vis.region.values.on.subject(subjects_dir, 'fsaverage6', 'Schaefer2018_400Parcels_7Networks_order',  lh, 
                             rh, makecmap_options = makecmap_options, "inflated", views="t4", draw_colorbar = T, rgloptions = rgloptions, rglactions = rglactions)
## Reasoning
#Test all 3 systems with Reasoning
main_long_only$matrix_reasoning_both <- ifelse(is.na(main_long_only$wppsi_matrix_valid),main_long_only$wisc_matrix_scaled,ifelse(main_long_only$wppsi_matrix_valid==0,NA,main_long_only$wppsi_matrix_scaled))
main_long_only$matrix_reasoning_both_raw <- ifelse(is.na(main_long_only$wppsi_matrix_valid),main_long_only$wisc_matrix_raw,ifelse(main_long_only$wppsi_matrix_valid==0,NA,main_long_only$wppsi_matrix_raw))
main_long_only$matrix_reasoning_both;main_long_only$matrix_reasoning_both_raw
lm_sys1to3<- lm(matrix_reasoning_both~age_scan+male+fd_mean_avg+avgweight+totalSizet+sys1to3, data=main_long_only)
lm_sys3to7<- lm(matrix_reasoning_both~age_scan+male+fd_mean_avg+avgweight+totalSizet+sys3to7, data=main_long_only)
lm_sys4to7<- lm(matrix_reasoning_both~age_scan+male+fd_mean_avg+avgweight+totalSizet+sys4to7, data=main_long_only)
summary(lm_sys1to3);lm.beta(lm_sys1to3)
summary(lm_sys3to7);lm.beta(lm_sys3to7)
summary(lm_sys4to7);lm.beta(lm_sys4to7)
visreg(lm_sys1to3, "sys1to3", main= "Matrix Reasoning (scaled score)", ylab="Scaled score", xlab="VIS to DA connectivity")
visreg(lm_sys3to7, "sys3to7", main= "Matrix Reasoning (scaled score)", ylab="Scaled score", xlab="DM to DA connectivity")
visreg(lm_sys4to7)

## Play
lm_sys3to7<- lm(litnum_freq_play_avg~age_scan+male+fd_mean_avg+avgweight+totalSizet+sys3to7+ses_composite, data=main_long_only)
lm_sys4to7<- lm(litnum_freq_play_avg~age_scan+male+fd_mean_avg+avgweight+totalSizet+sys4to7, data=main_long_only)
summary(lm(litnum_freq_play_avg~age_behav.x+male+ses_composite, data=main_replicate_unique))
#visreg(lm(litnum_freq_play_avg~age_behav.x+male+ses_composite, data=main_replicate_unique))
summary(lm_sys3to7);lm.beta(lm_sys3to7) #Not significant
summary(lm_sys4to7)
visreg(lm_sys3to7, "sys3to7", main="Frequency of play",ylab="", xlab="DMN to DAN")
visreg(lm_sys4to7, "litnum_freq_play_avg", main="DMN to VAN",ylab="", xlab="litnum_freq_play_avg")

# Power analyses ----------------------------------------------------------
library(pwr)
library(sensemakr)
http://web.pdx.edu/~newsomj/mvclass/ho_sample%20size.pdf
#convert f2 to betas

pwr.f2.test(u = 6,v = 86,sig.level = 0.05,f2 = 0.16)

pwr.r.test(n=92,sig.level = 0.05, r = 0.3)

library(effectsize)
cohens_f(lm_sys1to3)
