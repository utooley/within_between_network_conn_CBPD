export SUBJECTS_DIR=/cbica/projects/cbpd_main_data/CBPD_bids/derivatives/freesurfer_edits_t1

module load freesurfer/6.0.0.patched

cd $SUBJECTS_DIR/qdec
mri_glmfit-sim \
  --glmdir rh_area_agesext1rating_fwhm10 \
  --perm 1000 2.3 abs \
  --perm-resid \
  --cwp  0.05\
  --2spaces


# mri_glmfit.bin --y /cbica/projects/cbpd_main_data/CBPD_bids/derivatives/freesurfer_edits_t1/qdec/rh_areaagesext1rating_fwhm10/y.mgh --fsgd /cbica/projects/cbpd_main_data/CBPD_bids/derivatives/freesurfer_edits_t1/qdec/rh_areaagesext1rating_fwhm10/qdec.fsgd dods --glmdir rh_areaagesext1rating_fwhm10 --surf fsaverage lh --label /cbica/projects/cbpd_main_data/CBPD_bids/derivatives/freesurfer_edits_t1/fsaverage/label/lh.aparc.label --C /cbica/projects/cbpd_main_data/CBPD_bids/derivatives/freesurfer_edits_t1/qdec/rh_areaagesext1rating_fwhm10/contrasts/lh-Avg-Intercept-area.mat --C /cbica/projects/cbpd_main_data/CBPD_bids/derivatives/freesurfer_edits_t1/qdec/rh_areaagesext1rating_fwhm10/contrasts/lh-Avg-area-age_scan-Cor.mat --eres-save
#
#
# mri_glmfit.bin --y /gpfs/fs001/cbica/projects/cbpd_main_data/CBPD_bids/derivatives/freesurfer_edits_t1/qdec/rh_area_agesext1rating_fwhm10/y.mgh --fsgd /gpfs/fs001/cbica/projects/cbpd_main_data/CBPD_bids/derivatives/freesurfer_edits_t1/qdec/rh_area_agesext1rating_fwhm10/qdec.fsgd dods --glmdir /gpfs/fs001/cbica/projects/cbpd_main_data/CBPD_bids/derivatives/freesurfer_edits_t1/qdec/rh_area_agesext1rating_fwhm10 --surf fsaverage rh --label /gpfs/fs001/cbica/projects/cbpd_main_data/CBPD_bids/derivatives/freesurfer_edits_t1/fsaverage/label/rh.aparc.label --C /gpfs/fs001/cbica/projects/cbpd_main_data/CBPD_bids/derivatives/freesurfer_edits_t1/qdec/rh_area_agesext1rating_fwhm10/contrasts/rh-Avg-Intercept-area.mat --C /gpfs/fs001/cbica/projects/cbpd_main_data/CBPD_bids/derivatives/freesurfer_edits_t1/qdec/rh_area_agesext1rating_fwhm10/contrasts/rh-Avg-area-age_scan-Cor.mat --eres-save
