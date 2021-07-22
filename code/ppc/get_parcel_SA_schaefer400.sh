export SUBJECTS_DIR=/cbica/projects/cbpd_main_data/CBPD_bids/derivatives/freesurfer_edits_t3

annot_file_lh=/cbica/projects/cbpd_main_data/tools/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/fsaverage/label/lh.Schaefer2018_400Parcels_7Networks_order.annot
annot_file_rh=/cbica/projects/cbpd_main_data/tools/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/fsaverage/label/rh.Schaefer2018_400Parcels_7Networks_order.annot

for subjid in `find ${SUBJECTS_DIR} -iname "sub-CBPD*" -exec basename {} \;`;
do

#if [[ -f ${SUBJECTS_DIR}/${subjid}/label/lh.Schaefer2018_400Parcels_7Networks_order.annot ]]; then
#  echo it already exists
  #extract surface area
  mri_surf2surf --srcsubject fsaverage --trgsubject $subjid --hemi lh  --sval-annot ${annot_file_lh}  --tval ${SUBJECTS_DIR}/${subjid}/label/lh.Schaefer2018_400Parcels_7Networks_order.annot
  mri_surf2surf --srcsubject fsaverage --trgsubject $subjid --hemi rh  --sval-annot ${annot_file_rh}  --tval ${SUBJECTS_DIR}/${subjid}/label/rh.Schaefer2018_400Parcels_7Networks_order.annot

  declare -a array=(lh
  rh)
  for hemi in "${array[@]}";
  do
  annot_file=${SUBJECTS_DIR}/${subjid}/label/${hemi}.Schaefer2018_400Parcels_7Networks_order.annot

  #get thickness and surface area
  mri_segstats --annot ${subjid} ${hemi} ${annot_file} --i $SUBJECTS_DIR/${subjid}/surf/${hemi}.thickness --snr --sum ${SUBJECTS_DIR}/${subjid}/stats/${hemi}.schaefer400_7nets.thickness.stats
  mri_segstats --annot ${subjid} ${hemi} ${annot_file} --i $SUBJECTS_DIR/${subjid}/surf/${hemi}.area --accumulate --snr --sum ${SUBJECTS_DIR}/${subjid}/stats/${hemi}.schaefer400_7nets.surfarea.stats

  #when doing area or volume, use --accumulate flag to mri_segstats to get the total area or volume. Don't do with CT to get mean.
  done
#else
  # echo it doesnt exist
  # echo ${subjid}
#fi
done

#Must run this section with Python 2!
conda activate python2

subjlist=`find ${SUBJECTS_DIR} -iname "sub-CBPD*" -exec basename {} \;`
declare -a array=(lh
rh)
for hemi in "${array[@]}";
do
aparcstats2table --subjects ${subjlist} --hemi ${hemi} --parc schaefer400_7nets.surfarea --meas area --tablefile ${SUBJECTS_DIR}/schaefer400_${hemi}_surfarea_stats.txt --skip
asegstats2table --subjects ${subjlist} --statsfile ${hemi}.schaefer400_7nets.thickness.stats --meas mean --tablefile ${SUBJECTS_DIR}/schaefer400_${hemi}_thickness_stats.txt --skip
#for some reason aparcstats2table didn't work for thickness, so use asegstats2table
done
