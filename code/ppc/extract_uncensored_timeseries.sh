SCRIPTS_DIR=/cbica/home/tooleyu/projects/in_progress/within_between_network_conn_CBPD/code/ppc
pipeline=xcpEngine_gsr_spkreg_fd0.5dvars1.75_drpvls
subject_list=/cbica/home/tooleyu/projects/in_progress/within_between_network_conn_CBPD/data/subjectLists/n150_cross_sect_one_or_more_nonsleep_rest_at_least_130_vols.csv
SIMG=/cbica/projects/cbpd_main_data/tools/singularity/xcpEngine-100219.simg
data_dir=/cbica/projects/cbpd_main_data/CBPD_bids_crosssectional/derivatives/xcpEngine_gsr_spkreg_fd0.5dvars1.75_drpvls

for line in `cat ${subject_list}`
do
  sub=`echo $line | cut -d, -f1`
  echo ${sub}
  run=`echo $line | cut -d, -f2| tr -d '\r\n'` #or just \n if it isn't working
  echo ${run}
  newsub=${sub}
  if [[ -e ${data_dir}/${sub}/${run}/fcon/schaefer400 ]]; then
    echo 'we have fcon/schaefer400 for' ${sub} ${run}
    singularity exec --cleanenv $SIMG /xcpEngine/utils/roi2ts.R -i ${data_dir}/${sub}/${run}/regress/${sub}_${run}_uncensored.nii.gz -r ${data_dir}/${sub}/${run}/${sub}_${run}_atlas/${sub}_${run}_schaefer400.nii.gz -l  /xcpEngine/atlas/schaefer400x7/schaefer400x7NodeIndex.1D > ${data_dir}/${sub}/${run}/fcon/schaefer400/${sub}_${run}_schaefer400_uncensored_ts.1D
    #if there is no uncensored file found, then it means there were 0 volumes censored.
  elif [[ -e ${data_dir}/${sub}/${run}/fcon/schaefer400x7 ]]; then
    echo 'we have fcon/schaefer400x7 for' ${sub} ${run}
    singularity exec --cleanenv $SIMG /xcpEngine/utils/roi2ts.R -i ${data_dir}/${sub}/${run}/regress/${sub}_${run}_uncensored.nii.gz -r ${data_dir}/${sub}/${run}/${sub}_${run}_atlas/${sub}_${run}_schaefer400x7.nii.gz -l /xcpEngine/atlas/schaefer400x7/schaefer400x7NodeIndex.1D > ${data_dir}/${sub}/${run}/fcon/schaefer400x7/${sub}_${run}_schaefer400x7_uncensored_ts.1D
  fi
done


#for uncensored files that were not found, check if that's why and if so copy over the censored timeseries

for line in `cat ${subject_list}`
do
  sub=`echo $line | cut -d, -f1`
  echo ${sub}
  run=`echo $line | cut -d, -f2| tr -d '\r\n'` #or just \n if it isn't working
  echo ${run}
  newsub=${sub}
  if ! [[ -e ${data_dir}/${sub}/${run}/regress/${sub}_${run}_uncensored.nii.gz ]]; then
    echo 'there is no uncensored image for '${sub} ${run}
    if grep -q "1" ${data_dir}/${sub}/${run}/${sub}_${run}-nFlags.1D; then
      echo 'but volumes have been censored....'
    else
      echo 'and no volumes were censored'
      #copy over the timeseries and rename
      if [[ -e ${data_dir}/${sub}/${run}/fcon/schaefer400 ]]; then
        echo cp ${data_dir}/${sub}/${run}/fcon/schaefer400/${sub}_${run}_schaefer400_ts.1D ${data_dir}/${sub}/${run}/fcon/schaefer400/${sub}_${run}_schaefer400_uncensored_ts.1D
      elif [[ -e ${data_dir}/${sub}/${run}/fcon/schaefer400x7 ]]; then
        echo cp ${data_dir}/${sub}/${run}/fcon/schaefer400x7/${sub}_${run}_schaefer400x7_ts.1D  ${data_dir}/${sub}/${run}/fcon/schaefer400x7/${sub}_${run}_schaefer400x7_uncensored_ts.1D
      fi
    fi
fi
done
