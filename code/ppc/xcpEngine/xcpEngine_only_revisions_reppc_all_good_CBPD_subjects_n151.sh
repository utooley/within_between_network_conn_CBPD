#Do it by run!
SCRIPTS_DIR=/cbica/home/tooleyu//projects/in_progress/within_between_network_conn_CBPD/code/ppc
pipeline=xcpEngine_gsr_spkreg_fd0.25dvars1.75_drpvls
my_BIDS_dir=/cbica/projects/cbpd_main_data/CBPD_bids_crosssectional
subject_list=/cbica/home/tooleyu//projects/in_progress/within_between_network_conn_CBPD/data/subjectLists/n125_cross_sect_one_or_more_nonsleep_rest_10mm_max_RMS_perrun_at_least_130_vols.csv
subject_list_dir=/cbica/home/tooleyu//projects/in_progress/within_between_network_conn_CBPD/data/subjectLists
for line in `cat ${subject_list}`
do
  sub=`echo $line | cut -d, -f1`
  echo ${sub}
  run=`echo $line | cut -d, -f2| tr -d '\r\n'` #or just \n if it isn't working
  echo ${run}
  newsub=${sub}
  if [[ -f /cbica/projects/cbpd_main_data/CBPD_bids_crosssectional/derivatives/fmriprep/${sub}/func/${sub}_task-rest_${run}_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz ]]; then
    echo 'fmriprep exists for' ${sub} ${run}
    if [[ -e /cbica/projects/cbpd_main_data/CBPD_bids_crosssectional/derivatives/${pipeline}/${sub}/${run}/fcon/schaefer400x7/${sub}_${run}_schaefer400x7_network.txt ]]; then
      echo 'Xcpengine run already for run' ${run} 'sub' ${sub}
    elif [[ -e /cbica/projects/cbpd_main_data/CBPD_bids_crosssectional/derivatives/${pipeline}/${sub}/${run}/fcon/schaefer400/${sub}_${run}_schaefer400_network.txt ]]; then
      echo 'Xcpengine run already for run' ${run} 'sub' ${sub}
    else
      #echo rm /cbica/projects/cbpd_main_data/CBPD_bids_crosssectional/derivatives/${pipeline}/${sub}/${run}/ -R
      echo 'Put run' ${run} 'in list for xcpEngine for' ${sub}
      find ${my_BIDS_dir} -type f | grep "${sub}_task-rest_${run}_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz\$"| while read fname; do
      tmp=$(echo "$fname" | awk -F '_' '{print $7}' ) #this parses on underscores and pulls 'run-01'
      echo $tmp
      fname_mnt=$(echo "$fname" | sed -e 's|/data/|/mnt/|' )
      echo $fname_mnt
      echo ${sub},${tmp},${fname_mnt} >> ${subject_list_dir}/n125_fd0.25_pipeline_reprocess_revisions_for_xcpEngine.csv
      done;
    fi
  else
    echo 'No fMRIprep'
    #qsub -j y ${SCRIPTS_DIR}/fmriprep/fmriprep_cmd.sh sub-${newsub}
  fi
done


bash xcpEngine_batch_job_revisions

#Cleanup and rerun those that didn't finish
SCRIPTS_DIR=/cbica/home/tooleyu//projects/in_progress/within_between_network_conn_CBPD/code/ppc
pipeline=xcpEngine_gsr_spkreg_fd0.25dvars1.75_drpvls
my_BIDS_dir=/cbica/projects/cbpd_main_data/CBPD_bids_crosssectional
subject_list=/cbica/home/tooleyu//projects/in_progress/within_between_network_conn_CBPD/data/subjectLists/n125_cross_sect_one_or_more_nonsleep_rest_10mm_max_RMS_perrun_at_least_130_vols.csv
subject_list_dir=/cbica/home/tooleyu//projects/in_progress/within_between_network_conn_CBPD/data/subjectLists
for line in `cat ${subject_list}`
do
  sub=`echo $line | cut -d, -f1`
  echo ${sub}
  run=`echo $line | cut -d, -f2| tr -d '\r\n'` #or just \n if it isn't working
  echo ${run}
  newsub=${sub}
  if [[ -e /cbica/projects/cbpd_main_data/CBPD_bids_crosssectional/derivatives/${pipeline}/${sub}/${run}/norm/${sub}_${run}_normDice.txt ]]; then
    echo 'Xcpengine run already for run' ${run} 'sub' ${sub}
  else
    echo Xcpengine did not finish for ${run} 'sub' ${sub}
    rm /cbica/projects/cbpd_main_data/CBPD_bids_crosssectional/derivatives/${pipeline}/${sub}/${run}/ -Rf
    echo 'Put run' ${run} 'in list for xcpEngine for' ${sub}
    find ${my_BIDS_dir} -type f | grep "${sub}_task-rest_${run}_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz\$"| while read fname; do
    tmp=$(echo "$fname" | awk -F '_' '{print $7}' ) #this parses on underscores and pulls 'run-01'
    echo $tmp
    fname_mnt=$(echo "$fname" | sed -e 's|/data/|/mnt/|' )
    echo $fname_mnt
    echo ${sub},${tmp},${fname_mnt} >> ${subject_list_dir}/nxx_fd0.25_pipeline_reprocess_revisions_for_xcpEngine_cleanup3.csv
    done;
  fi
done
