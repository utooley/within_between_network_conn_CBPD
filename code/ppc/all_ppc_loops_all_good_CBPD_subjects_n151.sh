declare -a array=(sub-CBPD00022
sub-CBPD00183
sub-CBPD00233
sub-CBPD00483
sub-CBPD00743
sub-CBPD00773
sub-CBPD00783
sub-CBPD00952
sub-CBPD00962
sub-CBPD01252
sub-CBPD01262
sub-CBPD01602
sub-CBPD01622
sub-CBPD0175
sub-CBPD0178
sub-CBPD0179
sub-CBPD0182
sub-CBPD0183
sub-CBPD0185
sub-CBPD0186
sub-CBPD0187
sub-CBPD0188
sub-CBPD0191
sub-CBPD0192
sub-CBPD0193
sub-CBPD0194
sub-CBPD0195
sub-CBPD0196
sub-CBPD0197
sub-CBPD0198
sub-CBPD0199
sub-CBPD0201
sub-CBPD0202
sub-CBPD0204
sub-CBPD0208
sub-CBPD0209
sub-CBPD0210
sub-CBPD0211
sub-CBPD0212)

#remove old data from people who had sleep or runs that can now be included.
for dir in `find . -maxdepth 1 -type d -iname 'xcpEngine_*' -printf "%f\n"`;
do
echo $dir
cd $dir
for sub in "${array[@]}"
do
rm ${sub} -R
done
cd ..
done

#copy over new data from the clean CBPD bds directory.

subject_list=/data/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/subjectLists/n151_all_cbpd_one_or_more_nonsleep_rest_at_least_130_vols.txt
julia_BIDS_dir=/data/picsl/mackey_group/BPD/CBPD_bids
my_BIDS_dir=/data/picsl/mackey_group/CBPD/CBPD_bids

for sub in `cat ${subject_list}`
do
  echo ${sub}
  newsub=`echo ${sub} | tr -d _`
  if [[ -f ${julia_BIDS_dir}/sub-${newsub}/anat/sub-${newsub}_run-01_T1w.json ]]; then
  echo it exists
else
  echo it doesnt exist
  if [[ ${sub} == *_2 ]]; then
   echo 'its longitudinal T2'
   echo $newsub
   cp -r ${my_BIDS_dir}/sub-${sub:0:8}/ses-02/. ${julia_BIDS_dir}/sub-${newsub}/
 elif [[ ${sub} == *_3 ]]; then
     echo 'its longitudinal T3'
     echo $newsub
     cp -r ${my_BIDS_dir}/sub-${sub:0:8}/ses-03/. ${julia_BIDS_dir}/sub-${newsub}/
else
  echo 'its not longitudinal'
  cp -r ${my_BIDS_dir}/sub-${sub}/ses-01/. ${julia_BIDS_dir}/sub-${newsub}/
    fi
  fi
done

#remove accidental copies?
for sub in `find . -maxdepth 2 -iname 'ses-*'| sed -e 's|./||' | sed -e 's|/ses-01||'`

#taking out unnecessary ses-xx in filenames that is breaking BIDS validation, run on CBPD_bids
#Session 1
find . -not -path "./derivatives*" -type f -name 'sub-CBPD*_ses-01_*' | while read FILE; do
  echo ${FILE}
  newfile="$(echo ${FILE} | sed -e 's|_ses-01||')" ;
  echo ${newfile}
  mv "${FILE}" "${newfile}" ;
done
#session 2/3
find . -not -path "./derivatives*" -type f -name 'sub-CBPD*_ses-03_*' | while read FILE; do
  echo ${FILE}
  #add a 2 onto subject name, remove ses-02
  newfile="$(echo ${FILE} | sed -e 's|_ses-03||' | sed -e 's|[a-z0-9]/sub-CBPD....|&3|')" ;
  echo $newfile
  echo mv "${FILE}" "${newfile}"
done

#copy new freesurfer
for sub in `cat ${subject_list}`
do
  echo ${sub}
  newsub=`echo ${sub} | tr -d _`
  export SUBJECTS_DIR=${julia_BIDS_dir}/derivatives/freesurfer
  if grep -q RUNTIME_HOURS ${SUBJECTS_DIR}/sub-${newsub}/scripts/recon-all.done; then
    echo 'Freesurfer is already run for' ${sub}
  else
    echo 'Freesurfer not in directory for' ${sub}
    if [[ ${sub} == *_2 ]]; then
     echo 'its longitudinal T2'
     echo $newsub
     cp -r ${my_BIDS_dir}/derivatives/freesurfer_t2/sub-${sub:0:8}/. ${julia_BIDS_dir}/derivatives/freesurfer/sub-${newsub}/
   elif [[ ${sub} == *_3 ]]; then
       echo 'its longitudinal T3'
       echo $newsub
       cp -r ${my_BIDS_dir}/derivatives/freesurfer_t3/sub-${sub:0:8}/.  ${julia_BIDS_dir}/derivatives/freesurfer/sub-${newsub}/
  else
    echo 'its not longitudinal'
    cp -r ${my_BIDS_dir}/derivatives/freesurfer_t1/sub-${sub}/.  ${julia_BIDS_dir}/derivatives/freesurfer/sub-${newsub}/
      fi
  fi
done

#run MRIQC and fMRIprep and xcpEngine (once the first two are done)
subject_list=/data/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/subjectLists/n151_all_cbpd_one_or_more_nonsleep_rest_at_least_130_vols.txt
julia_BIDS_dir=/data/picsl/mackey_group/BPD/CBPD_bids
my_BIDS_dir=/data/picsl/mackey_group/CBPD/CBPD_bids
SCRIPTS_DIR=/data/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/code/ppc
subject_list_dir=/data/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/subjectLists
for sub in `cat ${subject_list}`
do
  echo ${sub}
  newsub=`echo ${sub} | tr -d _`
if [ -e ${julia_BIDS_dir}/derivatives/mriqc_fd_2_mm/sub-${newsub}_run-01_T1w.html ]; then
  echo 'MRIQC already run for' ${sub}
  echo ''
else
  echo 'Submitting MRIQC job for' ${sub}
  #qsub -j y ${SCRIPTS_DIR}/mriqc/mriqc_cmd.sh ${newsub}
fi
if [[ -e /data/picsl/mackey_group/BPD/CBPD_bids/derivatives/fmriprep/sub-${newsub}/func/sub-${newsub}_task-rest_run-01_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz ]]; then
  echo 'fmriprep exists for' ${sub}
  #check if xcpEngine is already run for this pipeline and this person.
  last_run=$(find ${julia_BIDS_dir}/sub-${newsub}/func/ -name "*task-rest_run-*" -type f| cut -d- -f5 | cut -d_ -f1 |sort -n | tail -n1)
  last_run=${last_run:1}
  echo $last_run
  for run in $(seq 1 $last_run); #loop through the runs
  do
  if [[ -e /data/picsl/mackey_group/BPD/CBPD_bids/derivatives/xcpEngine_nogsr_spkreg_fd0.5dvars1.75_drpvls/sub-${newsub}/run-0${run}/regress/sub-${newsub}_run-0${run}_residualised.nii.gz ]]; then
    echo 'Xcpengine run already for run' ${run} 'sub' ${sub}
  else
    echo rm /data/picsl/mackey_group/BPD/CBPD_bids/derivatives/xcpEngine_nogsr_spkreg_fd0.5dvars1.75_drpvls/sub-${newsub}/run-0${run}/ -R
    # echo 'Put run' ${run} 'in list for xcpEngine for' ${sub}
    # find ${julia_BIDS_dir} -type f | grep "${newsub}_task-rest_run-0${run}_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz\$"| while read fname; do
    # tmp=$(echo "$fname" | awk -F '_' '{print $5}' ) #this parses on underscores and pulls 'run-01'
    # fname_mnt=$(echo "$fname" | sed -e 's|/data/|/mnt/|' )
    # echo sub-${newsub},${tmp},${fname_mnt} >> ${subject_list_dir}/nxxx_cleanup_for_xcpEngine.csv
    # done;
  fi
done
else
  echo 'No fMRIprep'
  #qsub -j y ${SCRIPTS_DIR}/fmriprep/fmriprep_cmd.sh sub-${newsub}
fi
done

#Do it by run!
SCRIPTS_DIR=/data/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/code/ppc
pipeline=xcpEngine_gsr_spkreg_fd0.5dvars1.75_drpvls
subject_list=/data/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/subjectLists/n150_cross_sect_one_or_more_nonsleep_rest_at_least_130_vols.csv
for line in `cat ${subject_list}`
do
  sub=`echo $line | cut -d, -f1`
  echo ${sub}
  run=`echo $line | cut -d, -f2| tr -d '\r\n'` #or just \n if it isn't working
  echo ${run}
  newsub=${sub}
  if [[ -f /data/picsl/mackey_group/BPD/CBPD_bids/derivatives/fmriprep/${sub}/func/${sub}_task-rest_${run}_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz ]]; then
    echo 'fmriprep exists for' ${sub} ${run}
    if [[ -e /data/picsl/mackey_group/BPD/CBPD_bids/derivatives/${pipeline}/${sub}/${run}/fcon/schaefer400x7/${sub}_${run}_schaefer400x7_network.txt ]]; then
      echo 'Xcpengine run already for run' ${run} 'sub' ${sub}
    elif [[ -e /data/picsl/mackey_group/BPD/CBPD_bids/derivatives/${pipeline}/${sub}/${run}/fcon/schaefer400/${sub}_${run}_schaefer400_network.txt ]]; then
      echo 'Xcpengine run already for run' ${run} 'sub' ${sub}
    else
      echo rm /data/picsl/mackey_group/BPD/CBPD_bids/derivatives/${pipeline}/${sub}/${run}/ -R
      # echo 'Put run' ${run} 'in list for xcpEngine for' ${sub}
      # find ${julia_BIDS_dir} -type f | grep "${sub}_task-rest_${run}_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz\$"| while read fname; do
      # tmp=$(echo "$fname" | awk -F '_' '{print $5}' ) #this parses on underscores and pulls 'run-01'
      # fname_mnt=$(echo "$fname" | sed -e 's|/data/|/mnt/|' )
      # echo ${sub},${tmp},${fname_mnt} >> ${subject_list_dir}/nxxx_ursula_pipeline_cleanup_for_xcpEngine.csv
      # done;
    fi
  else
    echo 'No fMRIprep'
    #qsub -j y ${SCRIPTS_DIR}/fmriprep/fmriprep_cmd.sh sub-${newsub}
  fi
done
