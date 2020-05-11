declare -a array=(sub-CBPD0141
sub-CBPD0146
sub-CBPD0150
sub-CBPD0157
sub-CBPD0166
sub-CBPD00562
sub-CBPD00772
sub-CBPD0170
sub-CBPD0171)

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

#copy over new data from the clean CBPD bids directory.

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

#run MRIQC and fMRIprep
SCRIPTS_DIR=/data/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/code/ppc
for sub in `cat ${subject_list}`
do
  echo ${sub}
  newsub=`echo ${sub} | tr -d _`
if [ -e ${julia_BIDS_dir}/derivatives/mriqc_fd_2_mm/sub-${newsub}_run-01_T1w.html ]; then
  echo 'MRIQC already run with 2 mm threshold for' ${sub}
else
  echo 'Running MRIQC with 2 mm threshold for' ${sub}
  qsub -j y ${SCRIPTS_DIR}/mriqc/mriqc_cmd.sh ${newsub}
fi
if [[ -d /data/picsl/mackey_group/BPD/CBPD_bids/derivatives/fmriprep/sub-${newsub} ]]; then
echo 'fMRIprep already exists for' ${sub}
else
  echo 'No fMRIprep'
  qsub -j y ${SCRIPTS_DIR}/fmriprep/fmriprep_cmd.sh sub-${newsub}
fi
done

#run xcpEngine pipelines.
