MACKEY_HOME=/data/picsl/mackey_group/
#BIDS_folder=/data/picsl/mackey_group/BPD/niftis
BIDS_folder=/data/picsl/mackey_group/BPD/CBPD_bids
working_dir=/data/picsl/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/mriqc_working
subject=${1:4}
tools_dir=${MACKEY_HOME}/tools/singularity
output_dir=${BIDS_folder}/derivatives/mriqc_fd_2_mm
echo $subject

unset PYTHONPATH;
singularity run --cleanenv -B ${BIDS_folder}:/mnt ${tools_dir}/mriqc-0.14.2.simg \
/mnt/ /mnt/derivatives/mriqc_fd_2_mm \
participant \
-w /tmp \
--participant_label ${subject} \
--fd_thres 2 \
--no-sub \

#must bind to a folder that already exists in the container, and must point to data dir not subject dir
