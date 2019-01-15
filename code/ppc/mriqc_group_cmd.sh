MACKEY_HOME=/data/picsl/mackey_group/
#BIDS_folder=/data/picsl/mackey_group/BPD/niftis
BIDS_folder=/data/picsl/mackey_group/BPD/CBPD_bids
working_dir=/data/picsl/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/mriqc_working
tools_dir=${MACKEY_HOME}/tools/singularity
output_dir=${BIDS_folder}/derivatives/mriqc

unset PYTHONPATH;
singularity run --cleanenv -B ${BIDS_folder}:/mnt ${tools_dir}/mriqc-0.14.2.simg \
/mnt/ /mnt/derivatives/mriqc \
group \
-w /tmp \


#don't have the working dir be a mounted directory or you'll get errors in the logs about busy resources
#must bind to a folder that already exists in the container, and must point to data dir not subject dir
