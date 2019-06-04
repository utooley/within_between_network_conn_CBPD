MACKEY_HOME=/data/picsl/mackey_group/
#BIDS_folder=/data/picsl/mackey_group/BPD/niftis
BIDS_folder=/data/picsl/mackey_group/BPD/CBPD_bids
tools_dir=${MACKEY_HOME}/tools/singularity
output_dir=${BIDS_folder}/derivatives/mriqc_fd_.5_mm

unset PYTHONPATH;
singularity run --cleanenv -B ${BIDS_folder}:/mnt ${tools_dir}/mriqc-0.14.2.simg \
/mnt/ /mnt/derivatives/mriqc_fd_.5_mm \
group \
-w ${TMP} \


#don't have the working dir be a mounted directory or you'll get errors in the logs about busy resources
#must bind to a folder that already exists in the container, and must point to data dir not subject dir
