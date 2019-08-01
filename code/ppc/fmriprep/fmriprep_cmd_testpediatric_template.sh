#$ -j y
#$ -l h_vmem=25.1G,s_vmem=25G
#$ -o /data/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/output/qsub_output
#$ -q himem.q,all.q,basic.q,gpu.q

MACKEY_HOME=/data/picsl/mackey_group
BIDS_folder=/data/picsl/mackey_group/BPD/CBPD_bids
working_dir=/data/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/fmriprep_working
subject=CBPD0141
tools_dir=${MACKEY_HOME}/tools/singularity
output_dir=${BIDS_folder}/derivatives/

unset PYTHONPATH;
echo 'job is running'
#export SINGULARITYENV_TEMPLATEFLOW_HOME=/home/utooley/templateflow
export SINGULARITYENV_TEMPLATEFLOW_HOME=/home/CBLuser/templateflow
#singularity run --cleanenv -B /home/utooley/templateflow:/home/utooley/templateflow,${BIDS_folder}:/mnt ${tools_dir}/fmriprep-1-4-1rc5.simg \
singularity run --cleanenv -B /home/CBLuser/templateflow:/home/CBLuser/templateflow,${BIDS_folder}:/mnt ${tools_dir}/fmriprep-1-4-1rc5.simg \
/mnt/ /mnt/derivatives/fmriprep_test_2 \
participant \
-w /tmp/utooley${subject} \
--participant-label ${subject} \
--fs-license-file $HOME/license.txt \
--skull-strip-template MNIPediatricAsym:cohort-2 \
--output-spaces MNIPediatricAsym:cohort-2 T1w MNI152NLin2009cAsym fsaverage5 \
--nthreads 1 \



#--skip-bids-validation \
#--templace-resampling-grid native \
#must bind to a folder that already exists in the container, and must point to data dir not subject dir
