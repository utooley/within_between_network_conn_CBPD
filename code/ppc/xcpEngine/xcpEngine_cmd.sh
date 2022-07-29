#$ -j y
#$ -l h_vmem=200.1G,s_vmem=200.0G
#$ -o /cbica/home/tooleyu/projects/in_progress/within_between_network_conn_CBPD/output/qsub_output


# Adjust these so they work on your system
SNGL=/usr/bin/singularity
SIMG=/cbica/projects/cbpd_main_data/tools/singularity/xcpEngine-100219.simg
FULL_COHORT=/cbica/home/tooleyu/projects/in_progress/within_between_network_conn_CBPD/data/subjectLists/n1_fd0.25_pipeline_reprocess_revisions_for_xcpEngine_cleanup2.csv

echo 'temporary directory is ' $TMPDIR
#sleep $((RANDOM % 45))
$SNGL run --cleanenv --env R_PROFILE_USER=usr,R_ENVIRON_USER=usr -B $TMPDIR:/tmp $SIMG \
 -c ${FULL_COHORT} \
 -d /cbica/home/tooleyu/projects/in_progress/within_between_network_conn_CBPD/code/ppc/fc-36p-gsr-meancsfwm-spkreg-fd0.25-dv1.75-dropvols.dsn \
 -o /cbica/projects/cbpd_main_data/CBPD_bids_crosssectional/derivatives/xcpEngine_gsr_spkreg_fd0.25dvars1.75_drpvls \
 -i /cbica/projects/cbpd_main_data/CBPD_bids_crosssectional/derivatives/xcpengine_wd_2 \
 -t 2

# MACKEY_HOME=/data/picsl/mackey_group/
# #BIDS_folder=/data/picsl/mackey_group/BPD/niftis
# #BIDS_folder=/data/picsl/mackey_group/BPD/CBPD_bids
# project_dir=picsl/mackey_group//Ursula/projects/in_progress/within_between_network_conn_CBPD/
# subject=${1}
# tools_dir=${MACKEY_HOME}/tools/singularity
# output_dir=${BIDS_folder}/derivatives/xcpEngine_nogsr
#
# unset PYTHONPATH;
# singularity run --cleanenv -B /data:/mnt  \
#   ${tools_dir}/xcpEngine.simg \
#   -d /mnt/${project_dir}/code/ppc/test.dsn \
#   -c /mnt/${project_dir}/data/subjectLists/test_cohort_2subs.csv  \
#   -o /mnt/picsl/mackey_group/BPD/CBPD_bids/derivatives/ \
#   -t 1 \
