#$ -j y
#$ -l h_vmem=200.1G,s_vmem=200.0G
#$ -o /cbica/home/tooleyu/projects/in_progress/within_between_network_conn_CBPD/output/qsub_output
#$ -t 1-10

# Adjust these so they work on your system
SNGL=/usr/bin/singularity
SIMG=/cbica/projects/cbpd_main_data/tools/singularity/xcpEngine-100219.simg
FULL_COHORT=/cbica/home/tooleyu/projects/in_progress/within_between_network_conn_CBPD/data/subjectLists/n11_fd0.25_pipeline_reprocess_revisions_for_xcpEngine_cleanup2.csv

# Create a temp cohort file with 1 line
HEADER=$(head -n 1 $FULL_COHORT)
LINE_NUM=$( expr $SGE_TASK_ID + 1 )
LINE=$(awk "NR==$LINE_NUM" $FULL_COHORT)
TEMP_COHORT=${FULL_COHORT}.${SGE_TASK_ID}.csv
echo $HEADER > $TEMP_COHORT
echo $LINE >> $TEMP_COHORT

echo 'temporary directory is ' $TMPDIR
sleep 91
$SNGL run --cleanenv --env R_PROFILE_USER=usr,R_ENVIRON_USER=usr -B $TMPDIR:/tmp $SIMG \
  -c ${TEMP_COHORT} \
  -d /cbica/home/tooleyu/projects/in_progress/within_between_network_conn_CBPD/code/ppc/fc-36p-gsr-meancsfwm-spkreg-fd0.25-dv1.75-dropvols.dsn \
  -o /cbica/projects/cbpd_main_data/CBPD_bids_crosssectional/derivatives/xcpEngine_gsr_spkreg_fd0.25dvars1.75_drpvls \
  -i /cbica/projects/cbpd_main_data/CBPD_bids_crosssectional/derivatives/xcpengine_wd_2 \
  -t 2

