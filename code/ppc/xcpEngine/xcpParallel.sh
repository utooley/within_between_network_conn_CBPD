#$ -j y
#$ -l h_vmem=95.1G,s_vmem=95.0G
#$ -o /cbica/projects/cbpd_main_data/qsub_output
#$ -t 1-3

# Adjust these so they work on your system
SNGL=/usr/bin/singularity
SIMG=/cbica/projects/cbpd_main_data/tools/singularity/xcpEngine-100219.simg
FULL_COHORT=/cbica/home/tooleyu/projects/in_progress/within_between_network_conn_CBPD/data/subjectLists/n3_anne_cohort_piper_xcpEngine.csv

# Create a temp cohort file with 1 line
HEADER=$(head -n 1 $FULL_COHORT)
LINE_NUM=$( expr $SGE_TASK_ID + 1 )
LINE=$(awk "NR==$LINE_NUM" $FULL_COHORT)
TEMP_COHORT=${FULL_COHORT}.${SGE_TASK_ID}.csv
echo $HEADER > $TEMP_COHORT
echo $LINE >> $TEMP_COHORT

echo 'temporary directory is ' $TMPDIR
$SNGL run --cleanenv --env R_PROFILE_USER=usr,R_ENVIRON_USER=usr -B $TMPDIR:/tmp $SIMG \
  -c ${TEMP_COHORT} \
  -d /cbica/home/tooleyu/projects/in_progress/within_between_network_conn_CBPD/code/ppc/fc-36p-nogsr-meancsfwm-spkreg-dropvols.dsn \
  -o /cbica/projects/cbpd_main_data/CBPD_bids_crosssectional/derivatives/xcpEngine_nogsr_spkreg_fd1.25dvars2_drpvls \
  -i /cbica/projects/cbpd_main_data/CBPD_bids_crosssectional/derivatives/xcpengine_wd_2 \
  -t 2
