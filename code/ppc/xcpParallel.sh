#$ -j y
#$ -q all.q,gpu.q,himem.q
#$ -l h_vmem=95.1G,s_vmem=95.0G
#$ -o /data/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/output/qsub_output
#$ -t 1-3

# Adjust these so they work on your system
SNGL=/share/apps/singularity/2.5.1/bin/singularity
SIMG=/data/picsl/mackey_group/tools/singularity/xcpEngine-100219.simg
FULL_COHORT=/data/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/subjectLists/nxxx_ursula_pipeline_cleanup_for_xcpEngine.csv

# Create a temp cohort file with 1 line
HEADER=$(head -n 1 $FULL_COHORT)
LINE_NUM=$( expr $SGE_TASK_ID + 1 )
LINE=$(awk "NR==$LINE_NUM" $FULL_COHORT)
TEMP_COHORT=${FULL_COHORT}.${SGE_TASK_ID}.csv
echo $HEADER > $TEMP_COHORT
echo $LINE >> $TEMP_COHORT

echo 'temporary directory is ' $TMPDIR
$SNGL run --cleanenv -B /data:/mnt,$TMPDIR:/tmp $SIMG \
  -c /mnt${TEMP_COHORT#/data} \
  -d /mnt/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/code/ppc/fc-36p-nogsr-meancsfwm-spkreg-fd0.5-dv1.75-dropvols.dsn \
  -o /mnt/picsl/mackey_group//BPD/CBPD_bids/derivatives/xcpEngine_nogsr_spkreg_fd0.5dvars1.75_drpvls \
  -i /mnt/picsl/mackey_group//BPD/CBPD_bids/derivatives/xcpengine_wd_2 \
  -t 2

