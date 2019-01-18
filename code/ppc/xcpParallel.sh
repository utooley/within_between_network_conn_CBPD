#$ -j y
#$ -o /data/picsl/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/output/qsub_output
#$ -t 1-22

# Adjust these so they work on your system
SNGL=/share/apps/singularity/2.5.1/bin/singularity
SIMG=/data/picsl/mackey_group/tools/singularity/xcpEngine.simg
FULL_COHORT=/data/picsl/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/subjectLists/n47_second_batch_cohort_usable_t1_rest_1mm_outliers_10_2mm_11718.csv

# Create a temp cohort file with 1 line
HEADER=$(head -n 1 $FULL_COHORT)
LINE_NUM=$( expr $SGE_TASK_ID + 1 )
LINE=$(awk "NR==$LINE_NUM" $FULL_COHORT)
TEMP_COHORT=${FULL_COHORT}.${SGE_TASK_ID}.csv
echo $HEADER > $TEMP_COHORT
echo $LINE >> $TEMP_COHORT

$SNGL run -B /data:/mnt $SIMG \
  -c /mnt${TEMP_COHORT#/data} \
  -d /mnt/picsl/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/code/ppc/fc-36p-nogsr-meancsfwm.dsn \
  -o /mnt/picsl/mackey_group//BPD/CBPD_bids/derivatives/xcpEngine_nogsr_nospkreg \
  -i $TMPDIR

