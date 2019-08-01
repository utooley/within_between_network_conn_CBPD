#$ -j y
#$ -q all.q,basic.q,himem.q,gpu.q
#$ -l h_vmem=19.1G,s_vmem=19.0G
#$ -o /data/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/output/qsub_output
#$ -t 1-99

# Adjust these so they work on your system
SNGL=/share/apps/singularity/2.5.1/bin/singularity
SIMG=/data/picsl/mackey_group/tools/singularity/xcpEngine-060319.simg
FULL_COHORT=/data/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/data/subjectLists/n76_cohort_mult_runs_usable_t1_rest_1mm_outliers_10_2mm_80119.csv

# Create a temp cohort file with 1 line
HEADER=$(head -n 1 $FULL_COHORT)
LINE_NUM=$( expr $SGE_TASK_ID + 1 )
LINE=$(awk "NR==$LINE_NUM" $FULL_COHORT)
TEMP_COHORT=${FULL_COHORT}.${SGE_TASK_ID}.csv
echo $HEADER > $TEMP_COHORT
echo $LINE >> $TEMP_COHORT

$SNGL run --cleanenv -B /data:/mnt $SIMG \
  -c /mnt${TEMP_COHORT#/data} \
  -d /mnt/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/code/ppc/fc-36p-nogsr-meancsfwm-spkreg-dropvols.dsn \
  -o /mnt/picsl/mackey_group//BPD/CBPD_bids/derivatives/xcpEngine_nogsr_spkreg_fd1.25dvars2_drpvls \
  -i $TMPDIR

