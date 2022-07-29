#!/bin/bash
MACKEY_HOME=/cbica/projects/cbpd_main_data
project_dir=/cbica/home/tooleyu/projects/in_progress/within_between_network_conn_CBPD
sublist_dir=${project_dir}/data/subjectLists
tools_dir=${MACKEY_HOME}/tools/singularity
#tools_dir=/data/joy/BBL/applications/bids_apps #eventually check on that this hasn't changed.

FULL_COHORT=${sublist_dir}/n11_fd0.25_pipeline_reprocess_revisions_for_xcpEngine_cleanup2.csv

#FULL_COHORT=${sublist_dir}/n2_test.csv
NJOBS=`wc -l < $FULL_COHORT`

if [[ ${NJOBS} == 0 ]]; then
    echo 'you dont have enough lines in your csv file'
    exit 0
fi
#NJOBS=2

cat << EOF > xcpParallel.sh
#$ -j y
#$ -l h_vmem=200.1G,s_vmem=200.0G
#$ -o /cbica/home/tooleyu/projects/in_progress/within_between_network_conn_CBPD/output/qsub_output
#$ -t 1-${NJOBS}

# Adjust these so they work on your system
SNGL=/usr/bin/singularity
SIMG=${tools_dir}/xcpEngine-100219.simg
FULL_COHORT=${FULL_COHORT}

# Create a temp cohort file with 1 line
HEADER=\$(head -n 1 \$FULL_COHORT)
LINE_NUM=\$( expr \$SGE_TASK_ID + 1 )
LINE=\$(awk "NR==\$LINE_NUM" \$FULL_COHORT)
TEMP_COHORT=\${FULL_COHORT}.\${SGE_TASK_ID}.csv
echo \$HEADER > \$TEMP_COHORT
echo \$LINE >> \$TEMP_COHORT

echo 'temporary directory is ' \$TMPDIR
sleep $((RANDOM % 100))
\$SNGL run --cleanenv --env R_PROFILE_USER=usr,R_ENVIRON_USER=usr -B \$TMPDIR:/tmp \$SIMG \\
  -c \${TEMP_COHORT} \\
  -d ${project_dir}/code/ppc/fc-36p-gsr-meancsfwm-spkreg-fd0.25-dv1.75-dropvols.dsn \\
  -o ${MACKEY_HOME}/CBPD_bids_crosssectional/derivatives/xcpEngine_gsr_spkreg_fd0.25dvars1.75_drpvls \\
  -i ${MACKEY_HOME}/CBPD_bids_crosssectional/derivatives/xcpengine_wd_2 \\
  -t 2

EOF

qsub xcpParallel.sh
