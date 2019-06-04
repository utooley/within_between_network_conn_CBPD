#!/bin/sh

#Loop script for FMRIPREP
#$ -j y
#$ -t 1-15
#$ -l h_vmem=20.1G,s_vmem=20G
#$ -q all.q,basic.q,himem.q,gpu.q
#$ -o /data/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/output/qsub_output
unset PYTHONPATH;
URSULA_PROJ=/data/jux/mackey_group/Ursula/projects/in_progress
ERROR_DIR=/data/jux/mackey_group/Ursula/projects/in_progress/within_between_network_conn_CBPD/output/qsub_output
SUBLIST_DIR=${URSULA_PROJ}/within_between_network_conn_CBPD/data/subjectLists
SUB_FILE=${SUBLIST_DIR}/n15_add_subjs_usable_t1_norestqualfiltering_053019
SCRIPTS_DIR=${URSULA_PROJ}/within_between_network_conn_CBPD/code/ppc/mriqc

mapfile -t ARRAY < ${SUB_FILE}

LENGTH=${#ARRAY[@]}

i=`expr ${SGE_TASK_ID} - 1`

if [[ $i -ge ${LENGTH} ]]; then
 echo 'Array index > than number of elements'

else
 SUB=${ARRAY[$i]}

bash ${SCRIPTS_DIR}/mriqc_cmd.sh ${SUB}

fi
