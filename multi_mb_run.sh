#!/bin/bash
#$ -cwd
#$ -m ae
#$ -r y
#$ -N mrbayes
#$ -t 1-413:1
#$ -tc 40

DIR=/home/jgb/mbgenetrees/
nexfils=( ${DIR}*params.txt  )
filnum=$[SGE_TASK_ID-1]
parfil=${nexfils[$filnum]}

echo $SGE_TASK_ID
echo $parfil

cd ${DIR}
/opt/openmpi/bin/mpirun -np 1 /opt/bio/mrbayes/mb < ${parfil} > ${parfil}.log
