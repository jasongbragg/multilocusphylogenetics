#!/bin/bash
#$ -cwd
#$ -m ae
#$ -r y
#$ -N mraic
#$ -t 1-421:1
#$ -tc 40

DIR=/home/jgb/mraicrun/
nexfils=( ${DIR}*.phy  )
filnum=$[$SGE_TASK_ID-1]
parfil=${nexfils[$filnum]}

echo $SGE_TASK_ID
echo $parfil

MRAIC=/home/jgb/bin/mraic.pl
cd $DIR
perl ${MRAIC} $parfil
