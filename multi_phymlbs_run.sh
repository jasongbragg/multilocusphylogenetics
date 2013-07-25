#!/bin/bash
#$ -cwd
#$ -m ae
#$ -r y
#$ -N phyml 
#$ -t 1-413:1
#$ -tc 40

DIR=/home/jgb/phymlbs/
nexfils=( ${DIR}*.phy.phymlbsparams  )
filnum=$[$SGE_TASK_ID-1]
parfil=${nexfils[$filnum]}

echo $SGE_TASK_ID
echo $parfil

cd $DIR
/home/jgb/bin/phyml < $parfil
