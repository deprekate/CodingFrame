#!/bin/bash 

#$ -S /bin/bash
#$ -cwd

#$ -t 1-4374

SEEDFILE=list.tsv
SEED=$(awk "NR==$SGE_TASK_ID" $SEEDFILE)

#$ -o sge/out-$JOB_ID-$TASK_ID
#$ -e sge/err-$JOB_ID-$TASK_ID


python3 codingframe.py $SEED > genomes/aa/$SEED
