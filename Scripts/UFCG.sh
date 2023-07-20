#! /bin/bash

UFCG_ENV="/home/sebr/seb/Applications/miniconda3/envs/ufcg"
UFCG="/home/sebr/seb/Applications/UFCG/UFCG.jar"
CONDA="/mnt/gpfs/seb/Applications/miniconda3/bin/activate"

source $CONDA $UFCG_ENV

fail=1
cnt=1
while ((fail!=0))
do
	echo $cnt
	java -jar $UFCG profile -i $1 -t $2  -w $3 -o $4
	fail=$?
	((cnt+=1))
done
