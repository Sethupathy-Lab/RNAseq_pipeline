#!/bin/bash

for ARG in "$@"
do

BASE=`basename $ARG .fastq`
DIR=`dirname $ARG`

bsub -o qc.log fastqc $ARG

done

