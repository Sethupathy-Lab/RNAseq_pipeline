#!/bin/bash
#shell script for running FastQC to quality check fastq files
#usage: sh qc.sh <fastq files>
#output: kure log file for each fastq file, .zip & .html FastQC output for each fastq file

mkdir -p FastQC/logs

for ARG in "$@"
do

BASE=`basename $ARG .fastq`
DIR=`dirname $ARG`

bsub -o $DIR/FastQC/logs/qc.$BASE.log fastqc -o $DIR/FastQC $ARG

done

