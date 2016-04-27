#!/bin/bash
#calls next.sh to process each aligned bam file given as argument
#usage: sh nextOuter.sh */alignments.bam               ###sara's usage: */*[0-9B]_alignments.bam 
#outputs: 

for ARG in "$@"
do

   DIR=`dirname $ARG`
   LogNM=$DIR/next.kure.log
   bsub -o $LogNM sh next.sh $ARG
done


