#!/bin/bash

for ARG in "$@"
do

   DIR=`dirname $ARG`
   LogNM=$DIR.next.log
   bsub -o $LogNM sh next.sh $ARG
done


