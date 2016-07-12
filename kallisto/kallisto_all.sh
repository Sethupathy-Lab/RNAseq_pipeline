#!/bin/bash

###############################################################################
#
#  Run kallisto alignment and quantification for all paired ended samples
#
###############################################################################

module load kallisto

# Edit these parameters accordingly
transcriptome=rn5
projdir=/proj/seth_lab/projects/RNASeq/UCD_T2DM_LIVER_RAT
index=/proj/seth_lab/projects/genome/$transcriptome/${transcriptome}.kallisto.idx

cd ${projdir}
ls *_R1_*.fastq.gz>files.txt
mkdir results

while read line
do 
 f1=$line
 f2=${f1//R1/R2}
 dname=`echo ${f1}|cut -d '_' -f 1`
 mkdir ${projdir}/results/${dname}
 mkdir ${projdir}/results/${dname}/kallisto  
 out=${projdir}/results/${dname}/kallisto
 echo 'running '${f1}':'${f2}
 bsub -o ${projdir}/results/${dname}.log kallisto quant -b 30 -i ${index} -o ${out} ${projdir}/${f1} ${projdir}/${f2}
done<${projdir}/files.txt
