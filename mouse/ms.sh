#!/bin/bash

for ARG in "$@"
do

BASE=`basename $ARG .fastq`
Abase=`basename $BASE _R1_001`
A2=${Abase}_R2_001.fastq
DIR=`dirname $ARG`
OUTPUT_DIR=$DIR/$Abase
MAPSPLICE_DIR=/nas02/apps/mapsplice-2.1.4/src/MapSplice-v2.1.4
REF_GENOME=/proj/seq/data/MM10_UCSC/Sequence/Chromosomes/
BOWTIE_INDEX=/proj/seq/data/MM10_UCSC/Sequence/BowtieIndex/genome

mkdir -p $OUTPUT_DIR


bsub -o ms.log -n 4 -R "span[hosts=1]" "python $MAPSPLICE_DIR/mapsplice.py  -1 $ARG  -2 $A2  -c $REF_GENOME  -x $BOWTIE_INDEX  -p 4 --fusion --bam  --qual-scale phred33  -o $OUTPUT_DIR 2>$OUTPUT_DIR/$Abase.log"

done
#python mapsplice_multi_thread.py --fusion --all-chromosomes-files hg19_M_rCRS/hg19_M_rCRS.fa --pairend -X 8 -Q fq --chromosome-files- dir hg19_M_rCRS/chromosomes --Bowtieidx hg19_M_rCRS/ebwt/humanchridx_M_rCRS -1 working/prep_1.fastq -2 working/prep_2.fastq -o SAMPLE_BARCODE 2> working/mapsplice.log
