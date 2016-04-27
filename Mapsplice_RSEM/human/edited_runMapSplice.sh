#!/bin/bash
#shell script to run create directories for each samples and run MapSplice2 to align fastq files
#usage: sh runMapSplice.sh <fastq files>   ##for TCGA use *_1.fastq
#outputs: alignments.bam, junctions.txt, insertions.txt, deletions.txt, fusions_raw.txt, fusions_candidates.txt, stats.txt

for ARG in "$@"
do

BASE=`basename $ARG .fastq`
Abase=`basename $BASE _001.fastq` #make sure endings same if using this script for your files
A2=${Abase}2_001.fastq
DIR=`dirname $ARG`
OUTPUT_DIR=$DIR/$Abase
MAPSPLICE_DIR=/nas02/apps/mapsplice-2.1.4/src/MapSplice-v2.1.4
REF_GENOME=/proj/seq/data/HG19_UCSC/Sequence/Chromosomes/ #make sure ref_genome and bowtie_index are correct for you
BOWTIE_INDEX=/proj/seq/data/HG19_UCSC/Sequence/BowtieIndex/genome

mkdir -p $OUTPUT_DIR

echo $BASE
echo $Abase
echo $A2

bsub -o $OUTPUT_DIR/ms.kure.log -n 4 -R "span[hosts=1]" "python $MAPSPLICE_DIR/mapsplice.py  -1 $ARG  -2 $A2  -c $REF_GENOME  -x $BOWTIE_INDEX  -p 4 --fusion-non-canonical --bam  --qual-scale phred33  -o $OUTPUT_DIR 2>$OUTPUT_DIR/$Abase.ms.log"

done

#TCGA usage
#python mapsplice_multi_thread.py --fusion --all-chromosomes-files hg19_M_rCRS/hg19_M_rCRS.fa --pairend -X 8 -Q fq --chromosome-files- dir hg19_M_rCRS/chromosomes --Bowtieidx hg19_M_rCRS/ebwt/humanchridx_M_rCRS -1 working/prep_1.fastq -2 working/prep_2.fastq -o SAMPLE_BARCODE 2> working/mapsplice.log


#sara selitsky's perl script to do the same thing:
#$mapSpliceDir = "/nas02/apps/mapsplice-2.1.4/src/MapSplice-v2.1.4";
#$reference = "/proj/seq/data/HG19_UCSC/Sequence/Chromosomes/";
#$bowtieIndex = "/proj/seq/data/HG19_UCSC/Sequence/BowtieIndex/genome";
#
#@file_1 = @ARGV;
#
#for($i=0;$i<scalar(@file_1);$i++){
#	$nameSub1 = $file_1[$i];
#	$nameSub2 = $file_1[$i];
#	$nameSub2 =~ s/_R1_/_R2_/;
#	$nameDir = $nameSub1;
#	$nameDir =~ s/.fastq//;
#	$nameDir =~ s/_R1_001//;
#	`mkdir -p $nameDir`;
#
#	`bsub -n 4 -R "span[hosts=1]" "python $mapSpliceDir/mapsplice.py -1 $nameSub1 -2 $nameSub2 -c $reference -x $bowtieIndex -p 4 --non-canonical --fusion-non-canonical --bam --qual-scale phred33 -o $nameDir 2>$nameDir/$nameDir.log"`;
#}
