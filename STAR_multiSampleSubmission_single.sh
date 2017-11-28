#!/bin/bash

# Parameters that need to be changed
genome=/home/pr46_0001/projects/genome/gencode.vM10

# Make output directories
rmdir STAR_OUT
rmdir SALMON_OUT
rmdir STAR_Salmon_logs

mkdir STAR_OUT
mkdir SALMON_OUT
mkdir STAR_Salmon_logs

# Loop to submit samples to STAR
for fastq in *.fastq.gz; do 

# Getting information for job submission
  READ1=$fastq
  DIR=$(dirname $fastq)
  SAMP=${fastq//fastq.gz/}
  OUT_DIR=${DIR}/STAR_OUT/$SAMP
  LOGS=${DIR}/STAR_Salmon_logs

# Print out information about job submission
  echo Base name is: $SAMP
  echo Read 1 is: $READ1
  echo Output directory is: $OUT_DIR
  echo Logs are in: $LOGS

# Job submission
/programs/STAR_2.4.2a/STAR \
  --runThreadN 8 \
  --genomeDir $genome \
  --readFilesIn $READ1 \
  --outFilterMismatchNmax 20 \
  --outSAMtype BAM Unsorted \
  --quantMode TranscriptomeSAM \
  --outFileNamePrefix STAR_OUT/$SAMP \
  --outSAMunmapped Within \
  --readFilesCommand zcat \
  > STAR_Salmon_logs/${SAMP}.out

 echo
done

