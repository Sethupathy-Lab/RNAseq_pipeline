#!/bin/bash

# Parameters that need to be changed
transcriptome=/home/pr46_0001/projects/genome/gencode.vM10/GRCm38.p4.transcripts.fa
LOGS=STAR_Salmon_logs

for alignment in STAR_OUT/*Aligned.toTranscriptome.out.bam; do 

# Get info for job submission
  NAME=$(basename $alignment Aligned.toTranscriptome.out.bam)
  OUT_DIR=SALMON_OUT/${NAME}

# Print info for job submission
  echo Base name is: $NAME
  echo Output directory is: $OUT_DIR

# Job submission (The -l will have to be changed according to type of sample prep)
  /home/pr46_0001/shared/bin/salmon quant \
  -t $transcriptome \
  -l A \
  -p 8 \
  -a STAR_OUT/${NAME}Aligned.toTranscriptome.out.bam \
  -o $OUT_DIR

  echo
done

