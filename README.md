# RNAseq analysis pipeline

## Our current preferred pipeline involves STAR/Salmon alignment and quantification followed by DESeq2 for differential expression analysis. (Differential expression may be changing to Sleuth in coming months)

#### Building indexes
Need to update how to build STAR indexes

#### Running STAR on paired-end sequences

The STAR_multisample_single.sh is for single-read sequencing, whereas the STAR_multisample_paired.sh is for paired-end sequencing. Both are *very* similar, and you can look at both to see the differences, but below is the STAR_multisample_paired.sh script:

```
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
```

To run:

1. Move the STAR_multisubmission script to directory containing sequencing files

2. Make sure that your fastq files end in the extension `.fastq.gz`. Our STAR command expects our fastqs to be gzipped, so if you don't have gzipped fastqs, run `for x in *.fastq; do gzip $x; done`

3. Edit the genome location line for the species you are aligning to. In the example above we are using mouse, but we have human set up as well.

mouse location: `/home/pr46_0001/projects/genome/gencode.vM10`
human location: `/home/pr46_0001/projects/genome/GRCh38.p7`

4. Run STAR_multisubmission_(single or paired).sh (eg for single):
```
bash STAR_multisubmission_single.sh
```

