# RNAseq analysis pipeline
Last updated: Nov 2017

### Our current preferred pipeline involves STAR/Salmon alignment and quantification followed by DESeq2 for differential expression analysis. (Differential expression may be changing to Sleuth in coming months)

### Preamble

While the code below will get you going on your RNAseq analysis, I would greatly encourage looking at program documentation or tutorials for the software used.

[STAR documentation](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)  
[Salmon documentation](http://salmon.readthedocs.io/en/latest/salmon.html)  
[DESeq2 documentation](http://bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf)  
[DESeq2 tutorial](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)  

#### Building indexes
Need to update how to build STAR indexes

#### Running STAR

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

#### Running Salmon

The Salmon script can infer whether single or paired-end reads were used, so the same script can be used. The script looks like this:
```
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
```

To run:

1. Move Salmon_MultiSampleSubmission.sh to directory containing fastqs

2. Edit the location transcriptome fasta. We are using mouse, but human is also available.  
mouse transcriptome: `/home/pr46_0001/projects/genome/gencode.vM10/GRCm38.p4.transcripts.fa`  
human transcriptome: `/home/pr46_0001/projects/genome/GRCh38.p7/gencode.v25.transcripts.fa`

3. Run SALMON_MultiSampleSubmission.sh
```
bash SALMON_MultiSampleSubmission.sh
```

Once STAR/Salmon is done running, we have our counts. We can now transfer the counts to our local machine to run the differential analysis.

