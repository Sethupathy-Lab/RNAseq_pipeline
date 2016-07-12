# RNAseq analyses pipelines

## Current preferred RNAseq mapping/alignment: kallisto

#### Building indexes
Kallisto needs to index the transcriptome prior to psuedo-aligment and quantification.
To index the transcriptome, we must:
```
module load kallisto
bsub kallisto index -i transcriptome.kallisto.idx transcriptome.fa

where:
  -i transcriptome.kallisto.idx = name of kallisto index file (output)
  transcriptome.fa = name of transcriptome (input)
```

#### Running kallisto on paired-end sequences

Running kallisto is pretty straight forward, but to make it a bit easier, there is a bash script.  The script needs to be edited in the 'Edit these parameters accordingly' section.  This includes stating what the species transcriptome is, where the sequencing files are located, and where the kallisto index is located.

To run:
```
module load kallisto

*Edit the kallisto_all.sh script*

bash kallisto_all.sh
```
