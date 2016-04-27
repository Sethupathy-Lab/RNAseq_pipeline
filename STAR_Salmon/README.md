# STAR/Salmon RNAseq analyses 

Both STAR and Salmon are modules available on the cluster.  The environment for this pipeline can be set-up by sourcing RNAseq_env.sh from the cluster: '''source RNAseq_env.sh'''

STAR is a fast read aligner.  First, a STAR-indexed genome needs to be generated.  An example call is below:

bsub -n 8 -R "span[hosts=1]" STAR --runMode genomeGenerate \
--runThreadN 8 \
--genomeFastaFiles /proj/seth_lab/projects/genome/star_salmon_dev/mm9/mm9.fa \
--genomeDir /proj/seth_lab/projects/genome/star_salmon_dev/mm9 \
--sjdbGTFfile /proj/seth_lab/projects/genome/star_salmon_dev/mm9/mm9.gtf \
--sjdbOverhang 100


