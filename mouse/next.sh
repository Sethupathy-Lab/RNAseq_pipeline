#!/bin/bash

for ARG in "$@"
do
DIR=`dirname $ARG`
cd $DIR
mkdir working
samtools sort alignments.bam sorted_genome_alignments
samtools flagstat sorted_genome_alignments.bam > sorted_genome_alignments.bam.flagstat
samtools index sorted_genome_alignments.bam
genomeCoverageBed -split -bg -ibam sorted_genome_alignments.bam -g /proj/seth_lab/projects/genome/mm10  > $DIR.bedGraph
/proj/seth_lab/projects/bin/bedGraphToBigWig $DIR.bedGraph /proj/seth_lab/projects/genome/mm10.chromSizes $DIR.bW

perl /proj/seth_lab/projects/bin/ubu/src/perl/sort_bam_by_reference_and_name.pl --input sorted_genome_alignments.bam --output working/sorted_by_chr_read.bam --temp-dir . --samtools /nas02/apps/samtools-0.1.19/bin/samtools > working/sorted_by_chr_read.log 2> working/sorted_by_chr_read.log

java -Xms3G -Xmx3G -jar /proj/.test/roach/miRNA/bin/ubu.jar sam-xlate --bed /proj/seth_lab/projects/genome/mm10/mm10.ucsc.bed --in working/sorted_by_chr_read.bam --out working/transcriptome_alignments.bam --order /proj/seth_lab/projects/genome/mm10/mouse_mm10_125.transcripts.fa --xgtags --reverse > working/genome_to_transcriptome.log 2> working/genome_to_transcriptome.log

java -Xmx512M -jar /proj/.test/roach/miRNA/bin/ubu.jar sam-filter --in working/transcriptome_alignments.bam --out working/transcriptome_alignments_filtered.bam --strip-indels --max-insert 10000 --mapq 1 > working/sam_filter.log 2> working/sam_filter.log

/proj/seth_lab/projects/bin/RSEM/rsem-calculate-expression --bam --estimate-rspd -p 1 working/transcriptome_alignments_filtered.bam /proj/seth_lab/projects/genome/mm10/mouse_mm10_125 rsem > working/rsem.log 2> working/rsem.log

perl /proj/seth_lab/projects/bin/ubu/src/perl/strip_trailing_tabs.pl --input rsem.isoforms.results --temp working/orig.isoforms.results > working/trim_isoform_tabs.log 2> working/trim_isoform_tabs.log

mv rsem.genes.results working/orig.genes.results; sed /^uc0/d working/orig.genes.results > rsem.genes.results

perl /proj/seth_lab/projects/bin/ubu/src/perl/quartile_norm.pl -c 5 -q 75 -t 1000 -o rsem.genes.normalized_results -skip 1 -also 1 -also 2 -also 3 -also 4 -also 6 -also 7 rsem.genes.results

perl /proj/seth_lab/projects/bin/ubu/src/perl/quartile_norm.pl -c 5 -q 75 -t 300 -o rsem.isoforms.normalized_results -skip 1 -also 1 -also 2 -also 3 -also 4 -also 6 -also 7 -also 8 rsem.isoforms.results

#java -Xmx512M -jar /proj/.test/roach/miRNA/bin/ubu.jar sam-junc --junctions splice-junctions.txt --in sorted_genome_alignments.bam --out junction_quantification.txt > junction_quantification.log 2> junction_quantification.log

#coverageBed -split -abam sorted_genome_alignments.bam -b composite_exons.bed | perl normalizeBedToolsExonQuant.pl composite_exons.bed > bt.exon_quantification.txt 2> bt_exon_quantification.log

#rm alignments.bam logs/* working/rg_alignments.bam working/sorted_by_chr_read.bam working/transcriptome_alignments.bam working/transcriptome_alignments_filtered.bam working/prep_1.fastq working/prep_2.fastq > working/cleanup.log

#do
#
#LogNM=$DIR.next.log
#bsub -o $LogNM sh next.sh $ARG
done
