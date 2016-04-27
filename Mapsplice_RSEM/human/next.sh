#!/bin/bash
#shell script that is called by nextOuter.sh
#this script processes aligned reads in bam file, runs data thru RSEM, and further processes data
#NOTE: THIS SCRIPT IS NOT MEANT TO BE RUN ALONE
#for Sara's version of the script (with sample names rather than normal MapSplice output) look in FL-HCC-LR directory of RNASeq

for ARG in "$@"
do
DIR=`dirname $ARG`
BASE=`basename $ARG`    # _alignments.bam  ##add if bam file is: name_alignments.bam
echo "$BASE"   
cd $DIR
mkdir working

echo "samtools sort"
samtools sort $BASE sorted_genome_alignments
echo "samtools flagstat"
samtools flagstat sorted_genome_alignments.bam > sorted_genome_alignments.bam.flagstat
echo "samtools index"
samtools index sorted_genome_alignments.bam
echo "make bedGraph"
genomeCoverageBed -split -bg -ibam sorted_genome_alignments.bam -g /proj/seth_lab/projects/genome/hg19_TCGA  > $BASE.bedGraph
echo "make bigWig"
/proj/seth_lab/projects/bin/bedGraphToBigWig $BASE.bedGraph /proj/seth_lab/projects/genome/hg19.chromSizes $BASE.bW

echo "sort by chromosome, then read id"
perl /proj/seth_lab/projects/bin/ubu/src/perl/sort_bam_by_reference_and_name.pl --input sorted_genome_alignments.bam --output working/sorted_by_chr_read.bam --temp-dir . --samtools /nas02/apps/samtools-0.1.19/bin/samtools > working/sorted_by_chr_read.log 2> working/sorted_by_chr_read.log

echo "translate to transcriptome coordinates"
java -Xms3G -Xmx3G -jar /proj/.test/roach/miRNA/bin/ubu.jar sam-xlate --bed /proj/seth_lab/projects/genome/hg19_TCGA/unc_hg19.bed --in working/sorted_by_chr_read.bam --out working/transcriptome_alignments.bam --order /proj/seth_lab/projects/genome/hg19_TCGA/hg19_M_rCRS_ref.transcripts.fa --xgtags --reverse > working/genome_to_transcriptome.log 2> working/genome_to_transcriptome.log

echo "filter indels, large inserts, zero mapping quality from transcriptome bam"
java -Xmx512M -jar /proj/.test/roach/miRNA/bin/ubu.jar sam-filter --in working/transcriptome_alignments.bam --out working/transcriptome_alignments_filtered.bam --strip-indels --max-insert 10000 --mapq 1 > working/sam_filter.log 2> working/sam_filter.log

echo "RSEM"
/proj/seth_lab/projects/bin/RSEM/rsem-calculate-expression --paired-end --bam --estimate-rspd -p 1 working/transcriptome_alignments_filtered.bam /proj/seth_lab/projects/genome/hg19_TCGA/hg19_M_rCRS_ref rsem > working/rsem.log 2> working/rsem.log

echo "strip trailing tabs from rsem.isoforms.results"
perl /proj/seth_lab/projects/bin/ubu/src/perl/strip_trailing_tabs.pl --input rsem.isoforms.results --temp working/orig.isoforms.results > working/trim_isoform_tabs.log 2> working/trim_isoform_tabs.log

echo "prune isoforms from gene quant file"
mv rsem.genes.results working/orig.genes.results; sed /^uc0/d working/orig.genes.results > rsem.genes.results

#echo "normalize gene quant"
#perl /proj/seth_lab/projects/bin/ubu/src/perl/quartile_norm.pl -c 5 -q 75 -t 1000 -o ${BASE}_rsem.genes.normalized_results -skip 1 -also 1 -also 2 -also 3 -also 4 -also 5 -also 6 ${BASE}_rsem.genes.results

#echo "normalize isoform quant"
#perl /proj/seth_lab/projects/bin/ubu/src/perl/quartile_norm.pl -c 7 -q 75 -t 300 -o ${BASE}_rsem.isoforms.normalized_results -skip 1 -also 1 -also 2 -also 3 -also 4 -also 5 -also 6 -also 8 rsem.isoforms.results

#echo "junction counts"
#java -Xmx512M -jar /proj/.test/roach/miRNA/bin/ubu.jar sam-junc --junctions splice-junctions.txt --in sorted_genome_alignments.bam --out junction_quantification.txt > junction_quantification.log 2> junction_quantification.log

#echo "exon counts"
#coverageBed -split -abam sorted_genome_alignments.bam -b composite_exons.bed | perl normalizeBedToolsExonQuant.pl composite_exons.bed > bt.exon_quantification.txt 2> bt_exon_quantification.log

#echo "cleanup large intermediate output"
#rm alignments.bam logs/* working/sorted_by_chr_read.bam working/transcriptome_alignments.bam working/transcriptome_alignments_filtered.bam > working/cleanup.log

#do
#
#LogNM=$DIR.next.log
#bsub -o $LogNM sh next.sh $ARG
done
