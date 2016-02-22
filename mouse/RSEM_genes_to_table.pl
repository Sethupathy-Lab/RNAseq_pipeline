#!/usr/bin/perl
#RSEM_genes_to_table.pl
#use this script to convert RSEM gene level output (rsem.genes.results) into table format
#usage: perl RSEM_genes_to_table.pl <all rsem.genes.results to be combined> > output.txt
#for all rsem.genes.results input files use */rsem.genes.results and run from parent directory
#not necessary to put all files within same parent directory. script modified to allow RSEM outputs to be combined from different experiments (aka sorted in different parent directories)
#inputs: all the rsem.genes.results files that will be combined into one table
#outputs: a table with all the sample names, gene names, and expected count for each gene in each sample


use warnings;
use File::Basename;

die "usage: <rsem.genes.results>" unless @ARGV > 0;
my @files = @ARGV;

my %counts;

print "Gene";

foreach my $file (@files) {
   open (IN, $file) or die "Can't open $file\n";
#   my @name = split('/', $file); #split file name based on '/', first thing split should be the directory of the input file which should specific the sample
   my ($filename, $directory) = fileparse($file);
   my @sampleID;
   my $position;
   if ($filename =~ m/_rsem/) { #if sample is named in prefix of filename, get sample name from there
      @sampleID = split ("_rsem", $filename);
      $position = 0;
   } else { #sample is named by parent directory
      @sampleID = split(/\//, $directory);
      $position = scalar(@sampleID) - 1;
   }
   print "\t$sampleID[$position]"; #leading tab accounts for gene names

   while (<IN>) {
      chomp;
      my @results = split(/\t/, $_);
      my $gene = $results[0];
      my $expected_count = $results[4];

      if ($gene eq "gene_id") { #makes sure header line is not read in
	 next;
      }

      if (exists $counts{$gene}) {
	 $counts{$gene} = $counts{$gene} . "\t" . $expected_count;
      }
      if (!exists $counts{$gene}) {
	 $counts{$gene} = $expected_count;
      }
   }
}

print "\n";
foreach my $gene_name (keys %counts) {
   print "$gene_name\t$counts{$gene_name}\n";
}
