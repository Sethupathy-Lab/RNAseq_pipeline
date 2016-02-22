#!/usr/bin/perl -w 
open (REG,'register') or die();
while (<REG>) {
   chomp();
   my($A,$B) = split(/\t/,$_);
   $reg{$A}=$B;
}
close(REG);


foreach $f(@ARGV) {
   open (FILE,$f) or die ();
   $header=<FILE>;
   while (<FILE>) {
      chomp();
      @parts=split(/\t/,$_);
      @transcripts=split(",",$parts[0]);
      %transSymbs=();
      foreach $t(@transcripts) {
	 $transSymbs{$t}=1;
      }
      @vals=keys(%transSymbs);
      $symbList=join(",",@vals);
      $Count{$f}{$symbList} = $parts[4];
      $symbs{$symbList}=1;
   }
   close(FILE);
}

foreach $f (@ARGV) {
   @fparts=split("/",$f);
   print "\t$reg{$fparts[0]}";
}
print "\n";
foreach $f (@ARGV) {
   @fparts=split("/",$f);
   @nm=split("_",$fparts[0]);
   print "\t$nm[0]";
}
print "\n";

foreach $s(keys(%symbs)) {
   print "$s";
   foreach $f(@ARGV) {
      if (exists($Count{$f}{$s})){
      print "\t$Count{$f}{$s}" ;
      }
      else {
	 print "\tNA";
      }

   }
   print "\n";
}
