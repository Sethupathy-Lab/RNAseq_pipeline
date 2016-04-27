#!/usr/bin/perl -w

sub commify {
   local $_  = shift;
   1 while s/^([-+]?\d+)(\d{3})/$1,$2/;
   return $_;
}


$header='';
foreach $file (@ARGV) {

   open (FILE,$file) or die();
   while (<FILE>) {
      chomp();
      unless (($_ =~ /\#/)||(length($_)==0)) {
	 @parts = split(/\t/,$_);
	 if (scalar(@parts)>=2) {
	    $key = join(':',$header,$parts[0]);
	    $data{$file}{$key} = $parts[1];
	    if (scalar(@parts==3)) {
	       $perc{$file}{$key}=$parts[2];
	       $perKeys{$key}=1;
	    }
	    if (!exists($keylist{$key})) {
	       push(@keyArr,$key);
	       $keylist{$key} =1;
	    }

	 }
	 else{
	    $header=$_;
	 }
      }
   }
   close(FILE);
}

print "#MapSplice Summary\n<table border=1>\n<tr><td></td>\n";
foreach $f (@ARGV) {
   @fparts = split("_",$f);
   print "<td>$fparts[0]</td>";
}
print "\n</tr><tr><td>Mate pair quality files</td>";
foreach $f (@ARGV) {
   @fparts = split("/",$f);
   $p1="Qual/".$fparts[0]."_R1_001_fastqc.html";
   $p2="Qual/".$fparts[0]."_R2_001_fastqc.html";
   print "<td><a href=$p1> 1 </a>\n";
   print "&nbsp;&nbsp;<a href=$p2> 2 </a>\n</td>";

}
print "\n</tr>";
foreach $k(@keyArr){
   print "<tr><td>$k</td>";
   foreach $f (@ARGV) {
#<td style="mso-number-format:"mm\/dd\/yyyy"> 10/01/2011 </td>

      $formatted = commify($data{$f}{$k});
      print "<td>$formatted</td>";
   }
   if (exists($perKeys{$k})){
      print "</tr><tr>\n<td></td>";
      foreach $f (@ARGV) {
	 print "<td>$perc{$f}{$k}</td>";
      }
   }
   print "</tr>\n";
}

