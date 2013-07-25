
use strict;
use warnings;

my $phymlbsdir = "/home/jgb/anchtag/eugongylusDNA/prelim/Alignments/phymlbs/";
my @bstreefiles = <$phymlbsdir*phyml_boot_trees.txt>;
# print @bstreefiles;
my $out = $phymlbsdir . "genetreebstrees.out"; 
open OUT, ">$out" or die "can't open a treefile";

foreach my $t (@bstreefiles)
{


   open TRE, "<$t" or die "can't open a treefile";
   my @tre = <TRE>;
   chomp(@tre);
   my $c = 1;
   foreach my $b (@tre)
   {
       if ($c > 1) {print OUT "\t";}
       print OUT $b;
       $c=$c+1;
   }
   print OUT "\n";
   close TRE;
} 

close OUT;
