##########################################################
#
#   parseMrAIC.pl
#   Summer, 2012-2013
#   Jason Bragg
#   Takes 1. a list of exon names 
#         2. mraic output
#   Writes a list of models and trees
#
###########################################################

use warnings;
use strict;

my $dir = "/home/jgb/anchtag/eugongylusDNA/prelim/Alignments/";


# files with list of taxon names, exons
my $exfil  = $dir . "ex.alnfilt.txt";
open EX,  "<$exfil" or die "cannot open list of aligned exons";
my @exphy  = <EX>;
chomp(@exphy);

my $mraicsumm = $dir . "mraic.alnfilt.out";
open OUT,  ">$mraicsumm" or die "cannot open outfile";

foreach my $ex (@exphy)
{
  my $mraicfil = $dir . $ex . "_ng.phy.MrAIC.txt";
  my $model = getbestmodname($mraicfil);
  my $treefil = $dir . $ex ."_ng.phy.AICc-". $model .".tre";
  my $tree  = getbesttree($treefil);

  print OUT $ex ."\t". $model ."\t". $tree ."\n";
}

close EX;
close OUT;


sub getbestmodname {
   my $mraicfil = $_[0];   
   my ($nst, $rates, $statefreqpr) = "null";
   print "\n\n" . $mraicfil . "\n\n";
   open MAF,  "<$mraicfil" or die "cannot open the MrAIC file";


   my $model = "null";
   while (my $line = <MAF>)
   {
       chomp($line);
       if ($line =~ m/Minimum AICc model: /)
       {
           # print "\n\n" .$line ."\n\n" ;
           $model = $line;
           $model =~ s/Minimum AICc model: //;          
       
       }
   }
   close MAF;
  
   return $model;

}


sub getbesttree {
   my $treefil = $_[0]; 
   open TREE,  "<$treefil" or die "cannot open the tree file";
   my $tree = <TREE>;
   chomp($tree);
   return $tree;
}

