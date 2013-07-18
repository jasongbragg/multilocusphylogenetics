##########################################################
#
#   runMrAIC.pl
#   Summer, 2012-2013
#   Jason Bragg
#   Takes 1. a list of taxon names 
#         2. a directory (list) of sequence alignments
#   Writes a list of exons satisfying criteria
#
###########################################################

use warnings;
use strict;

my $phydir = "/home/jgb/anchtag/eugongylusDNA/prelim/Alignments/";
my $mraic  = "/home/jgb/software/mraic/MrAIC144/mraic.pl";

# files with list of taxon names, exons
my $exfil  = $phydir . "ex.alnfilt.txt";
open EX,  "<$exfil" or die "cannot open list of aligned exons";
my @exphy  = <EX>;
chomp(@exphy);


foreach my $ex (@exphy)
{
    my $exphyfil = $phydir . $ex . "_ng.phy";
    system( "perl $mraic $exphyfil"  );
}


