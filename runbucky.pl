use strict;
use warnings;

my $mbgtdir = "/home/jgb/anchtag/eugongylusDNA/prelim/Alignments/mbgenetrees/out/";
my $outstem = "/home/jgb/anchtag/eugongylusDNA/prelim/Alignments/mbgenetrees/out/buckyrun";


system("bucky $mbgtdir*.in");
