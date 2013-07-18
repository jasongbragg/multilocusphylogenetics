##########################################################
#
#   Summer, 2012-2013
#   Jason Bragg
#
###########################################################

use warnings;
use strict;

my $dir = "/home/jgb/anchtag/eugongylusDNA/prelim/Alignments/";
my $phymldir = "/home/jgb/phymlbs/";

# files with list of taxon names, exons
my $mraicsumm = $dir . "mraic.alnfilt.out";
open MODSUM,  "<$mraicsumm" or die "cannot open summary of MrAIC output";
my @exmod  = <MODSUM>;
chomp(@exmod);

my $numrep = 1000;

foreach my $ex (@exmod)
{

  my ($loc, $model, $tree) = split(/\t/,$ex);
  my $modargs     = getphymlargs($model, $numrep);
  my $phyfil      = $phymldir . $loc . "_ng.phy";
  my $phymlparamfil = $dir . $loc . "_ng.phy.phymlbsparams";
  my $pfilcmd = $phyfil . "\n" . $modargs ; 

  open  PHYPAR, ">$phymlparamfil" or die "could not make a new parameter file";
  print PHYPAR $pfilcmd;
  close PHYPAR; 

  # system("phyml < $phymlparamfil");
}

close MODSUM;



        
sub getphymlargs {

my ($model, $numrep) = ($_[0], $_[1]) ;
my $params;
my $bsparams = "+\n+\nB\n" . $numrep ."\nY\nY\n";
if ($model eq "JC69") 	{	$params	=	"+\nM\nM\nM\nM\nM\nR\n".$bsparams."Y\n";	}		# JC69
if ($model eq "JC69I") 	{	$params	=	"+\nM\nM\nM\nM\nM\nV\nY\nR\n".$bsparams."Y\n";	}		# JC69I
if ($model eq "JC69G") 	{	$params	=	"+\nM\nM\nM\nM\nM\n".$bsparams."Y\n";		}		# JC69G
if ($model eq "JC69IG") {	$params	=	"+\nM\nM\nM\nM\nM\nV\nY\n".$bsparams."Y\n";	}		# JC69IG
if ($model eq "F81") 	{	$params	=	"+\nM\nM\nM\nM\nM\nM\nM\nF\nR\n".$bsparams."Y\n";}		# F81
if ($model eq "F81I") 	{	$params	=	"+\nM\nM\nM\nM\nM\nM\nM\nF\nV\nY\nR\n".$bsparams."Y\n";}	# F81I
if ($model eq "F81G") 	{	$params	=	"+\nM\nM\nM\nM\nM\nM\nM\nF\n".$bsparams."Y\n";	}		# F81G
if ($model eq "F81IG") 	{	$params	=	"+\nM\nM\nM\nM\nM\nM\nM\nF\nV\nY\n".$bsparams."Y\n";}	# F81IG
if ($model eq "K2P") 	{	$params	=	"+\nM\nM\nM\nM\nM\nM\nT\nY\nR\n".$bsparams."Y\n";}		# K2P
if ($model eq "K2PI") 	{	$params	=	"+\nM\nM\nM\nM\nM\nM\nT\nY\nR\nV\nY\n".$bsparams."Y\n";}	# K2PI
if ($model eq "K2PG") 	{	$params	=	"+\nM\nM\nM\nM\nM\nM\nT\nY\n".$bsparams."Y\n";}		# K2PG
if ($model eq "K2PIG") 	{	$params	=	"+\nM\nM\nM\nM\nM\nM\nT\nY\nV\nY\n".$bsparams."Y\n";}	# K2PIG
if ($model eq "HKY") 	{	$params	=	"+\nF\nT\nY\nR\n".$bsparams."Y\n";		}		# HKY
if ($model eq "HKYI") 	{	$params	=	"+\nF\nT\nY\nR\nV\nY\n".$bsparams."Y\n";	}		# HKYI
if ($model eq "HKYG") 	{	$params	=	"+\nF\nT\nY\n".$bsparams."Y\n";			}		# HKYG
if ($model eq "HKYIG") 	{	$params	=	"+\nF\nT\nY\nV\nY\n".$bsparams."Y\n";		}		# HKYIG
if ($model eq "SYM") 	{	$params	=	"+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n012345\n1.00\n1.00\n1.00\n1.00\n1.00\n1.00\nR\n".$bsparams."Y\n";}	# SYM
if ($model eq "SYMI") 	{	$params	=	"+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n012345\n1.00\n1.00\n1.00\n1.00\n1.00\n1.00\nR\nV\nY\n".$bsparams."Y\n";}	# SYMI
if ($model eq "SYMG") 	{	$params	=	"+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n012345\n1.00\n1.00\n1.00\n1.00\n1.00\n1.00\n".$bsparams."Y\n";}	# SYMG
if ($model eq "SYMIG") 	{	$params	=	"+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n012345\n1.00\n1.00\n1.00\n1.00\n1.00\n1.00\nV\nY\n".$bsparams."Y\n";}	# SYMIG
if ($model eq "GTR") 	{	$params	=	"+\nM\nM\nM\nF\nR\n".$bsparams."Y\n";		}		# GTR
if ($model eq "GTRI") 	{	$params	=	"+\nM\nM\nM\nF\nR\nV\nY\n".$bsparams."Y\n";	}		# GTRI
if ($model eq "GTRG") 	{	$params	=	"+\nM\nM\nM\nF\n".$bsparams."Y\n";		}		# GTRG
if ($model eq "GTRIG") 	{	$params	=	"+\nM\nM\nM\nF\nV\nY\n".$bsparams."Y\n";	}		# GTRIG

return $params;
}
