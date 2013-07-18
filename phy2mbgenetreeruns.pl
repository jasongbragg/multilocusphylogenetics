##########################################################
#
#  phy2mbgenetreeruns.pl
#  Summer, 2012-2013
#  Jason Bragg
#  Takes 1. a list of sequence names, 
#        2. a list of taxon names, 
#        3. a directory of sequence alignments
#  writes concatenated alignment files 
#
###########################################################

use warnings;
use strict;
use Bio::SimpleAlign;
use Bio::AlignIO; 
use Bio::LocatableSeq;


# nominate directory 
my $dir = "/home/jgb/anchtag/eugongylusDNA/prelim/Alignments/";
my $mbgtdir = "/home/jgb/anchtag/eugongylusDNA/prelim/Alignments/mbgenetrees/";
my $rundir = "/home/jgb/mbgenetrees/";

# files with list of taxon names, exons
my $exfil = $dir . "ex.alnfilt.txt";

# open files
open EX,  "<$exfil" or die "cannot open list of aligned exons";

# put into arrays
my @ex  = <EX>;

foreach my $exon (@ex) {

  chomp($exon);
  my $mraicfil = $dir . $exon . "_ng.phy.MrAIC.txt";
  my $phyfil   = $dir . $exon . "_ng.phy";

  my ($nst, $rates, $statefreqpr) = getbestmodel($mraicfil);

  my $expfil    = $mbgtdir . $exon . "_ng_params.txt";
  my $nexfil    = $mbgtdir . $exon . "_ng.nex";
  my $nexfil_nm = $rundir . $exon . "_ng.nex";

  makenexfil($phyfil,$nexfil);
  makembparfil($nexfil_nm, $expfil, $nst, $rates, $statefreqpr);

#  system("mb <$expfil> $logfil");
}


######################################################################################
### Subroutines

### extract the best model from MrAIC outfile
sub getbestmodel {
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

   ($nst, $rates, $statefreqpr) = modelmap($model);
   return ($nst, $rates, $statefreqpr);

}


### make a parameter file for mrbayes
sub makembparfil {

   my ($nexfil_nm, $expfil, $nst, $rates, $statefreqpr) = ($_[0], $_[1], $_[2], $_[3], $_[4]);
   open EXP, ">$expfil" or die "cannot open exon file";

   print EXP<<RUN;
   set autoclose=yes nowarn=yes 
   exe $nexfil_nm
   lset nst=$nst rates=$rates
   prset statefreqpr=$statefreqpr
   mcmc samplefreq=100 diagnfreq=10000 nrun=4 nchain=2 stopval=0.01 file=$nexfil_nm.out
   sump burnin=1000
   sumt burnin=1000
   quit
RUN

   close EXP;
}


### write a nex file, from the phy file

sub makenexfil {

    my ($phyfil, $nexfil) = ($_[0], $_[1]);

    my $seqio = Bio::AlignIO->new(
                             -file   => $phyfil,
                             -format => 'phylip',
                             );
    my $aln = $seqio->next_aln();


    my $nexout = Bio::AlignIO->new(-file     => ">$nexfil",
                                   -format   => 'nexus',
                                   -show_symbols => 0,
                                   -show_endblock => 0,
                                   -idlength => 30);

    $nexout->write_aln($aln); 

}




### brute force mapping of MrAIC model
### to mrbayes input parameters

sub modelmap {
my $model = $_[0]; 
my $nst; my $rates; my $statefreqpr;

# HKY
if ($model eq "HKY")    {$nst=2; $rates="equal";    $statefreqpr="Dirichlet(1,1,1,1)";}
if ($model eq "HKYI")   {$nst=2; $rates="propinv";  $statefreqpr="Dirichlet(1,1,1,1)";}
if ($model eq "HKYG")   {$nst=2; $rates="gamma";    $statefreqpr="Dirichlet(1,1,1,1)";}
if ($model eq "HKYIG")  {$nst=2; $rates="invgamma"; $statefreqpr="Dirichlet(1,1,1,1)";}

# K2P
if ($model eq "K2P")    {$nst=2; $rates="equal";    $statefreqpr="fixed(equal)";}
if ($model eq "K2PI")   {$nst=2; $rates="propinv";  $statefreqpr="fixed(equal)";}
if ($model eq "K2PG")   {$nst=2; $rates="gamma";    $statefreqpr="fixed(equal)";}
if ($model eq "K2PIG")  {$nst=2; $rates="invgamma"; $statefreqpr="fixed(equal)";}

# F81
if ($model eq "F81")    {$nst=1; $rates="equal";    $statefreqpr="Dirichlet(1,1,1,1)";}
if ($model eq "F81I")   {$nst=1; $rates="propinv";  $statefreqpr="Dirichlet(1,1,1,1)";}
if ($model eq "F81G")   {$nst=1; $rates="gamma";    $statefreqpr="Dirichlet(1,1,1,1)";}
if ($model eq "F81IG")  {$nst=1; $rates="invgamma"; $statefreqpr="Dirichlet(1,1,1,1)";}

# JC69
if ($model eq "JC69")   {$nst=1; $rates="equal";    $statefreqpr="fixed(equal)";}
if ($model eq "JC69I")  {$nst=1; $rates="propinv";  $statefreqpr="fixed(equal)";}
if ($model eq "JC69G")  {$nst=1; $rates="gamma";    $statefreqpr="fixed(equal)";}
if ($model eq "JC69IG") {$nst=1; $rates="invgamma"; $statefreqpr="fixed(equal)";}

# GTR
if ($model eq "GTR")    {$nst=6; $rates="equal";    $statefreqpr="Dirichlet(1,1,1,1)";}
if ($model eq "GTRI")   {$nst=6; $rates="propinv";  $statefreqpr="Dirichlet(1,1,1,1)";}
if ($model eq "GTRG")   {$nst=6; $rates="gamma";    $statefreqpr="Dirichlet(1,1,1,1)";}
if ($model eq "GTRIG")  {$nst=6; $rates="invgamma"; $statefreqpr="Dirichlet(1,1,1,1)";}

# SYM
if ($model eq "SYM")    {$nst=6; $rates="equal";    $statefreqpr="fixed(equal)";}
if ($model eq "SYMI")   {$nst=6; $rates="propinv";  $statefreqpr="fixed(equal)";}
if ($model eq "SYMG")   {$nst=6; $rates="gamma";    $statefreqpr="fixed(equal)";}
if ($model eq "SYMIG")  {$nst=6; $rates="invgamma"; $statefreqpr="fixed(equal)";}

return ($nst, $rates, $statefreqpr);
}



