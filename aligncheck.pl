##########################################################
#
#   aligncheck.pl
#   Summer, 2012-2013
#   Jason Bragg
#   Takes 1. a list of taxon names 
#         2. a directory (list) of sequence alignments
#   Writes a list of exons satisfying criteria
#
###########################################################

use warnings;
use strict;
use Bio::SimpleAlign;
use Bio::AlignIO; 
use Bio::LocatableSeq;


# directory for output files
my $outdir = "/home/jgb/skinktrans/exsaiph/";

# directory containing alignments (muscle)
my $alndir = "/home/jgb/skinktrans/out/";

# files with list of taxon names, exons
my $taxfil = "/home/jgb/skinktrans/tax_exsaiph.txt";
my $exfil  = "/home/jgb/skinktrans/allloci.txt";

# open files
open TAX, "<$taxfil" or die "cannot open list of taxa";
open EX,  "<$exfil" or die "cannot open list of aligned exons";

# put into arrays
my @tax = <TAX>;
my @ex  = <EX>;
my $numtax = scalar(@tax);

# initialize alignment hashes
my %tmpaln;
my %cataln;
my $catlst;

# process arrays a little...
chomp(@tax);
chomp(@ex);
@ex = grep(/\S/, @ex);
%cataln = map { $_ => "" } @tax;
$catlst = "";


# initialize the outfiles
my $exfiltfil = $outdir . "alnfilt.txt";
my $exinfofil = $outdir . "alninfo.txt";

open EXFILT, ">$exfiltfil" or die "cannot open list of taxa";
open EXINFO, ">$exinfofil" or die "cannot open list of aligned exons";
print EXINFO "exon" ."\t". "keep" ."\t". "alnlen" ."\t". "pcid" . "\t". "var" ."\t". "pic" ."\t". "len_ng" ."\t". "pcid_ng" . "\t". "var_ng" ."\t". "pic_ng" ."\t". "outgaps"  ."\t". "ingaps" . "\n";

# set some parameter values
my $lenthresh_ng  = 100;  # minimum length of sequences (leading lagging gaps removed)
my $pcidthresh_ng = 20;   # minimum percentage identity across alignment
my $ingpthresh_ng = 20;    # maximum percentage of internal gaps in any sequence
my $chngthresh_ng = 60;   # maximum percentage change in any single aligned sequence


######################################################
# loop through alignments, open, 
# test, list if satisfactory
######################################################

my %alltax;
my %exkeep;

foreach my $exon (@ex) {

    print "up to locus: $exon ...\n";
    
    $alltax{$exon} = "yes";
    $exkeep{$exon} = "yes";

    %tmpaln = map { $_ => "" } @tax;
    chomp($exon);

    my $alnfil = $alndir . $exon . ".fa.muscle";

    # open the alignment, store in $aln 
    my $seqio = Bio::AlignIO->new(
                             -file   => $alnfil,
                             -format => 'FASTA',
                             );
    my $aln = $seqio->next_aln();


    # loop through seqs in the alignment
    # catch seqs for nominated taxa in 'tmpaln'
    foreach my $tax (@tax) {
        chomp($tax);
        my $seq = $aln->get_seq_by_id($tax);

        if (defined $seq) {
             my $s = $seq->seq();
             my $i = $seq->id();
             $tmpaln{$i} = $s;   
             #dbg print "collect: " . $i . "\t" . $tmpaln{$i} . "\n";
        }
        else {  
            $alltax{$exon} = "no";
            print "Sequence missing, exon: $exon, taxon: for $tax\n";
            }
    }


   
    # make two alignment objects
    # with gaps, without gaps
    # only include alignments with all taxa

    if ($alltax{$exon} eq "yes") {

         my $exaln = Bio::SimpleAlign->new();
         my $exaln_ng;
 	 my $exaln_nog;
         foreach my $tax (@tax) {
                 
             my ($taxabbr, $dummy) = split(/_/,$tax,2);
             my $lseq = new Bio::LocatableSeq(-seq => $tmpaln{$tax},
     		                              -id  => $tax );
             $exaln->add_seq(-SEQ=>$lseq);
         }

         # Make alignment object without gaps
         my $gapline = $exaln->gap_line();

         if ($gapline =~ m/\./) {                                # only continue if alignment contains non-gaps
         my $firstnongap = index($gapline, ".") + 1;    
         my $lastnongap  = rindex($gapline, ".") + 1;

         $exaln_nog = $exaln->slice($firstnongap,$lastnongap);   # keep internal gaps
         $exaln_ng  = $exaln->remove_columns(['gaps']);          # remove all gaps

         my $outgaps = $exaln->length()  - $exaln_nog->length();
         my $ingaps  = $exaln_nog->length() - $exaln_ng->length();

         $exaln_ng = $exaln_ng->remove_columns(['weak']);       # remove ambigs

         my $phyfil = $outdir . $exon . ".phy"; 
         my $phyout = Bio::AlignIO->new(-file     => ">$phyfil",
                                        -format   => 'phylip',
                                        -idlength => 30);
         $phyout->write_aln($exaln);

         my $phyfil_ng = $outdir . $exon . "_ng.phy"; 
         my $phyout_ng = Bio::AlignIO->new(-file     => ">$phyfil_ng",
                                           -format   => 'phylip',
                                           -idlength => 30);
         $phyout_ng->write_aln($exaln_ng);
 

##############################################################################
## infoalign code -- remove
#
#         my $infoalnfil = $phyfil . ".ia.out";
#         my $infoalnfil_ng = $phyfil_ng . ".ia.out";
#
#         system("infoalign $phyfil $infoalnfil");
#         system("infoalign $phyfil_ng $infoalnfil_ng");
#
#        ### collect information
#         open ALNINFO, "<$infoalnfil_ng" or die "cannot open list of taxa";
#         my @info = <ALNINFO>;
#         my %taxchg = map { $_ => 0 } @tax;
#         my %taxigp = map { $_ => 0 } @tax;
#
#         foreach my $inflin (@info) {
#             chomp($inflin);
#             if (substr($inflin,0,1) ne "#")
#             {
#                    my @infbits = split(/\t/,$inflin);
#                    $taxchg{$infbits[1]} = $infbits[9];
#                    $taxigp{$infbits[1]} = $infbits[5];
#             }
#
#         }
#############################################################################

         my $alnpcid    = $exaln->percentage_identity;
         my $alnpcid_ng = $exaln_ng->percentage_identity;
         my $alnlen     = $exaln->length();
         my $alnlen_ng  = $exaln_ng->length();
         my ($alnpic, $alnvar) = calcVariable($exaln, $alnlen);
         my ($alnpic_ng, $alnvar_ng) = calcVariable($exaln_ng, $alnlen_ng);

 #print "$alnpic\n";
 #print "$alnvar\n";
 #print "$alnpic_ng\n";
 #print "$alnvar_ng\n";

         if ($alnlen_ng < $lenthresh_ng)               {$exkeep{$exon} = "no";}
         if ($alnpcid_ng < 70 || $alnpcid_ng == 100)   {$exkeep{$exon} = "no";}

 #        foreach my $t (@tax){
 #          if ($taxchg{$t} > $chngthresh_ng)                       {$exkeep{$exon} = "no";}
 #          if ( ($taxigp{$t}/$alnlen_ng * 100) > $chngthresh_ng)   {$exkeep{$exon} = "no";}
 #        }

         if ($exkeep{$exon} eq "yes")
         {
              print EXFILT $exon . "\n";
         }

         print EXINFO $exon ."\t". 1 ."\t". $alnlen ."\t". $alnpcid . "\t". $alnvar ."\t". $alnpic ."\t". $alnlen_ng ."\t". $alnpcid_ng . "\t". $alnvar_ng ."\t". $alnpic_ng ."\t". $outgaps  ."\t". $ingaps . "\n";
         }

         else {   # if the alignment is all gaps
         print EXINFO $exon ."\t". 0 ."\t\t\t\t\t\t\t\t\t\t\n";

         }
     }

    else { # if alltax($exon) ne "yes"
         print EXINFO $exon ."\t". 0 ."\t\t\t\t\t\t\t\t\t\t\n"; 

    }

}

close TAX;
close EX;
close EXINFO;
close EXFILT;




sub calcVariable {
         my $exaln = $_[0];
         my $alnlen = $_[1];
         my @bases = ("A", "C", "G", "T");
         my %count;
         my $alnpic = 0;
         my $alnvar = 0;
         for (my $pos=1; $pos <= $alnlen; $pos++)
         {
            %count = map { $_ => 0 } @bases;
            foreach my $tmps ($exaln->each_seq) {
               my $base = $tmps->subseq($pos, $pos);
               $count{$base}++;
            }
             
            my $charpi = 0;
            my $charva = 0;
            if ($count{"A"}+$count{"C"}+$count{"G"}+$count{"T"} == $numtax) # note make <= if permitting gaps non-ACGT sites
            {
                 if ($count{"A"} > 1) { $charpi = $charpi+1; }
                 if ($count{"C"} > 1) { $charpi = $charpi+1; }
                 if ($count{"G"} > 1) { $charpi = $charpi+1; }
                 if ($count{"T"} > 1) { $charpi = $charpi+1; }
                 if ($charpi > 1)   { $alnpic = $alnpic+1; }

                 if ($count{"A"} > 0) { $charva = $charva+1; }
                 if ($count{"C"} > 0) { $charva = $charva+1; }
                 if ($count{"G"} > 0) { $charva = $charva+1; }
                 if ($count{"T"} > 0) { $charva = $charva+1; }
                 if ($charva > 1)   { $alnvar = $alnvar+1; }

            }

         }

 return ($alnpic, $alnvar);

}



__END__



