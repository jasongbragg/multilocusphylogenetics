##########################################################
#
#  phy2concat.pl
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


# directory containing alignments (muscle)
my $alndir = "/home/jgb/anchtag/eugongylusDNA/prelim/Alignments/";

# files with list of taxon names, exons
my $taxfil = "taxabbr.txt";
my $exfil  = "ex.alnfilt.txt";


# open files
open TAX, "<$taxfil" or die "cannot open list of taxa";
open EX,  "<$exfil" or die "cannot open list of aligned exons";

# put into arrays
my @tax = <TAX>;
my @ex  = <EX>;

# initialize alignment hashes
my %tmpaln;
my %cataln;

my $tmpfirst;
my $tmplast;
my $catlst;
my $catlen = 0;

chomp(@tax);
%cataln = map { $_ => "" } @tax;
$catlst = "";

# open files
my $concatoutfil = "concat/concatenated.out";
open CATOUT, ">$concatoutfil" or die "cannot open list of taxa";


######################################################
# loop through the list of exons
# open phylip file
#
# get sequences for each nominated taxon
#         -- add sequence to concatenated alignment
#         -- sample SNPs for SNAPP
######################################################

foreach my $exon (@ex) {

    %tmpaln = map { $_ => "" } @tax;
    chomp($exon);

    my $alnfil = $alndir . $exon . "_ng.phy";

    # open the alignment, store in $aln 
    my $seqio = Bio::AlignIO->new(
                             -file   => $alnfil,
                             -format => 'phylip',
                             );
    my $aln = $seqio->next_aln();

    # loop through seqs in the alignment
    # catch seqs for nominated taxa in tmpaln
    foreach my $tax (@tax) {
        chomp($tax);
        my $seq = $aln->get_seq_by_id($tax);

        if (defined $seq) {
             my $s = $seq->seq();
             my $i = $seq->id();
             $tmpaln{$i} = $s;   
             #dbg print "collect: " . $i . "\t" . $tmpaln{$i} . "\n";
        }
    }

    # check that nothing crazy has happened
    my $isflush = $aln->is_flush;
    my $alen    = $aln->length();

    # if satisfied, do concatenation
    if ($catlen > 0  && $isflush == 1){
              foreach my $t (@tax) {
              $cataln{$t} = $cataln{$t} . $tmpaln{$t};
            #dbg print $cataln{$t};
              }
        $tmpfirst = $catlen+1;
        $tmplast  = $catlen+$alen;
        print CATOUT $exon ."\t". $tmpfirst ."\t". $tmplast ."\n";
        $catlen=$catlen+$alen;
    }

    if ($catlen == 0 && $isflush == 1){
              foreach my $t (@tax) {
              $cataln{$t} = $tmpaln{$t};
              }

        $tmpfirst = 1;
        $tmplast  = $alen;
        print CATOUT $exon ."\t". $tmpfirst ."\t". $tmplast ."\n";
        $catlen = $alen;
    }
}


### make concat alignments

my $caln = Bio::SimpleAlign->new();
foreach my $tax (@tax) {
   #print $cataln{$tax};
   my $cseq = new Bio::LocatableSeq(-seq => $cataln{$tax},
     				    -id  => $tax);
   $caln->add_seq(-SEQ=>$cseq);
}

my $catphyfil = "concatenated.phy"; 
my $phycat    = Bio::AlignIO->new(-file     => ">$catphyfil",
                                  -format   => 'phylip',
                                  -idlength => 30);
$phycat->write_aln($caln);

my $catnexfil = "concatenated.nex";
my $nexcat = Bio::AlignIO->new(-file     => ">$catnexfil",
                                   -format   => 'nexus',
                                   -show_symbols => 0,
                                   -show_endblock => 0,
                                   -idlength => 30);

$nexcat->write_aln($caln); 

   
close TAX;
close EX;
close CATOUT;
