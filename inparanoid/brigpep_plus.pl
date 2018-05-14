#!/usr/bin/perl

# genes2brigpep.pl -- Erich Schwarz <emsch@its.caltech.edu>, 3/19/2008.
# Purpose: add briggsae gene names to header lines of brigpep2 FASTA.

use strict;
use warnings;
use Ace;

unless ($#ARGV == 2) {
    die 'Usage [example]: ./genes2brigpep.pl brigpep2',
        ' 21633.C_briggsae.dat',
        ' /usr/local/acedb/elegans_WS180',
        ' > brigpep2gid',
        "\n",
        ;
}

my $brig_fasta    = $ARGV[0];
my $brig_dat      = $ARGV[1];
my $database      = $ARGV[2];

my $input          = q{};

my $bprot          = q{};
my %bprots2bdata   = ();

my $uniprot        = q{};
my %seen           = ();
my $cds            = q{};
my %cbgs2data      = ();
my %uniprots2cdses = ();
my %ncbi_ids2cdses = ();
my $ncbi_id        = q{};

open (my $BRIG_FASTA, "$brig_fasta") or die "Can't open $brig_fasta: $!";
while ($input = <$BRIG_FASTA>) {
    chomp $input;
    if ($input =~ /\A > (\S+) /x) {
        $bprot = $1;
    }
    elsif ( ($input !~ /\A >/x) and ($input =~ /[a-zA-Z]/x) ) {
        $input =~ tr/[^a-zA-Z]//;
        $input =~ tr/[a-z]/[A-Z]/;
        $bprots2bdata{$bprot}->{"seq"} .= $input;
    }
}
close $BRIG_FASTA;

open (my $BRIG_DAT, "$brig_dat") or die "Can't open $brig_dat: $!";
while ($input = <$BRIG_DAT>) { 
my $pattern1 = '\A GN \s+ Name=                    (CBG\d+)    ';
my $pattern2 = '\A GN .+  ORFNames=                (CBG\d+) ;  ';
my $pattern3 = '\A DR \s+ WormBase; \s+ WBGene\d+; (CBG\d+) \. ';
    if ($input =~ /\A AC \s+ (\S+) ; /x) { 
        $uniprot = "UNIPROT:" . $1;
        if ( $seen{$uniprot} ) { 
             die "Observed Uniprot acc. $uniprot twice\n";
        }
        $seen{$uniprot} = 1;
    }

    elsif ( $input =~ / integrated \s into \s UniProtKB \/ Swiss-Prot \. /x ) {
        $uniprot =~ s/UNIPROT:/SW:/;
    }
    elsif ( $input =~ / integrated \s into \s UniProtKB \/ TrEMBL \. /x ) {
        $uniprot =~ s/UNIPROT:/TR:/;
    }

    # N.B.: somewhat unsystematic storage of CBG\d+ data.
    # Most common form is "GN   Name=CBG".
    # Others -- "GN   ORFNames=" ; "DR   WormBase..." ; "ORFNames=CBG"!

    elsif (    ( $input =~ / $pattern1 /xms ) 
            or ( $input =~ / $pattern2 /xms ) 
            or ( $input =~ / $pattern3 /xms ) ) {
        $cds = $1;
        if ( ( $cbgs2data{$cds}->{"uniprot"} ) and 
             ( $uniprots2cdses{ $cbgs2data{$cds}->{"uniprot"} } ne $cds ) ) { 
            die "Mapped ", , " to ", $cbgs2data{$cds}->{"uniprot"}, " already\n";
        }
        $cbgs2data{$cds}->{"uniprot"} = $uniprot;
        $uniprots2cdses{$uniprot} = $cds;
    }

    elsif ($input =~ /\A DR \s+ (\w+;\s+)* (\w+?\d+\.\d) ; /x ) {
        $ncbi_id = "protein_id:" . $2;
        if ( ( $ncbi_ids2cdses{$ncbi_id} ) and ( $ncbi_ids2cdses{$ncbi_id} ne $cds ) ) { 
            die "Observed protein ID $ncbi_id twice\n";
        }
        $ncbi_ids2cdses{$ncbi_id} = $cds;
        $cbgs2data{$cds}->{"ncbi_id"} = $ncbi_id;
    }
}
close $BRIG_DAT;

# Igor Antoschechkin says (10/17/2006):
# AQL gets confused, so make variables unambiguous.

my $aql_query = 'select p'                        # $ref->[0]
                . ', c'                           # $ref->[1]
                . ', g'                           # $ref->[2]
                . ' from p in class Protein '
                . ', c in p->Corresponding_CDS'
                . ', g in c->Gene'
                . ' where p->Species like "Caenorhabditis briggsae"'
    ;

my $db = Ace->connect(-path=>$database) 
    or die "No local database: ",Ace->error;

my @objects = $db->aql($aql_query);

foreach my $ref (@objects) { 
    $bprot = $ref->[0];
    $cds  = $ref->[1];
    my $gene = $ref->[2];
    $bprot =~ s/\A BP://x;
    $bprots2bdata{$bprot}->{"cds"}   = $cds;
    $bprots2bdata{$bprot}->{"gene"}  = $gene;
}

foreach $bprot (sort keys %bprots2bdata) { 
    if ( $bprots2bdata{$bprot}->{"seq"} ) {
        my $gene = $bprots2bdata{$bprot}->{"gene"};
        $cds  = $bprots2bdata{$bprot}->{"cds"};

        print '>';
        print $bprot;

        print '    ';
        print $gene;
        print '; ';
        print $cds;

        print '    ';
        if ( $cbgs2data{$cds}->{"uniprot"} ) {
            print $cbgs2data{$cds}->{"uniprot"};
        }

        print '    ';
        if ( $cbgs2data{$cds}->{"ncbi_id"} ) {
            print $cbgs2data{$cds}->{"ncbi_id"};
        }

        print "\n";    

        my $fullseq = $bprots2bdata{$bprot}->{"seq"};
        my @seqs 
            = unpack("a60" x (length($fullseq)/60 + 1), $fullseq );
        foreach my $aa60 (@seqs) {
            print "$aa60\n" if ($aa60 =~ /\S/);
        }
    }
}

