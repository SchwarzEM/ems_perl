#!/usr/bin/perl

# ali_splices2genes.pl -- Erich Schwarz <emsch@its.caltech.edu>, 2/19/2008.
# Purpose: given CDS->WBGene list, map Ali-Solexa splices (@ CDS names).

use strict;
use warnings;

my $chromosome = 'I|II|III|IV|V|X';

unless ($#ARGV == 1) { 
    die "Format: ./ali_splices2genes.pl [namefile] [splicefile]\n";
}

my $namefile   = $ARGV[0];
my $splicefile = $ARGV[1];

my %cds2wbgene = ();
my %wbgene_count = ();

open my $NAMES, "<", "$namefile" 
    or die "Can't open name file $namefile\n";

while ( my $input = <$NAMES> ) { 
    chomp $input;

    # Sample input line:
    # WBGene00000001  aap-1   Y110A7A.10      aap-1 encodes ...

    if ( $input =~ / ( WBGene\d+ )  # $1
                     .+? \s         # \s forces \w+ to lack '-'
                     ( \w+\.\w+ )   # $2; '\w' is to block '-'
                   /xms ) {
        my ($wbgene, $cds);
        ($wbgene, $cds) = ($1, $2);
        $cds =~ s/[a-z]\z//;
        $cds2wbgene{$cds} = $wbgene;
    }
}

close $NAMES;

open my $SPLICES, "<", "$splicefile" 
    or die "Can't open splice file $splicefile\n";

while ( my $input = <$SPLICES> ) { 
    chomp $input;

    # Sample input line:
    # chrIII  4753789 4753886 B0393.1-7F1 ...

    if ( $input =~ / chr(?: $chromosome )  # not $1 !
                     \s+ \d+ \s+ \d+ \s+
                     ( \w+\.\w+ )          # e.g. "B0393.1"
#                    \- \S+
                   /xms ) { 
        my $cds;
        $cds = $1;
        $cds =~ s/[a-z]\z//;
        if ( exists $cds2wbgene{$cds} ) { 
            my $wbgene = $cds2wbgene{$cds};
            $wbgene_count{$wbgene}++;
        }
        if (! exists $cds2wbgene{$cds} ) { 
            warn "Couldn't find WBGene name for CDS name $cds\n";
        }
    }
}

close $SPLICES;

my @genes = sort 
    { $wbgene_count{$b} <=> $wbgene_count{$a} } keys %wbgene_count;

foreach my $gene (@genes) { 
    print "$gene\t$wbgene_count{$gene}\n";
}

