#!/usr/bin/env perl

# count_fimo_hits_15oct2011.pl -- Erich Schwarz <emsch@caltech.edu>, 10/15/2011.  LEGACY VERSION, kept around for reproducibility of older work.
# Purpose: given a fimo.txt output, count hits per WBGene.  Note that this is specialized; a more general version would just get sequences.

use strict;
use warnings;

my $wbgene      = q{};
my %wbgene2hits = ();

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A \S+ \t (WBGene\d+) /xms ) { 
        $wbgene = $1;
        $wbgene2hits{$wbgene}++;
    }
}

foreach my $wbgene1 ( sort keys %wbgene2hits ) { 
    print "$wbgene1\t$wbgene2hits{$wbgene1}\n";
}

