#!/usr/bin/env perl

# count_fimo_hits.pl -- Erich Schwarz <emsch@caltech.edu>, 10/15/2011.
# Purpose: given a fimo.txt output, count hits per sequence.

use strict;
use warnings;

my $seq      = q{};
my %seq2hits = ();

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A \S+ \t (\S+) /xms ) { 
        $seq = $1;
        $seq2hits{$seq}++;
    }
}

foreach my $seq1 ( sort keys %seq2hits ) { 
    print "$seq1\t$seq2hits{$seq1}\n";
}

