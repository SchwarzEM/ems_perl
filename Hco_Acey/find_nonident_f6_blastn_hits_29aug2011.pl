#!/usr/bin/env perl

# find_nonident_f6_blastn_hits_29aug2011.pl -- Erich Schwarz <emsch@caltech.edu>, 8/29/2011.
# Purpose: filter out self-hits from a format-6 BlastN report (useful for finding possible allelism).

use strict;
use warnings;

while (my $input = <>) { 
    chomp $input;
    my $seq1 = q{};
    my $seq2 = q{};
    if ( $input =~ /\A (\S+) \t (\S+) \t /xms ) { 
        $seq1 = $1;
        $seq2 = $2;
        if ( $seq1 ne $seq2 ) { 
            print "$input\n";
        }
    }
    else {
        warn "Misformatted line: $input\n";
    } 
}

