#!/usr/bin/env perl

# pseudomask_N_dna.pl -- Erich Schwarz, orig. 9/25/2003, update 11/10/2011, <emsch@caltech.edu>.
# Purpose: pseudomask 'n' or 'N' nucleotides as 'A'.

use strict;
use warnings;

while (my $input = <>) {
    chomp $input;
    if ( $input !~ /\A > /xms ) {
        $input =~ s/n/A/g;
        $input =~ s/N/A/g;
    }
    print "$input\n";
}

