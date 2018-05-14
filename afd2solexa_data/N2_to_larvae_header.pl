#!/usr/bin/env perl

# N2_to_larvae_header.pl -- Erich Schwarz <emsch@its.caltech.edu>, 7/23/2011.
# Purpose: fix anachronistic "N2" in header lines of tables for LC paper.

use strict;
use warnings;

my $reading_first_line = 1;

while (my $input = <>) { 
    chomp $input;
    if ($reading_first_line) { 
        $input =~ s/N2/larvae/g;
        $reading_first_line = 0;
    }
    print "$input\n";
}


