#!/usr/bin/env perl

# filter_zero_table_rows.pl -- Erich Schwarz <emsch@its.caltech.edu>, 9/23/2008.
# Purpose: hack to weed out zero-data rows in r^2 table.

use strict;
use warnings;

while (my $input = <>) { 
    chomp $input;
    if ( $input !~ /\A \S+ (\t 0.00)+ \z /xms ) { 
        print "$input\n";
    }
}

