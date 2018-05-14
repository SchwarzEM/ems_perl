#!/usr/bin/env perl

# filter_21oct2009_5.pl -- get a non-exon ncRNA subset (from pre-filtered "Transcript \"[^\"\s]+\" lines from huge GFF files).

use strict;
use warnings;

while (my $input = <>) { 
    chomp $input;
    if ( $input !~ / \A \S+ \t \S+ \t (exon) /xms ) { 
        print "$input\n";
    }
} 

