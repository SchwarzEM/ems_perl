#!/usr/bin/env perl

# filter_21oct2009_6.pl -- get a non-exon ncRNA subset (from pre-filtered "Transcript \"[^\"\s]+\" lines from huge GFF files).

use strict;
use warnings;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ / \A \S+ \t \S+ \t (exon|five_prime_UTR|coding_exon|three_prime_UTR) \t /xms ) { 
        print "$input\n";
    }
} 

