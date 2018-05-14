#!/usr/bin/env perl

# filter_21oct2009_7.pl -- filter exons etc. and comments from mGene GFF3s; also select AUGUSTUS 5'/3'-UTRs.

use strict;
use warnings;

while (my $input = <>) { 
    chomp $input;
    if ( ( $input =~ / \A \# /xms ) 
             or 
         ( $input =~ / \A \S+ \t \S+ \t 
                       (five_prime_UTR|5'\-UTR|CDS|three_prime_UTR|3'\-UTR) \t /xms ) ) { 
        print "$input\n";
    }
} 

