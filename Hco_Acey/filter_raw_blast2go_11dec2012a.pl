#!/usr/bin/env perl

use strict;
use warnings;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A (\S+) \- mRNA \- \d+ \s+ (GO:\d+) \b /xms ) { 
        my $gene    = $1;
        my $go_term = $2;
        print "$gene\t$go_term\n";
    }
}

