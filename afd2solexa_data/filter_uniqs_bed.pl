#!/usr/bin/perl

# filter_uniqs_bed.pl -- Erich Schwarz <emsch@its.caltech.edu>, 2/18/2008.
# Purpose: filter out, and print, only the non-obvious stuff from uniqs.bed.

use strict;
use warnings;

my $chromosome  = 'I|II|III|IV|V|X';

while (my $input = <>) { 
    chomp $input;
    if ( $input !~ /\A chr\w+ \s+ \d+ \s+ \d+ /xms ) { 
        print "$input\n";
    }
    if ( $input =~ /\A chr(?:$chromosome) \s+ (\d+) \s+ (\d+) /xms ) { 
        my $i;
        my $j;
        ($i, $j) = ($1, $2);
        if ( ($j - $i) != 24 ) { 
            print "$input\n";
        }
    }
}

