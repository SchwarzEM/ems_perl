#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A (WBGene\d+\S+ \t \S+) \t (\S+) \t Housekeeping \z/xms ) { 
        my $data1 = $1;
        my $data2 = $2;
        print "$data1\t\t$data2\n";
    }
    else {
        print "$input\n";
    }
}


