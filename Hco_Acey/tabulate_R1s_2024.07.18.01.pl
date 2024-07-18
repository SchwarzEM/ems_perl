#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
        
while (my $input = <>) {
    chomp $input;
    # Sample input:
    # .../22_STAT6.KO_G14_trm.fastq.gz
    if ( $input =~ /\A \S+ \/ (\d+)_ (\S+?) _trm\.fastq\.gz \z/xms ) {
        my $number = $1;
        my $prefix = $2;
        print "$prefix" . 'rep_' . "$number\t$input\n";
    }
    else {
        die "Cannot parse input: $input\n";
    }
}
