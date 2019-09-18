#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \s+ \d+ \s+ \d+ \s+ ([ACGT]+) \s (\S.+\S) \s* \z/xms ) {
        my $name   = $1;
        my $seq    = $2;
        my $header = $3; 
        print ">$name  $header\n";
        print "$seq\n";
    }
    else {
        die "Cannot parse: $input\n";
    }
}

