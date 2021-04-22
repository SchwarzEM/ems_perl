#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;
    # sample input line: s0183.g942.t2
    if ( $input =~ /\A (s\d+ \. g\d+) \. t\d+ \z/xms ) {
        my $gene = 'Acey_' . $1;
        print "$gene\n";
    }
    else {
        die "Cannot parse input: $input\n";
    }
}
