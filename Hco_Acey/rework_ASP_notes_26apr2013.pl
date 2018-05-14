#!/usr/bin/env perl

use strict;
use warnings;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A > ASP\-s(\d+\.g\d+) /xms ) { 
        my $serial_number = $1;
        my $output = 'Acey_2012.08.05_' . $serial_number;
        print "$output\t [aligned]\n";
    }
    elsif ( $input =~ /\A > (ASP\-\d) \s+ i\.e\.: \s+ ASP-s(\d+\.g\d+) /xms ) { 
        my $synonym       = $1;
        my $serial_number = $2;
        my $output = 'Acey_2012.08.05_' . $serial_number;
        print "$output\t [$synonym; aligned]\n";
    }
}


