#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t/xms ) {
        my $seq = $1;
        $input =~ s/TSEBRA_/$seq./g;
        if ( $input =~ /TSEBRA/xms ) {
            die "Failed to fully correct input: $input\n";
        }
        print "$input\n";
    }
    else {
        die "Cannot parse input: $input\n";
    }
}

