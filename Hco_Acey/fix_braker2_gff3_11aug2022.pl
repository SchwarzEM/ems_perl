#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t/xms ) {
        my $seq = $1;
        $input =~ s/file_1_file_1_j/$seq./g;
        $input =~ s/file_1_file_1_/$seq./g;
        $input =~ s/file_1_file_1_m1-/$seq./g;
        if ( $input =~ /file_1/xms ) {
            die "Failed to fully correct input: $input\n";
        }
        print "$input\n";
    }
    else {
        die "Cannot parse input: $input\n";
    }
}
