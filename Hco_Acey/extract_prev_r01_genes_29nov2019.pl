#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;
    if ( ( $input =~ /\t Intestinal_R01 \t/xms ) or ( $input =~ /\t Immunoregulated_R01 \t/xms ) ) {
        print "$input\n";
    }
}

