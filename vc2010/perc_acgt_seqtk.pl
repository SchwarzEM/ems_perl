#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \t (\S+) \t (\S+) \t (\S+)\z/xms ) {
        my $s = $1;  # $s for 'sum total'

        my $a = $2;
        my $c = $3;
        my $g = $4;
        my $t = $5;

        my $perc_a = ($a/$s);
        my $perc_c = ($c/$s);
        my $perc_g = ($g/$s);
        my $perc_t = ($t/$s);

        print "$perc_a\t$perc_c\t$perc_g\t$perc_t\n";
    }
    else {
        die "Cannot parse input: $input\n";
    }
}
