#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A > ([^\|\s]+ \| [^\|\s]+ \| ([^\|\s]+) \/ (\d+ \- \d+) ) \s* \z/xms ) {
        my $header   = $1;
        my $protein  = $2;
        my $residues = $3;
        my $mod_header = '>' . $protein . q{_} . "$residues\t$header";
        print "$mod_header\n";
    }
    elsif ( $input =~ /\A > /xms ) {
        die "Can't parse header: $input\n";
    }
    else {
        print "$input\n";
    }
}


