#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $header = 1;

while (my $input = <>) {
    chomp $input;
    if ($header) {
        print "$input\n";
        $header = 0;
    }
    elsif ( $input !~ /\A (?: [^\t]* \t){14} [^\t]* \z/xms ) {
        die "Can't parse input line: $input\n";
    }
    else {
        my @data     = split /\t/, $input;
        my $awc_rpkm = $data[1];
        my $gpcr     = $data[8];
        if ( ( $awc_rpkm >= 0.1 ) and ( $gpcr ) and ( $gpcr eq '7TM_GPCRs' ) ) {
            print "$input\n";
        }
    }
}

