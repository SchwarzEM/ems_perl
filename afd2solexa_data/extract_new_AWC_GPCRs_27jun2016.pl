#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use List::Util qw(max);

my $header = 1;

while (my $input = <>) {
    chomp $input;
    if ($header) {
        print "$input\n";
        $header = 0;
    }
    elsif ( $input !~ /\A (?: [^\t]* \t){27} [^\t]* \z/xms ) {
        die "Can't parse input line: $input\n";
    }
    else {
        my @data       = split /\t/, $input;
        my @minTPMs    = ( $data[5], $data[8], $data[10],  $data[12],  $data[14],  $data[16] );
        my $max_minTPM = max(@minTPMs);
        if ( ( $max_minTPM >= 0.1 ) and ( $data[21] ) and ( $data[21] eq '7TM_GPCRs' ) ) {
            print "$input\n";
        }
    }
}

