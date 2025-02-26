#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while ( my $input = <> ) {
    chomp $input;
    if ( $input !~ /\A [>] /xms ) {
        print "$input\n";
    }
    # Sample input:
    # >WashU_chrI|g2440|g2440.t1|1|7281
    elsif ( $input =~ /\A [>] ( (\S+) \| \S+ \| (\S+) \| \d+ \| \d+ ) \s* \z/xms ) {
        my $annot = $1;
        my $seq   = $2;
        my $id    = $3;
        my $output = '>' . "$seq.$id  $annot";
        print "$output\n";
    }
    else {
        die "Cannot parse input: $input\n";
    }
}

