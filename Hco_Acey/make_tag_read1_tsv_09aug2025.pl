#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while ( my $input = <> ) {
    chomp $input;

    if ( $input =~ /\A (?: \/ [^\s\/]+ ){6} \/ ( [^\s\/]+) \/ ( [^\s\/]+) \/ /xms ) {
        my $tag_01 = $1;
        my $tag_02 = $2;
        my $tag = $tag_01. '_' . $tag_02;
        print "$tag\t$input\n";
    }

    else {
        die "Cannot parse input: $input\n";
    }
}
