#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my %seen = ();

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A \S+ \z/xms ) {
        if (! $seen{$input} ) {
            print "$input\n";
            $seen{$input} = 1;
        }
    }
    else {
        die "Cannot parse: $input\n";
    }
}
