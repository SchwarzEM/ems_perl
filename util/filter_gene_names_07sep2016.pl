#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my %seen = ();

while (my $input = <>) {
    chomp $input;
    if (! $seen{$input} ) {
        print "$input\n";
        $seen{$input} = 1;
    }
}

