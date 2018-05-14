#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my %seen = ();

while (my $input = <>) {
    chomp $input;
    my $first_col_val = q{};
    if ( $input =~ /\A ([^\t]+) \t/xms ) { 
        $first_col_val = $1;
    }
    if (! exists $seen{$first_col_val} ) { 
        print "$input\n";
    }
    $seen{$first_col_val} = 1;
}


