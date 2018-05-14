#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;
    my @infields = split /\t/, $input;
    my @outfields = ();
    foreach my $infield (@infields) {
        if ( ( $infield =~ /\A \S+ \z/xms ) or ( $infield !~ /\A .+ \z/xms ) ) { 
            $infield = q{="} . $infield . q{"};
        }
        push @outfields, $infield;
    }
    my $output = join "\t", @outfields;
    print "$output\n";
}

      
