#!/usr/bin/env perl

use strict;
use warnings;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A>/xms ) { 
        if ( $input =~ /\A > (NECAME_(\d{5})) \s* /xms ) {
            my $cds   = $1;
            my $label = $2;
            my $gene  = 'Nam-ASPR-g' . $label;
            print ">$gene\n";
        }
        else {
            die "Can't parse header: $input\n";
       }
    }
    else { 
        print "$input\n";
    }
}

