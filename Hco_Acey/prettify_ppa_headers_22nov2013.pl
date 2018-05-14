#!/usr/bin/env perl

use strict;
use warnings;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A>/xms ) { 
        if ( $input =~ /\A > PPA(\d+) /xms ) { 
            my $gene = $1;
            $gene = 'Ppa_' . $gene;
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

