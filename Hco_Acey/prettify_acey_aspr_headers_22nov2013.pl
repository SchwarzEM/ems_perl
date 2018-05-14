#!/usr/bin/env perl

use strict;
use warnings;

# >ASP-s0001.g198
# >ASPR-s0002.g551

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A>/xms ) { 
        if ( $input =~ /\A > (Acey_2012\.08\.05_(\d+\.g\d+))\.t\d+ /xms ) { 
            my $cds   = $1;
            my $label = $2;
            my $gene  = 'ASPR-s' . $label;
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

