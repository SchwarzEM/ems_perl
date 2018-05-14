#!/usr/bin/env perl

use strict;
use warnings;

my %synonyms = (
    'HCOI00305500' => 'Hco_Hc24',
    'HCOI01318800' => 'Hco_Hc40',
);

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A>/xms ) { 
        if ( $input =~ /\A > (HCOI0(\d+))\.t\d+ \b/xms ) { 
            my $cds   = $1;
            my $label = $2;
            my $gene = 'Hco_' . $label;
            if ( exists $synonyms{$cds} ) {
                $gene = $synonyms{$cds};
            }
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

