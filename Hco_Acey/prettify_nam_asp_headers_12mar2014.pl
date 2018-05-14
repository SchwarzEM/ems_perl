#!/usr/bin/env perl

use strict;
use warnings;

my %synonyms = ( 
    'NECAME_09334' => 'Nam-ASP-1',
    'NECAME_11281' => 'Nam-ASP-2',
    'NECAME_15644' => 'Nam-ASP-6',
    'NECAME_01333' => 'Nam-ASP-7',
);

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A>/xms ) { 
        if ( $input =~ /\A > (NECAME_(\d{5})) \s* /xms ) {
            my $cds   = $1;
            my $label = $2;
            my $gene  = 'Nam-ASP-g' . $label;
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

