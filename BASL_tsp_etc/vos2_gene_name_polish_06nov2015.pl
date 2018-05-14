#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my %synonyms = (
    'Arabidopsis_thaliana.1' => 'At1g09710',
    'Arabidopsis_thaliana.2' => 'At1g58220_vos2',
    'Pyrus_x' => 'Pyrus_x_bretschneideri',
);

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A > (\S+) (.+) \z/xms ) { 
        my $gene = $1;
        my $text = $2;

        if ( exists $synonyms{$gene} ) {
            $gene = $synonyms{$gene};
            print ">$gene$text\n";
        }
        else { 
            print "$input\n";
        }
    }
    else {
        print "$input\n";
    }
}


