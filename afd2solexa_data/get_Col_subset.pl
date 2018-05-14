#!/usr/bin/perl

use strict;
use warnings;

my %colosimo = ();

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A (WBGene\d+) \z/xms ) { 
        $colosimo{$1} = 1;
    }
    elsif ( ( $input =~ / (WBGene\d+) /xms ) 
            and ( $colosimo{$1} ) ) {
        print "$input\n";
    }
}

