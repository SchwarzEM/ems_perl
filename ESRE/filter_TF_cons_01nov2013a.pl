#!/usr/bin/env perl

use strict;
use warnings;

my $header = q{};

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A WBGene\d+\S* \t [^\t]* \t ([^\t]*) \t ([^\t]*) \t /xms ) {
        my $flank      = $1;
        my $cons_flank = $2; 
        if (     ( ( $flank =~ /\d+/xms ) or ( $cons_flank =~ /\d+/xms ) ) 
             and ( ( $input =~ /homo_sapiens/xms ) and ( $input =~ /drosophila_melanogaster/xms ) ) ) {
            print "$input\n";
        }
    }
}

