#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t ([^\t]*) \z/xms ) { 
        my $gene  = $1;
        my $ofind = $2;
        my $print = 1;

        if ( $gene eq 'Gene' ) { 
            $print = 0;
        }
        elsif ( $ofind =~ /\A OG\d+ [(] \d+ [ ] genes [,] (\d+) [ ] taxa [)] [:] /xms ) {
            my $taxcount = $1;
            if ( $taxcount >= 2 ) {
                $print = 0;
            }
        }

        if ( $print ) {
            print "$input\n";
        }
    }
}

