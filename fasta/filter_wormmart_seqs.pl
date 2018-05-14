#!/usr/bin/env perl

# filter_wormmart_seqs.pl -- Erich Schwarz <emsch@its.caltech.edu>, 10/26/2010.
# Purpose: censor unwanted text in WormMart-produced FASTAs, i.e., ">WBGene00000011|abc-1 / Sequence unavailable".

use strict;
use warnings;

my $print = 0;
my $line1 = q{};
my $line2 = q{};

while ( my $input = <> ) { 
    if ( $input =~ /\A > \S+ /xms ) {
        $line1 = $input;
        $line2 = q{};
        $print = 0;
    }
    else { 
        if ($print) { 
            print $input;
        }
        else { 
            if ( $input =~ / \A Sequence \s unavailable /xms ) {
                $line1 = q{};
                $line2 = q{};
                $print = 0;
            }
            else { 
                $print = 1;
                $line2 = $input;
                print $line1;
                print $line2;
            }
        }
    }
}

