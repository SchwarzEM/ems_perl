#!/usr/bin/env perl

# censor_emptyheadered_fasta.pl -- Erich Schwarz <emsch@caltech.edu>, 9/30/2012.
# Purpose: filter out (from file or stream) parts of a FASTA file with bogus, empty 'headers' like this (e.g., from a SAM/BAM unmapped read output): ">     ".

use strict;
use warnings;

my $print = 0;

while (my $input = <>) { 
    if ( $input =~ /\A > /xms ) { 
        if ( $input =~ /\A > \S+ /xms ) {
            $print = 1;
        }
        else { 
            $print = 0;
        }
    }
    print $input if ( $print == 1 );
}

