#!/usr/bin/env perl

# maker_gene2tx_table.pl -- Erich Schwarz <emsch@caltech.edu>, 11/18/2012.
# Purpose: given a set of MAKER predictions for protein-coding gene transcripts, make a "gene \t transcript" table (e.g., for use by RSEM).

use strict;
use warnings;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A > /xms ) {
        if ( $input =~ /\A > ( (\S+) [-]mRNA[-]\d+) \b/xms ) { 
            my $tx   = $1;
            my $gene = $2;
            print "$gene\t$tx\n"; 
        }
        else { 
            die "Can't parse input: $input\n";
        }
    }
} 
