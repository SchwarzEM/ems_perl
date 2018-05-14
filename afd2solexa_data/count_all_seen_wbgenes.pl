#!/usr/bin/env perl

# count_all_seen_wbgenes.pl -- Erich Schwarz <emsch@its.caltech.edu>, 9/24/2008.
# Purpose: get nonredundant count of all WBGenes in file stream.

use strict;
use warnings;

my %seen_wbgenes = ();

while (my $input = <>) { 
    if ($input =~ /\A .* (WBGene\d+\S*) /xms ) { 
        my $wbgene = $1;
        $seen_wbgenes{$wbgene} = 1;
    }
}

foreach my $gene (sort keys %seen_wbgenes) { 
    print "$gene\n";
}
