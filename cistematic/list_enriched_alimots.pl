#!/usr/bin/env perl

# list_enriched_alimots.pl -- Erich Schwarz <emsch@its.caltech.edu>, 11/2/2008.
# Purpose: get sorted nonredundant list of Ali motifs from consolidated enriched list.

use strict;
use warnings;

my %seen = ();

while (my $input = <>) { 
    if ($input =~ / (\S+_meme_\S+) /xms ) { 
        my $mot = $1;
        $seen{$mot} = 1;
    }
}

foreach my $motif (sort keys %seen) { 
    print "$motif\n";
}

