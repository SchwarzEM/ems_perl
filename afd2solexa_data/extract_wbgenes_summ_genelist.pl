#!/usr/bin/perl

# extract_wbgenes_summ_genelist.pl -- Erich Schwarz <emsch@its.caltech.edu>, 2/22/2008.
# Purpose: just get the full gene list from a summary.

use strict;
use warnings;

my %seen_gene = (); 

while (my $input = <>) { 
    if ($input =~ /\A ( WBGene\d+\S+ ) \s /xms) { 
       $seen_gene{$1} = 1;
    }
}

foreach my $gene (sort keys %seen_gene) { 
    print "$gene\n";
}

