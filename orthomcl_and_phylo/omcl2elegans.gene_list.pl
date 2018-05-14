#!/usr/bin/env perl

# omcl2elegans.gene_list.pl -- Erich Schwarz <emsch@its.caltech.edu>, 2/28/2010.
# Purpose: given a gene-centric OrthoMCL, get the elegans genes.

use strict;
use warnings;

my $gene       = q{};
my %worm_genes = ();

while (my $input = <>) { 
    chomp $input;
    while ( $input =~ / (WBGene\d+) \(elegans\) /gxms ) { 
        my $gene = $1;
        $worm_genes{$gene} = 1;
    }
}

foreach my $worm_gene ( sort keys %worm_genes ) { 
    print "$worm_gene\n";
}

