#!/usr/bin/env perl

# list_ps1010_genes_in_omcl.pl -- Erich Schwarz <emsch@its.caltech.edu>, 2/3/2010.
# Purpose: extract a nonredundant, ordered list of the Feb. 2010 PS1010 genes in an gene-centric OrthoMCL or slice thereof.

use strict;
use warnings;

my %gene_list    = ();
my $ps1010_gene = q{};

while (my $input = <>) { 
    while ( $input =~ / (ps1010rel4_g\d+) \(ps1010\) /gxms ) { 
        $ps1010_gene = $1;
        $gene_list{$ps1010_gene} = 1;
    }
}      

foreach my $gene (sort keys %gene_list) { 
    print "$gene\n";
}

