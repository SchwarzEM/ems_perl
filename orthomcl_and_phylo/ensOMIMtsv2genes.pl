#!/usr/bin/env perl

# ensOMIMtsv2genes.pl -- Erich Schwarz <emsch@its.caltech.edu>, 2/28/2010.
# Purpose: get ordered genelist from ensembl_to_OMIM TSV made by EnsMart.

# Sample input lines:
# 
# Human gene:     Human protein:  OMIM g: OMIM disease:
# 
# ENSG00000218497 ENSP00000394086         
# ENSG00000148828 ENSP00000357501 609032  
# ENSG00000184895 ENSP00000372547 480000  306100

use strict;
use warnings;

my $gene         = q{};
my %morbid_genes = ();

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ / \A (ENSG\d+) \t ENSP\d+ \t \d+ \t \d+ /xms ) { 
        $gene = $1;
        $morbid_genes{$gene} = 1;
    }
}

foreach my $morbid_gene (sort keys %morbid_genes) { 
    print "$morbid_gene\n";
}

