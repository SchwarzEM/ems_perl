#!/usr/bin/env perl

# pick_strict_omcls.pl -- Erich Schwarz <emsch@its.caltech.edu>, 3/4/2010.
# Purpose: given gene-centric OrthoMCL file/stream, print only the strict orthology groups.

use strict;
use warnings;

my $gene_no   = 0;
my $taxon_no  = 0;

while (my $input = <>) { 
    chomp $input;
    $gene_no   = 0;
    $taxon_no  = 0;

    # Sample input line:
    # ORTHOMCL0(467 genes,6 taxa):      Bm1_00320(b_malayi) [etc.]

    if ( $input =~ / \A ORTHOMCL\d+ \( (\d+) [ ] genes, (\d+) [ ] taxa \): \s+ \S.+\S \s* \z /xms ) { 
        $gene_no   = $1;
        $taxon_no  = $2;
        if ( ( $gene_no == $taxon_no ) and ( $gene_no > 0 ) ) { 
            print "$input\n";
        }
    }

    # Enforce successful parsing.
    else { 
        die "Can't parse input line: $input\n";
    }
}

