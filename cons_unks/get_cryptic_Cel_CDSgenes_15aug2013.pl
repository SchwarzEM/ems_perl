#!/usr/bin/env perl

use strict;
use warnings;

# NCBI Tax. ID 6239  = caenorhabditis_elegans

my %seen = ();

while (my $input = <>) { 
    chomp $input;
    while ( $input =~ / G = ([^=]+) : T = ( 6239 ) \D /xmsg ) {
        my $cds_gene = $1;
        $seen{$cds_gene} = 1;
    }
}

my @genes = sort keys %seen;

foreach my $cds_genes (@genes) { 
    print "$cds_genes\n";
}

