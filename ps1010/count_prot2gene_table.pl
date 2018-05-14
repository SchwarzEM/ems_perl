#!/usr/bin/env perl

# count_prot2gene_table.pl -- Erich Schwarz <emsch@its.caltech.edu>, 10/24/2009.
# Purpose: count proteins and genes in an output of prot2gene_table.pl.

use strict;
use warnings;

my %prots = ();
my %genes = ();

while (my $input = <>) { 
    chomp $input;
    if ( $input !~ /\A \S+ \t \S+ \z/xms ) { 
        die "Can't parse: $input!\n";
    }
    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) {
        my $prot = $1;
        my $gene = $2;
        $prots{$prot} = 1;
        $genes{$gene} = 1;
    }
}

my $prot_count = scalar(keys %prots);
my $gene_count = scalar(keys %genes);

print "Genes -- $gene_count; proteins -- $prot_count.\n";



