#!/usr/bin/env perl

use strict;
use warnings;

my $data_ref;

while (my $input = <>) {
    if ( $input =~ /\A \S+ \t (\S+) /xms ) { 
        my $long_gene = $1;
        my $mod_gene  = $long_gene;
        if ( $mod_gene =~ /\A ([^|]+) \| /xms ) { 
            $mod_gene = $1;
        }
        $data_ref->{'mod_gene'}->{$mod_gene}->{'long_gene'}->{$long_gene} = 1;
    }
}

my @mod_genes = sort keys %{ $data_ref->{'mod_gene'} };
foreach my $mod_gene (@mod_genes) { 
    my @long_genes = sort keys %{ $data_ref->{'mod_gene'}->{$mod_gene}->{'long_gene'} };
    my $long_gene_count = @long_genes;
    if ( $long_gene_count >= 2 ) { 
        foreach my $long_gene (@long_genes) {
            print "$mod_gene\t$long_gene\n";
        }
    }
}


