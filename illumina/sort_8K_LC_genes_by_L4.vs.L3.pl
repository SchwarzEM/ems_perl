#!/usr/bin/env perl

use strict;
use warnings;

my @genes_8K;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A WBGene\d+\S+ \t (\S+) \t (\S+) \t /xms ) { 
        my $rpkm_l3 = $1;
        my $rpkm_l4 = $2;    
        if ( ( $rpkm_l3 > 0 ) or ( $rpkm_l4 > 0 ) ) { 
            push @genes_8K, $input;
        }
    }
}

@genes_8K = sort { get_ratio_l4_l3($b) <=> get_ratio_l4_l3($a) } @genes_8K;

foreach my $gene_8K (@genes_8K) { 
    print "$gene_8K\n";
}

sub get_ratio_l4_l3 { 
    my $_input = $_[0];
    if ( $_input !~ /\A WBGene\d+\S+ \t \S+ \t \S+ \t /xms ) {
        die "Can't parse $_input\n";
    }
    if ( $_input =~ /\A WBGene\d+\S+ \t (\S+) \t (\S+) \t /xms ) {
        my $_rpkm_l3 = $1;
        my $_rpkm_l4 = $2;
        if ( $_rpkm_l3 == 0 ) { 
            $_rpkm_l3 = 0.01;
        }
        if ( $_rpkm_l4 == 0 ) {
            $_rpkm_l4 = 0.01;
        }
        my $ratio_l4_l3 = ($_rpkm_l4 / $_rpkm_l3);
        return $ratio_l4_l3;
    }
}

