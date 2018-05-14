#!/usr/bin/env perl

# genlist2speccounts.pl -- Erich Schwarz <emsch@its.caltech.edu>, 4/25/2008.
# Purpose: ./genlist2speccounts.pl var_genelists/afd13_aser1_plml123_genes.txt uniqs_sums/afd1.uniqs.txt > genenames_w_counts.txt

use strict;
use warnings;

my %seen = ();
my %name = ();

while (my $input = <>) { 
    if ($input =~ /\A ( ( WBGene\d+ ) \| \S+ ) \s* /xms ) { 
        my $full_gene_name = q{};
        my $wbgene_id = q{};
        $full_gene_name = $1;
        $wbgene_id = $2;
        $name{$wbgene_id} = $full_gene_name;
        $seen{$wbgene_id} = 1;
    }
    if ($input =~ /\A (WBGene\d+) ( \s+ \d+ .* ) \z /xms) { 
        my $wbgene_id = q{};
        my $data      = q{};
        $wbgene_id = $1;
        $data = $2;
        if ( $seen{$wbgene_id} ) { 
            print $name{$wbgene_id}, $data;
        }
    }
}

