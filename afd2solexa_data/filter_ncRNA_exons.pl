#!/usr/bin/env perl

# filter_ncRNA_exons.pl -- Erich Schwarz <emsch@its.caltech.edu>, 10/10/2008.
# Purpose: extract only ncRNA genes from a WormMart table of WS190 exon coords.

use strict;
use warnings;

my $genes_to_reject_ref;
my $exons_to_keep_ref;

my %terms_to_reject = ( 'coding'            => 1,
                        'Coding_pseudogene' => 1, 
                        'RNA_pseudogene'    => 1,
                        'pseudogene'        => 1,
                      );

while (my $input = <>) { 
    chomp $input;
    # Sample input:
    # WBGene00000002  aat-1   F27C8.1 coding  IV      9598963 9601672
    if ( $input =~ / \A (WBGene\d+) 
                     \t [^\t]* \t \S+ \t (\S+) 
                     \t [IVX]{1,2} \t \d+ \t \d+ \z 
                   /xms ) { 
        my $gene_id = $1;
        my $description = $2;
        if ( (! exists $genes_to_reject_ref->{$gene_id} ) 
             and ( $description =~ /\S/xms ) 
             and ( exists $terms_to_reject{$description} ) 
           ) { 
            $genes_to_reject_ref->{$gene_id} = 1;
            if ( exists $exons_to_keep_ref->{$gene_id} ) { 
                delete $exons_to_keep_ref->{$gene_id};
            }
        }
        if ( (! exists $genes_to_reject_ref->{$gene_id} ) 
              and ($description ne 'coding')              ) {
            $exons_to_keep_ref->{$gene_id}->{$input} = 1;
        }
    }
}

foreach my $gene (sort keys %{ $exons_to_keep_ref } ) { 
    foreach my $exon (sort keys %{ $exons_to_keep_ref->{$gene} } ) { 
        print "$exon\n"
    }
}

