#!/usr/bin/env perl

# swap_gtf_gene_tx_names.pl -- Erich Schwarz <ems394@cornell.edu>, 7/11/2013.
# Purpose: switch order of gene_id and transcript_id in GTF column 9; motivated by a desperate attempt to make a kludge script work.

use strict;
use warnings;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A ( (?:[^\t]* \t){8} ) ([^\t]+) \z /xms ) { 
        my $bulk_text   = $1;
        my $gene_tx_ids = $2;
        if (    ( $gene_tx_ids !~ /\A gene_id [^;]+       ; [ ] transcript_id [^;]+ ; \z/xms ) 
             and ( $gene_tx_ids !~ /\A transcript_id [^;]+ ; [ ] gene_id [^;]+      ; \z/xms ) ) {
            die "Can't parse gene/tx id: $input\n";
        }
        if ( $gene_tx_ids =~ /\A ([^;]+) ; [ ] ([^;]+) ; \z/xms ) {
            my $part1 = $1;
            my $part2 = $2;
            $gene_tx_ids = $part2 . q{; } . $part1 . q{;};
        }
        $input = $bulk_text . $gene_tx_ids;
    }
    print "$input\n";
}

