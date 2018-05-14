#!/usr/bin/env perl

# wmart2geneid_tab.pl -- Erich Schwarz <emsch@its.caltech.edu>, 12/2/2009.
# Purpose: convert typical WormMart GeneID .tsv to a more useful generic .tsv

# Some typical inputs:
# 
# WBGene00000010  aat-9   Y53H1C.1
# WBGene00000011  abc-1   
# WBGene00000141 
# WBGene00000263          F23H11.5

use strict;
use warnings;

my %gene_data = ();

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ / \A (WBGene\d+) \s+ (\S.+\S) \z /xms ) { 
        my $gene  = $1;
        my $other = $2;
        my $cds   = q{};
        my $cgc   = q{};
        my $id    = q{};
        if ( $other =~ / \A (\S+) \s+ (\S+) \z /xms ) { 
            $cgc = $1;
            $cds = $2;
            if ( ( $cgc =~ /\A [^\s\-]+ \- [^\s\-]+ \z/xms ) 
                 and ( $cds =~ /\A [^\s\.]+ \. [^\s\.]+ \z/xms ) ) { 
                $id = $gene . "\t" . $cds . "\t" . $cgc;
                $gene_data{$id} = 1;
            }
        }
        elsif ( $other =~ / \A (\S+) \z /xms ) {
            $cds = $1;
            # Don't end this regex with '\z', to allow isoforms:
            if ( $cds =~ /\A [^\s\-\.]+ \. [^\s\-\.]+ /xms ) {
                $id = "$gene\t$cds\t";
                $gene_data{$id} = 1;
            }
        }
        else { 
            warn "Couldn't parse: $input\n";
        }
    }
}

foreach my $gene_id (sort keys %gene_data ) { 
    print "$gene_id\n";
}

