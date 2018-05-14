#!/usr/bin/env perl

# tbmk2geneid_tab.pl -- Erich Schwarz <emsch@its.caltech.edu>, 2/14/2010.
# Purpose: convert typical TableMaker GeneID .tsv to a more useful generic .tsv

# Typical input lines:
# 
# "WBGene00007062"        "Y59A8B.9"      "ebp-3" "Caenorhabditis elegans"
# "WBGene00007063"        "2L52.1"                "Caenorhabditis elegans"
# "WBGene00007064"        "2RSSE.1"               "Caenorhabditis elegans"
# "WBGene00007065"        "3R5.1" "pot-3" "Caenorhabditis elegans"

use strict;
use warnings;

my %gene_data = ();

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ / \A \" (WBGene\d+) \" \t (\S+) \t (\S*) \t /xms ) { 
        my $gene = $1;
        my $cds  = $2;
        my $cgc  = $3;
        $cds =~ s/\"//g;
        $cgc =~ s/\"//g;
        my $id    = q{};
        if ( $cgc =~ /\S/xms ) { 
            $id = $gene . "\t" . $cds . "\t" . $cgc;
            $gene_data{$id} = 1;
        }
        else { 
            $id = $gene . "\t" . $cds . "\t";
            $gene_data{$id} = 1;
        }
    }
    else { 
        die "Couldn't parse: $input\n";
    }
}

foreach my $gene_id (sort keys %gene_data ) { 
    print "$gene_id\n";
}

