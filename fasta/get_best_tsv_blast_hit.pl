#!/usr/bin/env perl

# get_best_tsv_blast_hit.pl -- Erich Schwarz <ems394@cornell.edu>, 5/8/2013.
# Purpose: given a tabular Blast output, report the best hit for each query (assumed to be the first one).

use strict;
use warnings;

my $data_ref;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A (\S+) \s+ (\S+) \s+/xms ) { 
        my $query  = $1;
        my $result = $2;
        if (! exists $data_ref->{'query'}->{$query} ) { 
            $data_ref->{'query'}->{$query} = $result;
        }
    }
}

my @queries = sort keys %{ $data_ref->{'query'} };
foreach my $query (@queries) {
    print "$data_ref->{'query'}->{$query}\n";
}

