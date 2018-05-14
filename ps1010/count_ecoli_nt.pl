#!/usr/bin/env perl

# count_ecoli_nt.pl -- Erich Schwarz <emsch@its.caltech.edu>, 1/4/2009.
# Purpose: get quick count of how many nt a list of E. coli contigs from PS1010 assembly has.
# N.B.: totally depends on idiosyncracies of how Ali labels the contigs.
#       Also, designed to get this number from a FASTA file.

use strict;
use warnings;

my $nt_count;
my $contig_count;
my $total_nt;

while (my $input = <>) { 
    if ( $input =~ / \A (?: > ){0,1} NODE_\d+_length_ (\d+) /xms ) { 
        $nt_count = $1;
        $total_nt += $nt_count;
        $contig_count++;
    }
}

print "Total: $total_nt nt in $contig_count contigs.\n";

