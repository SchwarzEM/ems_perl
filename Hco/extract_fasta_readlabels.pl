#!/usr/bin/env perl

# extract_fasta_readlabels.pl -- Erich Schwarz <emsch@caltech.edu>, 3/29/2012.
# Purpose: given many FASTA Illumina reads, get a list of their initial labels "ABCD1234:5" from them.

use strict;
use warnings;

my %observed = ();

while (my $input = <>) { 
    if ( $input =~ /\A > ([^:\s]+ [:] [^:\s]+) [:] /xms ) { 
        my $label = $1;
        $observed{$label} = 1;
    }
}

foreach my $obs_label (sort keys %observed) { 
    print "$obs_label\n";
}

