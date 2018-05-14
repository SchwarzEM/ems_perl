#!/usr/bin/env perl

# dat2pmidlist.pl -- Erich Schwarz <emsch@its.caltech.edu>, 5/9/2008.
# Purpose: get ordered list of PMID numbers from an Integr8-derived *.dat file.

use strict;
use warnings;

my %pmids = ();

while (my $input = <>) { 
    if ($input =~ /RX\s+PubMed=(\d+)/xms) { 
        $pmids{$1} = 1;
    }
}

foreach my $ref (sort keys %pmids) { 
    print "$ref\n";
}

