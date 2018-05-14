#!/usr/bin/env perl

# list_blast_hitnames.pl -- Erich Schwarz <emsch@its.caltech.edu>, 10/13/2009.
# Purpose: get a nonredundant, sorted list of Blast hit names from a NCBI BLAST report.

use strict;
use warnings;

my %seen = ();

while ( my $input = <> ) { 
    if ($input =~ /\AQuery=\s+(\S+)/xms ) { 
        my $name = $1;
        $seen{$name} = 1;
    }
}

foreach my $name (sort keys %seen) { 
    print "$name\n";
}

