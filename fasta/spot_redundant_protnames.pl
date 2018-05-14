#!/usr/bin/env perl

# spot_redundant_protnames.pl -- Erich Schwarz <emsch@its.caltech.edu>, 7/21/2009.
# Purpose: list all names from a FASTA file that occur 2 or more times.

use strict;
use warnings;

my %protnames = ();

while (my $input = <>) { 
    if ( $input =~ /\A > (\S+) /xms ) { 
        my $id = $1;
        $protnames{$id} += 1;
    }
    elsif ( $input =~ /\A > /xms ) { 
        die "Bizarre header line: $input\n";
    }
}

foreach my $prot (sort keys %protnames) { 
    if ($protnames{$prot} > 1 ) { 
        print "$prot occurs $protnames{$prot} times\n";
    }
}


