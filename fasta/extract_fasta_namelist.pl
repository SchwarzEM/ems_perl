#!/usr/bin/env perl

# extract_fasta_namelist.pl -- Erich Schwarz <emsch@its.caltech.edu>, 9/30/2002; reworked 5/12/2008, and again on 11/2/2010.
# Purpose: convert ">text ..." to sorted list of "text".

use strict;
use warnings;

my %names = ();

while (my $input = <>) { 
    if ($input =~ /\A > (\S+) /xms) { 
        my $id = $1;
        if (exists $names{$id}) { 
            die "Redundant sequence name $id in header: $input\n";
        }
        $names{$id} = 1;
    }
}

foreach my $name (sort keys %names) { 
    print "$name\n";
}

