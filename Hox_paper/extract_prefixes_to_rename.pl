#!/usr/bin/perl

# extract_prefixes_to_rename.pl -- Erich Schwarz <emsch@its.caltech.edu>, 12/2/2007.
# Purpose: get self-ordered but nonrepeated list of prefixes to rename for final Kuntz et al. paper.

use strict;
use warnings;

my %seen = ();

while (my $input = <>) { 
    if ( $input =~ /\A ( (?:CB5161|PS1010) \S+ \.tfa ) /xms ) { 
        my $prefix = $1;
        if (! $seen{$prefix} ) {
            print "$prefix\n";
        }
        $seen{$prefix} = 1;
    }
}

