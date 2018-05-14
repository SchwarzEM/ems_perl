#!/usr/bin/env perl

# check_MincV1A1_protnames.pl -- Erich Schwarz <emsch@its.caltech.edu>, 1/27/2010.
# Purpose: check the headers of MincV1A1 for non-redundant protein names and gene loci.

# Sample input:
# >gnl|MincDB|prot:Minc05158 length:799 contig:MiV1ctg127 region:1761-9463 strand:-

use strict;
use warnings;

my $protname = q{};

my %seen     = ();

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A > gnl\|MincDB\|prot: 
                       (\w+) 
                       \s+ length:\d+ \s+ 
                       contig:\w+ \s+ region:\d+\-\d+ \s+ strand:[+-]
                       \s* \z /xms ) { 
        $protname = $1;
        if ( $seen{$protname} ) { 
            die "Protein redundancy, in:  $input\n";
        }
        $seen{$protname} = 1;
    }
    elsif ( $input =~ / \A > /xms ) { 
        die "Unparseable header:  $input\n";
    }
}

print "\nAll header lines are free of protein redundancies!\n\n";

