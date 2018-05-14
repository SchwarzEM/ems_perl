#!/usr/bin/env perl

use strict;
use warnings;

my %seen = ();

while (my $input = <>) { 
    chomp $input;

    # Note that '(3)' denotes that a disease gene has actually been cloned;
    #     other numbers mean various forms of association which are not conclusive.
    #     Only the *first* name after (3)| is the gene name; everything after "," is a synonym or a disease name.

    if ( $input =~ /\A .+? \( 3 \) \s* (\S+) /xms ) { 
        my $gene = $1;
        $gene =~ s/[,]\z//;
        $seen{$gene} = 1;
    }
}

my @genes = grep { /\S/ } sort keys %seen;
foreach my $gene (@genes) {
    print "$gene\n";
}

