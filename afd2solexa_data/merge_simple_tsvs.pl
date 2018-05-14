#!/usr/bin/env perl

# merge_simple_tsvs.pl -- Erich Schwarz <emsch@its.caltech.edu>, 2/23/2010.
# Purpose: given a slew of individual one-off gene annotation TSVs, merge them into a single unified table; assumes "xxx." for each annot.

use strict;
use warnings;
use List::MoreUtils qw(uniq);

my $gene;
my $annot;
my $gene_data_ref;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A ([^\t]+) \t \" ([^\"\t]+) \. \" /xms ) { 
        $gene  = $1;
        $annot = $2;
        push @{ $gene_data_ref->{$gene} }, $annot;
    }
}

foreach my $gid (sort keys %{ $gene_data_ref } ) { 
    my @annots = uniq @{ $gene_data_ref->{$gid} };
    @annots = sort @annots;
    my $annot_text = join "; ", @annots;
    print "$gid\t\"$annot_text.\"\n";
}
