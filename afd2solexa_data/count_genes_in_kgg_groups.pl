#!/usr/bin/env perl

# count_genes_in_kgg_groups.pl -- Erich Schwarz <emsch@its.caltech.edu>, 7/21/2011.
# Purpose: get count of genes in each numbered group of a *.kgg description of a k-clustering.

use strict;
use warnings;

my $group        = q{};
my %groups2count = ();
my @group_list   = ();
my %seen         = ();

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A WBGene\d+\S* \t (\d+) \s* \z /xms ) { 
        $group = $1;
        $groups2count{$group}++;
        if (! exists $seen{$group}) { 
            push @group_list, $group;
            $seen{$group} = 1;
        }
    }
}

foreach my $group_name (@group_list) { 
    print "$group_name\t$groups2count{$group_name}\n";
}


