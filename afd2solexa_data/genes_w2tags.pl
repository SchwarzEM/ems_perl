#!/usr/bin/perl

# solex_genes_w2tags.pl -- Erich Schwarz <emsch@its.caltech.edu>, 2/19/2008.
# Purpose: get various counts from a Solexa summary -- genes with 1 tag, 2+ tags; unique vs. ambiguous.

use strict;
use warnings;

my %genes_unitagged  = (); 
my %ambivs_unitagged = ();
my %genes_tagged     = ();
my %ambivs_tagged    = ();

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ / \A ( WBGene \d+ \S+ ) .+? \s+ 1 \z  /xms ) { 
        my $g1tag = $1;
        if ( $g1tag =~ /\|/xms ) {
            $ambivs_unitagged{$g1tag} = 1;
        }
        else { 
            $genes_unitagged{$g1tag} = 1;
        }
    }
    elsif ( $input =~ / \A ( WBGene \d+ \S+ ) .+? \s+ \d+ \z  /xms ) {
        my $gNtag = $1;
        if ( $gNtag =~ /\|/xms ) { 
            $ambivs_tagged{$gNtag} = 1;
        }
        else { 
            $genes_tagged{$gNtag} = 1;
        }
    }
}

my $tagged_count      = keys %genes_tagged;
my $ambiv_count       = keys %ambivs_tagged;
my $unitagged_count   = keys %genes_unitagged;
my $unitagambiv_count = keys %ambivs_unitagged;

print "Genes with >= 2 tags:             $tagged_count\n";
print "Ambiv. id. genes with >= 2 tags:  $ambiv_count\n";
print "Genes with 1 tag:                 $unitagged_count\n";
print "Ambiv. id. genes with 1 tag:      $unitagambiv_count\n";

