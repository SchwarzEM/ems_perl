#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $header = "Gene\tC0_orthogroup\tAcey_C01_orthologs";

my %seen = ();

while (my $input = <> ){
    chomp $input;
    if ( $input =~ /\A (OG\d+) \t (\S[^\t]*\S) \z/xms ) {
        my $og_id   = $1;
        my $og_list = $2;
        my @genes   = split '; ', $og_list;
        foreach my $gene (@genes) {
            if ( exists $seen{$gene} ) {
                warn "Redundant entry for: $gene\n";
            }
            print "$header\n" if $header;
            $header = q{};
            print "$gene\t$og_id\t$og_list\n";
            $seen{$gene} = 1;
        }
    }
}

