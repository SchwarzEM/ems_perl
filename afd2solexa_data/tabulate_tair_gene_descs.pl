#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $header = "Gene\tAnnotation\n";

my %seen = ();

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A (AT\S+) \s+ protein_coding \s+ (.*) /xms ) { 
        my $gene = $1;
        my $desc = $2;
        $desc =~ s/\t/; /g;
        $desc =~ s/\A\s+//;
        $desc =~ s/\s+\z//;
        $desc =~ s/; ;/;/g;
        $gene =~ s/\.d+\z//;
        if ( $seen{$gene} ) { 
            die "Redundant gene annotated in: $input\n";
        }
        print $header if $header;
        $header = q{};
        print "$gene\t$desc\n";
    }
}

        
