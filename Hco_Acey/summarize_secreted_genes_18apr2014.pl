#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my %gene2annot = ();

my $header = "Gene\tOther Secreted Groups";

while (my $input = <>) { 
    chomp $input;
    if ( ( $input !~ /\A Gene /xms ) and ( $input =~ /\A (\S+) \t (.*) \z/xms ) ) { 
        my $gene   = $1;
        my $text = $2;
        $text =~ s/\A[\t]+//;
        $text =~ s/[\t]+\z//;
        $text =~ s/[\t]+/\t/g;
        my @annots = split "\t", $text;
        if (@annots) {
            @annots = grep { $_ =~ /\S/xms }  sort @annots;
            my $summary = join '; ', @annots;
            $gene2annot{$gene} = $summary;
        }
    }
}

my @genes = sort keys %gene2annot;

foreach my $gene (@genes) {
    print "$header\n" if $header;
    $header = q{};
    print "$gene\t$gene2annot{$gene}\n";
}


