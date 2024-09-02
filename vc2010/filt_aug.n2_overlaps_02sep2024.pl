#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $genelist = q{};
my $overlaps = q{};

my %seen = ();

$genelist = $ARGV[0] if $ARGV[0];
$overlaps = $ARGV[1] if $ARGV[1];

if ( (! $genelist ) or (! $overlaps ) ) {
    die "Format: filt_aug.n2_overlaps_02sep2024.pl [gene list] [overlap table] > [gene list-filtered overlap table]\n";
}

open my $GENELIST, '<', $genelist;
while (my $gene = <$GENELIST> ) {
    chomp $gene;
    if ( $gene =~ /\A \S+ \z/xms ) {
        $seen{$gene} = 1;
    }
    else {
        die "From gene list $genelist, cannot parse: $gene\n";
    }
}
close $GENELIST;

open my $OVERLAPS, '<', $overlaps;
while (my $input = <$OVERLAPS> ) {
    chomp $input;
    if ( $input =~ /\A \S+ \t ([^\t]+) \z/xms ) {
        my $gene_text = $1;
        my $print = 0;
        my @genes = split '; ', $gene_text;
        foreach my $gene (@genes) {
            if ( exists $seen{$gene} ) {
                $print = 1;
            }
        }
        if ( ( $print ) or ( $input =~ /\AGene\t/xms ) ) {
            print "$input\n";
        }
    }
    else {
        die "From overlaps table $overlaps, cannot parse: $input\n";
    }
}
close $OVERLAPS;

