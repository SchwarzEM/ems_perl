#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $main_annots = $ARGV[0];
my $raw_q_values = $ARGV[1];

my %gene2q_value = ();

my $header = "Gene\tAWC_q_value\n";

open my $MAIN, '<', $main_annots;
while (my $input = <$MAIN>) {
    chomp $input;
    if ( $input =~ /\A (WBGene\S+) /xms ) {
        my $gene = $1;
        my $q_value = q{>=0.1};
        # The default value, which remains in place *unless* a specific annotation
        $gene2q_value{$gene} = $q_value;
    }
}
close $MAIN;

open my $QVALS, '<', $raw_q_values;
while (my $input = <$QVALS>) {
    if ( $input =~ /\A (WBGene\S+) \t (?: [^\t]* \t){5} (\S+) /xms ) { 
        my $gene    = $1;
        my $q_value = $2;
        $gene2q_value{$gene} = $q_value;
    }
}
close $QVALS;

my @genes = sort keys %gene2q_value;
foreach my $gene (@genes) {
    print $header if $header;
    $header = q{};
    my $q_value = $gene2q_value{$gene};
    print "$gene\t$q_value\n";
}


