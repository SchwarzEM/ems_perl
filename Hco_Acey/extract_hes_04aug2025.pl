#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $hes   = q{};
my $pairs = q{};

$hes   = $ARGV[0] if $ARGV[0];
$pairs = $ARGV[1] if $ARGV[1];

my %seen = ();

if ( (! $hes ) or (! $pairs ) ) {
    die "Format: extract_hes_04aug2025.pl [list of HES peptides] [gene-peptide pairs, TSV] > [subset of gene-HES peptide pairs, TSV]\n";
}

open my $HES, '<', $hes;
while ( my $peptide = <$HES> ) {
    chomp $peptide;
    if ( $peptide =~ /\A \S+ \z/xms ) {
        $seen{$peptide} = 1;
    }
    else {
        die "From HES list $hes, cannot parse: $peptide\n";
    }
}
close $HES;

open my $PAIRS, '<', $pairs;

while ( my $input = <$PAIRS> ) {
    chomp $input;
    if ( $input =~ /\A \S+ \t (\S+) \z/xms ) {
        my $peptide = $1;
        if ( exists $seen{$peptide} ) {
            print "$input\n";
        }
    }
    else {
        die "From gene-peptide pairs $pairs, cannot parse: $input\n";
    }
}
close $PAIRS;

