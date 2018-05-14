#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $prefix = q{};
my $DIGITS = q{};
my $list   = q{};
my $i      = 0;

$prefix = $ARGV[0] if $ARGV[0];
$DIGITS = $ARGV[1] if $ARGV[1];
$list   = $ARGV[2] if $ARGV[2];

if ( ( $prefix !~ /\A \S+ \z/xms ) or ( $DIGITS !~ /\A \d+ \z/xms ) or (! -r $list) ) {
    die "Format: map_genes_to_ncbi_ids.pl [NCBI locus prefix] [number of zero-padded digits] [gene name list] > [TSV of gene-to-locus]\n";
}

open my $LIST, '<', $list;
while (my $input = <$LIST>) {
    chomp $input;
    if ( $input =~ /\A \S+ \z/xms ) { 
        $i++;

        my $j = $i;
        $j = sprintf("%0${DIGITS}u", $j) or die "Can't zero-pad serial number $j\n";

        print "$input\t$prefix", '_', "$j\n";
    }
    else {
        die "In list $list, cannot parse: $input\n";
    }
}
close $LIST;
